package main

import (
	"database/sql"
	"fmt"
	_ "github.com/mattn/go-sqlite3"
	"github.com/tuneinsight/lattigo/v3/ring"
	"lattice-LE-go/LE"
	"lattice-LE-go/matrix"
	"math"
	"math/rand"
	"strconv"
	"time"
)

func main() {

	start := time.Now()
	le := LE.Setup(180143985094819841, 58, 1024, 4)
	end := time.Now()
	fmt.Println("Setting up the system takes: ", end.Sub(start))
	fmt.Println("Running system with 1000 users and d = 256")
	LaconicEncTEST(le, "trees/tree256.db", 1000)

	fmt.Println("___________________________________________________")
	start = time.Now()
	le2 := LE.Setup(180143985094819841, 58, 512, 4)
	end = time.Now()
	fmt.Println("Setting up the system takes: ", end.Sub(start))
	fmt.Println("Running system with 1000 users and d = 512")
	LaconicEncTEST(le2, "trees/tree512.db", 1000)
	fmt.Println("___________________________________________________")
	start = time.Now()
	le3 := LE.Setup(180143985094819841, 58, 1024, 4)
	end = time.Now()
	fmt.Println("Setting up the system takes: ", end.Sub(start))
	fmt.Println("Running system with 1000 users and d = 1024")
	LaconicEncTEST(le3, "trees/tree1024.db", 1000)

}

func LaconicEncTEST(le *LE.LE, dbPath string, registerTests int) {

	db, _ := sql.Open("sqlite3", dbPath)
	for i := 0; i <= 50; i++ {
		query := "CREATE TABLE tree_" + strconv.Itoa(i) + " (p1 BLOB, p2 BLOB, P3 BLOB, p4 BLOB, y_def BOOLEAN)"
		db.Exec(query)
	}

	rand.Seed(time.Now().UnixNano())
	idSet := make([]uint64, registerTests)
	mask := uint64(math.Pow(2, 50) - 1)
	for i := 0; i < registerTests; i++ {
		idSet[i] = ring.RandUniform(le.PRNG, mask, mask)
	}
	targetTests := registerTests
	pkTests := make([]*matrix.Vector, registerTests)
	skTests := make([]*matrix.Vector, registerTests)
	for i := 0; i < registerTests; i++ {
		pkTests[i], skTests[i] = le.KeyGen()

	}

	start := time.Now()
	for i := 0; i < registerTests; i++ {
		LE.Upd(db, idSet[i], 50, pkTests[i], le)
	}
	end := time.Now().Sub(start)
	fmt.Println("Creating the tree takes ", end)

	pp := LE.ReadFromDB(db, 0, 0, le).NTT(le.R)
	msg := matrix.NewRandomPolyBinary(le.R)
	//db.Close()
	encTime := make([]time.Duration, targetTests)
	witGenTime := make([]time.Duration, targetTests)
	decTime := make([]time.Duration, targetTests)

	for i := 0; i < targetTests; i++ {
		r := make([]*matrix.Vector, le.Layers+1)
		for j := 0; j < le.Layers+1; j++ {
			r[j] = matrix.NewRandomVec(le.N, le.R, le.PRNG).NTT(le.R)
		}
		e := le.SamplerGaussian.ReadNew()
		e0 := make([]*matrix.Vector, le.Layers+1)
		e1 := make([]*matrix.Vector, le.Layers+1)
		for j := 0; j < le.Layers+1; j++ {
			if j == le.Layers {
				e0[j] = matrix.NewNoiseVec(le.M2, le.R, le.PRNG, le.Sigma, le.Bound).NTT(le.R)
			} else {
				e0[j] = matrix.NewNoiseVec(le.M, le.R, le.PRNG, le.Sigma, le.Bound).NTT(le.R)
			}
			e1[j] = matrix.NewNoiseVec(le.M, le.R, le.PRNG, le.Sigma, le.Bound).NTT(le.R)
		}
		start = time.Now()
		c0, c1, c, d := LE.Enc(le, pp, idSet[i], msg, r, e0, e1, e)
		end = time.Now().Sub(start)
		encTime[i] = end

		time.Sleep(250 * time.Millisecond)
		start = time.Now()
		vec1, vec2 := LE.WitGen(db, le, idSet[i])
		end = time.Now().Sub(start)
		witGenTime[i] = end
		time.Sleep(250 * time.Millisecond)
		start = time.Now()
		m := LE.Dec(le, skTests[i], vec1, vec2, c0, c1, c, d)
		end = time.Now().Sub(start)
		decTime[i] = end
		if CorrecnessCheck(m, msg, le) == false {
			fmt.Println("Correctness Failed")
		}
	}
	encTimeAvg := time.Duration(0)
	decTimeAvg := time.Duration(0)
	witGenTimeAvg := time.Duration(0)
	for i := 0; i < targetTests; i++ {
		encTimeAvg += encTime[i]
		decTimeAvg += decTime[i]
		witGenTimeAvg += witGenTime[i]
	}
	fmt.Println("Encryption takes ", encTimeAvg/time.Duration(targetTests))
	fmt.Println("Getting Witness takes ", witGenTimeAvg/time.Duration(targetTests))
	fmt.Println("Decryption takes ", decTimeAvg/time.Duration(targetTests))

}

func CorrecnessCheck(p1 *ring.Poly, p2 *ring.Poly, le *LE.LE) bool {
	q14 := le.Q / 4
	q34 := (le.Q / 4) * 3
	p := le.R.NewPoly()
	for i := 0; i < le.R.N; i++ {
		if p1.Coeffs[0][i] < q14 || p1.Coeffs[0][i] > q34 {
			p.Coeffs[0][i] = 0
		} else {
			p.Coeffs[0][i] = 1
		}
	}
	b := p.Equals(p2)
	return b
}
