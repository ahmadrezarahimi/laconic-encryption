package LE

import (
	"database/sql"
	_ "github.com/mattn/go-sqlite3"
	"github.com/tuneinsight/lattigo/v3/ring"
	"lattice-LE-go/matrix"
	"strconv"
	"sync"
)

/*
TreeHash is a 2to1 hash function
inputs: vectors of ring elements v1,v2 in normal(not NTT) form
output: H(v1,v2)= A0*(-G^{-1}(v1) + A1*(-G^{-1}(v2))
*/
func TreeHash(vec1, vec2 *matrix.Vector, le *LE) (vec3 *matrix.Vector) {
	n := len(vec1.Elements)
	vec3 = matrix.NewVector(n, le.R)
	u0 := vec1.GInvMNTT(le.R)
	u1 := vec2.GInvMNTT(le.R)
	matrix.Add(le.A0NTT.MulVecRight(u0, le.R), le.A1NTT.MulVecRight(u1, le.R), vec3, le.R)
	vec3 = vec3.InvNTT(le.R)
	return vec3
}

/*
Upd algorithm from the paper, here db refers to aux in the paper where we assume that the function has RAM access to db
The fucntion performs registration of the public key to the corresponding leave of the tree
inputs:

	le: the public parameters for the system
	db: a pointer to the database db
	layer: In this case the layer is the last layer of the binary tree (layer = \ell)
	row: In our case is ind
	vec: the public key

output: the function has no output.

	it updates the tree layer by layer
*/
func Upd(db *sql.DB, row uint64, layer int, vec *matrix.Vector, le *LE) {
	WriteToDB(db, layer, row, vec)
	TreeUpd(le, db, layer, row, vec)
	for i := layer - 1; i > 0; i-- {
		v := ReadFromDB(db, i, row>>(layer-i), le)
		TreeUpd(le, db, i, row>>(layer-i), v)
	}
}

/*
TreeUpd algorithm takes a node in the tree and updates it parent.
it first finds out if a node has a left or right sibling, and then it computes the Hash of the node and its sibling and writes it to the database
inputs:

	le: the public parameters for the system
	db: a pointer to the database
	layer, row: the exact location of the node
	vec: the node

outputs:

	the function has no output
*/
func TreeUpd(le *LE, db *sql.DB, layer int, row uint64, vec *matrix.Vector) {
	b := row & 1
	var sibling uint64
	if b == 0 {
		sibling = row + 1
	} else {
		sibling = row - 1
	}
	vecSib := le.YDefault
	//Checks if the sibling exists
	if CheckNode(db, layer, sibling) == true {
		vecSib = ReadFromDB(db, layer, sibling, le)
	}

	var vecHash *matrix.Vector
	if b == 0 {
		vecHash = TreeHash(vec, vecSib, le)
	} else {
		vecHash = TreeHash(vecSib, vec, le)
	}
	WriteToDB(db, layer-1, row>>1, vecHash)

}

// EncWithRandomness This function generates all the needed randomness within the encryption algorithm
func EncWithRandomness(le *LE, pp *matrix.Vector, id uint64, m *ring.Poly) ([]*matrix.Vector, []*matrix.Vector, *matrix.Vector, *ring.Poly) {

	r := make([]*matrix.Vector, le.Layers+1)
	for j := 0; j < le.Layers+1; j++ {
		r[j] = matrix.NewRandomVec(le.N, le.R, le.PRNG).NTT(le.R)
	}
	e := le.SamplerGaussian.ReadNew()
	e0 := make([]*matrix.Vector, le.Layers+1)
	e1 := make([]*matrix.Vector, le.Layers+1)
	for j := 0; j < le.Layers+1; j++ {
		e0[j] = matrix.NewNoiseVec(le.M, le.R, le.PRNG, le.Sigma, le.Bound).NTT(le.R)
		e1[j] = matrix.NewNoiseVec(le.M, le.R, le.PRNG, le.Sigma, le.Bound).NTT(le.R)
	}
	msg := le.R.NewPoly()

	q2 := le.Q / 2
	le.R.MulScalar(m, q2, msg)

	c0 := make([]*matrix.Vector, le.Layers)
	c1 := make([]*matrix.Vector, le.Layers)
	wg := sync.WaitGroup{}
	for i := 0; i < le.Layers; i++ {
		ithBitInv := 1 - (id>>(le.Layers-i-1))&1
		if ithBitInv == 1 {
			//wg.Add(1)
			go CreateEncBlockIf(le, c0, c1, r[i], r[i+1], e0[i], e0[i+1], i, &wg)
		} else {
			//wg.Add(1)
			go CreateEncBlockElse(le, c0, c1, r[i], r[i+1], e0[i], e0[i+1], i, &wg)
		}
	}
	//wg.Wait()
	c := le.BNTT.MulVecLeft(r[le.Layers], le.R)

	matrix.Add(c, e0[le.Layers], c, le.R)

	p := matrix.Mul(pp, r[0], le.R)
	le.R.InvNTT(p, p)

	le.R.Add(p, e, p)

	d := le.R.NewPoly()
	le.R.Add(p, msg, d)

	return c0, c1, c, d
}

/*
Enc function performs the encryption
In this case the randomness required for the encryption as well as the noise are both given to the function as inputs
We refer to the paper for the pseudocode for the function and full intuitive description of it.
Inputs:

	 	le: public parameters
		pp: the root of the tree, which itself is a part of public parameter
		id: the identity that we are encrypting to
		m: the message which is a ring element with binary coefficients.
		e0,e1 : noise vectors
		r0,r1 : the random ring element vectors for encryption

outputs:

	Outputs the ciphertext in 4 blocks.
	the output can be compressed into one TODO for future versions.
*/
func Enc(le *LE, pp *matrix.Vector, id uint64, m *ring.Poly, r []*matrix.Vector, e0 []*matrix.Vector, e1 []*matrix.Vector, e *ring.Poly) ([]*matrix.Vector, []*matrix.Vector, *matrix.Vector, *ring.Poly) {

	msg := le.R.NewPoly()

	q2 := le.Q / 2
	le.R.MulScalar(m, q2, msg)

	c0 := make([]*matrix.Vector, le.Layers)
	c1 := make([]*matrix.Vector, le.Layers)
	/*
		TODO at the moment, we are not using wait groups, so we are letting the parallel code to run in several processes.
		However, we should be careful to wait for the encryption to finish before trying to do the decryption.
		This can be resolved by using the WaitGroups, but they'll make the performance slightly worse.
	*/
	wg := sync.WaitGroup{}
	for i := 0; i < le.Layers; i++ {
		ithBitInv := 1 - (id>>(le.Layers-i-1))&1
		if ithBitInv == 1 {
			//wg.Add(1)
			go CreateEncBlockIf(le, c0, c1, r[i], r[i+1], e0[i], e0[i+1], i, &wg)
		} else {
			//wg.Add(1)
			go CreateEncBlockElse(le, c0, c1, r[i], r[i+1], e1[i], e1[i+1], i, &wg)
		}
	}
	//wg.Wait()
	c := le.BNTT.MulVecLeft(r[le.Layers], le.R)
	matrix.Add(c, e0[le.Layers], c, le.R)
	p := matrix.Mul(pp, r[0], le.R)
	le.R.InvNTT(p, p)
	le.R.Add(p, e, p)
	d := le.R.NewPoly()
	le.R.Add(p, msg, d)
	return c0, c1, c, d
}

/*
CreateEncBlockIf
This is one of the functions that performs one block of ciphertext.
We made it in a seperate function inorder to make the code run in parallel.
*/
func CreateEncBlockIf(le *LE, c0 []*matrix.Vector, c1 []*matrix.Vector, r0, r1, e0, e1 *matrix.Vector, i int, wg *sync.WaitGroup) {
	ct00 := matrix.NewVector(le.M, le.R)
	matrix.Add(le.A0NTT.MulVecLeft(r0, le.R), le.GNTT.MulVecLeft(r1, le.R), ct00, le.R)
	ct01 := le.A1NTT.MulVecLeft(r0, le.R)
	matrix.Add(ct00, e0, ct00, le.R)
	matrix.Add(ct01, e1, ct01, le.R)
	c0[i] = ct00
	c1[i] = ct01
	//wg.Done()
}

/*
CreateEncBlockElse
This is one of the functions that performs one block of ciphertext.
We made it in a seperate function inorder to make the code run in parallel.
*/
func CreateEncBlockElse(le *LE, c0 []*matrix.Vector, c1 []*matrix.Vector, r0, r1, e0, e1 *matrix.Vector, i int, wg *sync.WaitGroup) {
	ct00 := le.A0NTT.MulVecLeft(r0, le.R)
	ct01 := matrix.NewVector(le.M, le.R)
	matrix.Add(le.A1NTT.MulVecLeft(r0, le.R), le.GNTT.MulVecLeft(r1, le.R), ct01, le.R)
	matrix.Add(ct00, e0, ct00, le.R)
	matrix.Add(ct01, e1, ct01, le.R)
	c0[i] = ct00
	c1[i] = ct01
	//wg.Done()
}

/*
Dec is the decryption function.
We refer to the paper for the pseudocode for the function and full intuitive description of it.
This function runs a subroutine DecParallel for the parallel run.
Inputs:

	le: public parameter
	sk: secret key
	vec1,vec2 : Vectors of witnesses, these vectors can be obtained from the WGen (WitGen) algorithm
	c0,c1,c,d: ciphertext

Outputs:

	outputs the message.
*/
func Dec(le *LE, sk *matrix.Vector, vec1 []*matrix.Vector, vec2 []*matrix.Vector, c0 []*matrix.Vector, c1 []*matrix.Vector, c *matrix.Vector, d *ring.Poly) *ring.Poly {

	ctd := le.R.NewPoly()
	ctd = matrix.Mul(c, sk, le.R)
	le.R.InvNTT(ctd, ctd)
	m := le.R.NewPoly()
	le.R.Sub(d, ctd, m)

	ctd1 := make([]*ring.Poly, le.Layers)
	// TODO we are using WaitGroups for the parallel decryption as the order matters.
	for i := 0; i < le.Layers; i++ {
		ctd1[i] = le.R.NewPoly()
	}
	wg := sync.WaitGroup{}
	for i := 0; i < le.Layers; i++ {

		wg.Add(1)
		go DecParallel(le, c0[i], c1[i], vec1[i], vec2[i], ctd1, i, &wg)
	}
	wg.Wait()
	for i := 0; i < le.Layers; i++ {
		le.R.Sub(m, ctd1[i], m)
	}
	return m
}

func DecParallel(le *LE, c0, c1, u0, u1 *matrix.Vector, ctd1 []*ring.Poly, i int, wg *sync.WaitGroup) {
	p := le.R.NewPoly()

	le.R.Add(matrix.Mul(c0, u0, le.R), matrix.Mul(c1, u1, le.R), p)
	le.R.InvNTT(p, p)
	ctd1[i] = p
	wg.Done()
}

/*
WitGen This is the function that generates the witnesses required for decryption
The function is referred to WGen in the paper.
The function runs two subroutines WitGenParLeft and WitGenParRight for parallelization purposes
Inputs:

	le: public parameters'
	db: a pointer to the database
	id: the location of the database that we want to generate the witness for.

Outputs:

	a co-path from id to the root of the tree
	a co-path is a path where all the siblings of each node in the path-node are also included in it.
*/
func WitGen(db *sql.DB, le *LE, id uint64) ([]*matrix.Vector, []*matrix.Vector) {
	vecLeft := make([]*matrix.Vector, le.Layers)
	vecRight := make([]*matrix.Vector, le.Layers)
	//TODO we are not using WaitGroups here, as there is no need for the witnesses to be generated in a specific order. However all the processes must be finished before performing the decryption.
	wg := sync.WaitGroup{}
	for i := le.Layers; i > 0; i-- {
		idInd := id >> (le.Layers - i)
		b := (id >> (le.Layers - i)) & 1
		if b == 0 {
			//wg.Add(1)
			go WitGenParLeft(le, db, i, idInd, vecLeft, vecRight, &wg)
		} else {
			//wg.Add(1)
			go WitGenParRight(le, db, i, idInd, vecLeft, vecRight, &wg)
		}
	}
	//wg.Wait()
	return vecLeft, vecRight
}
func WitGenParLeft(le *LE, db *sql.DB, i int, ind uint64, vecLeft, vecRight []*matrix.Vector, wg *sync.WaitGroup) {
	vecLeft[i-1] = ReadFromDB(db, i, ind, le).GInvMNTT(le.R)
	vecRight[i-1] = ReadFromDB(db, i, ind+1, le).GInvMNTT(le.R)
	//wg.Done()
}

func WitGenParRight(le *LE, db *sql.DB, i int, ind uint64, vecLeft, vecRight []*matrix.Vector, wg *sync.WaitGroup) {
	vecLeft[i-1] = ReadFromDB(db, i, ind-1, le).GInvMNTT(le.R)
	vecRight[i-1] = ReadFromDB(db, i, ind, le).GInvMNTT(le.R)
	//wg.Done()
}

//TODO move these functions to a utils pakcage

/*
CheckNode This function checks if a node is already filled in the tree or not.
If it does not find any node, it returns the default value from the public parameters
TODO complete the comments
*/
func CheckNode(db *sql.DB, layer int, rowid uint64) bool {
	query := "SELECT y_def FROM tree_" + strconv.Itoa(layer) + " WHERE rowid = ?"
	row := db.QueryRow(query, rowid)
	var exists bool
	row.Scan(&exists)
	return exists
}

/*
ReadFromDB This function takes a location of a node in the tree as input and outputs the corresponding node
TODO complete the comments
*/
func ReadFromDB(db *sql.DB, layer int, row uint64, le *LE) (vec *matrix.Vector) {

	if CheckNode(db, layer, row) == false {
		return le.YDefault
	}
	query1 := "SELECT p1 FROM tree_" + strconv.Itoa(layer) + " WHERE rowid = ?"
	query2 := "SELECT p2 FROM tree_" + strconv.Itoa(layer) + " WHERE rowid = ?"
	query3 := "SELECT p3 FROM tree_" + strconv.Itoa(layer) + " WHERE rowid = ?"
	query4 := "SELECT p4 FROM tree_" + strconv.Itoa(layer) + " WHERE rowid = ?"
	query5 := "SELECT y_def FROM tree_" + strconv.Itoa(layer) + " WHERE rowid = ?"

	row1 := db.QueryRow(query1, row)
	row2 := db.QueryRow(query2, row)
	row3 := db.QueryRow(query3, row)
	row4 := db.QueryRow(query4, row)
	row5 := db.QueryRow(query5, row)

	var p1B, p2B, p3B, p4B []byte
	var yDef bool
	row1.Scan(&p1B)
	row2.Scan(&p2B)
	row3.Scan(&p3B)
	row4.Scan(&p4B)
	row5.Scan(&yDef)

	vecB := make([][]byte, 4)
	vecB[0] = p1B
	vecB[1] = p2B
	vecB[2] = p3B
	vecB[3] = p4B

	vec = matrix.NewVector(4, le.R)
	vec.Decode(vecB)
	return vec

}

/*
WriteToDB this function takes a vertex as input and encodes it to bytes arrays and writes it to the corresponding location of the database
TODO complete the comment
*/
func WriteToDB(db *sql.DB, layer int, row uint64, vec *matrix.Vector) {
	vecB := vec.Encode()

	var query1 string
	if CheckNode(db, layer, row) == false {
		query1 = "INSERT INTO tree_" + strconv.Itoa(layer) + " (rowid, p1) VALUES (?,?)"
		db.Exec(query1, row, vecB[0])
	} else {
		query1 = "UPDATE tree_" + strconv.Itoa(layer) + " SET p1 = ? WHERE rowid=?"
		db.Exec(query1, vecB[0], row)
	}
	query2 := "UPDATE tree_" + strconv.Itoa(layer) + " SET p2 = ? WHERE rowid=?"
	query3 := "UPDATE tree_" + strconv.Itoa(layer) + " SET p3 = ? WHERE rowid=?"
	query4 := "UPDATE tree_" + strconv.Itoa(layer) + " SET p4 = ? WHERE rowid=?"
	query5 := "UPDATE tree_" + strconv.Itoa(layer) + " SET y_def = ? WHERE rowid=?"
	db.Exec(query2, vecB[1], row)
	db.Exec(query3, vecB[2], row)
	db.Exec(query4, vecB[3], row)
	db.Exec(query5, true, row)
}
