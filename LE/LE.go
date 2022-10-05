package LE

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
	"lattice-LE-go/matrix"
	"math"
)

type LE struct {
	Q     uint64
	N     int
	M     int
	M2    int
	D     int
	QBits int

	Layers int

	R *ring.Ring

	Sigma float64
	Bound int

	A0 *matrix.Matrix
	A1 *matrix.Matrix
	B  *matrix.Matrix

	YDefault *matrix.Vector

	A0NTT *matrix.Matrix
	A1NTT *matrix.Matrix
	BNTT  *matrix.Matrix

	A0t *matrix.Matrix
	A1t *matrix.Matrix
	Bt  *matrix.Matrix

	A0tNTT *matrix.Matrix
	A1tNTT *matrix.Matrix
	BtNTT  *matrix.Matrix

	G    *matrix.Matrix
	GNTT *matrix.Matrix

	PRNG            utils.PRNG
	Sampler         *ring.UniformSampler
	SamplerGaussian *ring.GaussianSampler
}

func Setup(Q uint64, qBits int, D int, N int) (LEParams *LE) {
	LEParams = new(LE)
	LEParams.Layers = 50
	LEParams.Sigma = math.Pow(2, 30)
	LEParams.Bound = math.MaxInt32
	LEParams.QBits = qBits
	LEParams.N = N
	LEParams.M = N * qBits
	LEParams.Q = Q
	LEParams.M2 = 512
	r, _ := ring.NewRing(D, []uint64{Q})
	LEParams.R = r
	prng, _ := utils.NewPRNG()
	LEParams.PRNG = prng
	LEParams.Sampler = ring.NewUniformSampler(prng, r)
	LEParams.SamplerGaussian = ring.NewGaussianSampler(prng, r, LEParams.Sigma, LEParams.Bound)
	LEParams.A0 = matrix.NewRandomMatrix(N, N*qBits, r, prng)
	LEParams.A1 = matrix.NewRandomMatrix(N, N*qBits, r, prng)
	LEParams.B = matrix.NewRandomMatrix(N, LEParams.M2, r, prng)

	LEParams.YDefault = matrix.NewRandomVec(N, r, prng)
	for i := 0; i < N; i++ {
		LEParams.YDefault.Elements[i].Coeffs[0][0] = 444
	}

	LEParams.A0NTT = LEParams.A0.NTT(N, N*qBits, r)
	LEParams.A1NTT = LEParams.A1.NTT(N, N*qBits, r)
	LEParams.BNTT = LEParams.B.NTT(N, N*qBits, r)

	LEParams.A0t = LEParams.A0.Transpose()
	LEParams.A1t = LEParams.A1.Transpose()
	LEParams.Bt = LEParams.B.Transpose()

	LEParams.A0tNTT = LEParams.A0NTT.Transpose()
	LEParams.A1tNTT = LEParams.A1NTT.Transpose()
	LEParams.BtNTT = LEParams.BNTT.Transpose()

	LEParams.G = matrix.NewMatrix(LEParams.N, LEParams.M, r)
	for i := 0; i < LEParams.N; i++ {
		for j := 0; j < qBits; j++ {
			LEParams.G.Elements[i][i*qBits+j].Coeffs[0][0] = uint64(math.Pow(2, float64(j)))
		}
	}
	LEParams.GNTT = LEParams.G.NTT(N, N*qBits, r)

	return
}
