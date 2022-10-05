package matrix

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
	"math"
	"math/rand"
	"time"
)

type Vector struct {
	Elements []*ring.Poly
}

/*
NewVector generates a new vector and fills it with empty polynomials of the ring.
the empty polynomial means that all coefficients are set to 0
*/
func NewVector(n int, r *ring.Ring) (vec *Vector) {
	vec = new(Vector)
	vec.Elements = make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		vec.Elements[i] = r.NewPoly()
	}
	return

}

/*
NewNoiseVec generates a new noise vector with the given parameters.
*/
func NewNoiseVec(n int, r *ring.Ring, prng utils.PRNG, sigma float64, bound int) (vec *Vector) {
	vec = new(Vector)
	vec.Elements = make([]*ring.Poly, n)
	s := ring.NewGaussianSampler(prng, r, sigma, bound)
	for i := 0; i < n; i++ {
		vec.Elements[i] = s.ReadNew()
	}
	return
}

/*
NewRandomVec generates a new random vector with the given parameters.
*/
func NewRandomVec(n int, r *ring.Ring, prng utils.PRNG) (vec *Vector) {
	vec = new(Vector)
	vec.Elements = make([]*ring.Poly, n)
	s := ring.NewUniformSampler(prng, r)
	for i := 0; i < n; i++ {
		vec.Elements[i] = s.ReadNew()
	}
	return
}

/*
NewRandomVecBinary generates a new random vector with the given parameters.
the coefficients of the polynomials in the vector are from {0,1}
*/
func NewRandomVecBinary(n int, r *ring.Ring) (vec *Vector) {
	vec = new(Vector)
	vec.Elements = make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		p := NewRandomPolyBinary(r)
		vec.Elements[i] = p
	}
	return
}

/*
NewRandomPolyBinary generates a new polynomial with random coefficients from {0,1}
TODO this function is not using any secure random number generator
*/
func NewRandomPolyBinary(r *ring.Ring) (poly *ring.Poly) {
	p := r.NewPoly()
	pCoeffs0 := make([][]uint64, 1)
	pCoeffs0[0] = make([]uint64, r.N)

	rand.Seed(time.Now().UnixNano())
	for j := 0; j < r.N; j++ {
		pCoeffs0[0][j] = uint64(rand.Intn(2))
	}
	p.SetCoefficients(pCoeffs0)
	return p
}

/*
NTT computes the NTT representation of each polynomial in the vector and outputs a new vector in NTT form
*/
func (vec1 *Vector) NTT(r *ring.Ring) (vec2 *Vector) {
	n := len(vec1.Elements)
	vec2 = NewVector(n, r)
	for i := 0; i < n; i++ {
		r.NTT(vec1.Elements[i], vec2.Elements[i])
	}
	return
}

/*
InvNTT computes the inverse NTT of each ring element of the vector
*/
func (vec1 *Vector) InvNTT(r *ring.Ring) (vec2 *Vector) {
	n := len(vec1.Elements)
	vec2 = NewVector(n, r)
	for i := 0; i < n; i++ {
		r.InvNTT(vec1.Elements[i], vec2.Elements[i])
	}
	return
}

/*
Add adds two vector, elementwise
*/
func Add(vec1, vec2, vec3 *Vector, r *ring.Ring) {
	n := len(vec1.Elements)
	for i := 0; i < n; i++ {
		r.Add(vec1.Elements[i], vec2.Elements[i], vec3.Elements[i])
	}
	return
}

/*
Mul computes the inner product of two vectors
*/
func Mul(vec1, vec2 *Vector, r *ring.Ring) (p *ring.Poly) {
	n := len(vec1.Elements)
	p = r.NewPoly()
	for i := 0; i < n; i++ {

		//fmt.Println(p.Coeffs[0][0])
		p1 := r.NewPoly()
		r.MulCoeffs(vec1.Elements[i], vec2.Elements[i], p1)
		r.Add(p, p1, p)
	}
	return p
}

/*
GMulRight Computes Gv = w
*/
func GMulRight(n int, vec1 *Vector, r *ring.Ring) (vec2 *Vector) {
	vec2 = NewVector(n, r)
	for k := 0; k < n; k++ {
		for i := 0; i < 58; i++ {
			p := r.NewPoly()
			pow2 := uint64(math.Pow(2, float64(i)))
			for j := range p.Coeffs[0] {
				p.Coeffs[0][j] = vec1.Elements[i].Coeffs[0][j] * pow2
			}
			r.Add(vec2.Elements[k], p, vec2.Elements[k])
		}
	}
	return
}

/*
GMulLeft Computes v^{T}G  =w^{T}
*/
func GMulLeft(n int, vec1 *Vector, r *ring.Ring) (vec2 *Vector) {
	vec2 = NewVector(58*n, r)
	for i := 0; i < n; i++ {
		for j := 0; j < 58; j++ {
			for k := 0; k < r.N; k++ {
				vec2.Elements[i*58+j].Coeffs[0][k] = ring.BRed(vec1.Elements[i].Coeffs[0][k],
					uint64(math.Pow(2, float64(j))), r.Modulus[0], r.BredParams[0])
			}
		}
	}
	return
}

/*
Sub subtracts two vectors
*/
func Sub(vec1, vec2, vec3 *Vector, r *ring.Ring) {
	n := len(vec1.Elements)
	for i := 0; i < n; i++ {
		r.Sub(vec1.Elements[i], vec2.Elements[i], vec3.Elements[i])
	}
	return
}

/*
GInv computes G^{-1}(x) = y on a given vector x
*/
func (vec1 *Vector) GInv(r *ring.Ring) (vec2 *Vector) {
	n := len(vec1.Elements)
	m := n * 58
	vec2 = NewVector(m, r)
	for i := 0; i < n; i++ {
		for j := 0; j < r.N; j++ {
			binCoeffs := CoeffToBin(vec1.Elements[i].Coeffs[0][j])
			for k := 0; k < 58; k++ {
				vec2.Elements[i*58+k].Coeffs[0][j] = binCoeffs[k]
			}
		}

	}
	return
}

/*
GInvMNTT computes -G^{-1}(x) = y on a given vector x, then outputs it in NTT form
*/
func (vec1 *Vector) GInvMNTT(r *ring.Ring) (vec2 *Vector) {
	n := len(vec1.Elements)
	m := n * 58
	vec2 = NewVector(m, r)
	for i := 0; i < n; i++ {
		for j := 0; j < r.N; j++ {
			binCoeffs := CoeffToBin(vec1.Elements[i].Coeffs[0][j])
			for k := 0; k < 58; k++ {
				if binCoeffs[k] == 0 {
					vec2.Elements[i*58+k].Coeffs[0][j] = 0
				} else {
					vec2.Elements[i*58+k].Coeffs[0][j] = r.Modulus[0] - 1
				}
			}
		}
	}
	for i := 0; i < m; i++ {
		r.NTT(vec2.Elements[i], vec2.Elements[i])
	}
	return vec2
}

/*
CoeffToBin takes a uint64 as input and outputs the 58 bit binary representation of it.
TODO must be written for general QBit lengths
*/
func CoeffToBin(b uint64) []uint64 {
	out := make([]uint64, 58)
	for i := 0; i < 58; i++ {
		out[i] = b % 2
		b = b / 2
	}
	return out
}

/*
Encode encodes a vector to byte array
*/
func (vec *Vector) Encode() [][]byte {
	n := len(vec.Elements)
	vecB := make([][]byte, n)
	for i := 0; i < n; i++ {
		vecB[i], _ = vec.Elements[i].MarshalBinary()
	}
	return vecB
}

/*
Decode decodes a vector from its byte array representation
*/
func (vec *Vector) Decode(vecB [][]byte) {
	n := len(vec.Elements)
	for i := 0; i < n; i++ {
		_ = vec.Elements[i].UnmarshalBinary(vecB[i])
	}
	return
}
