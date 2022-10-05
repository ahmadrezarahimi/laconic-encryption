package matrix

import (
	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/utils"
)

type Matrix struct {
	Elements [][]*ring.Poly // Dimension-2 slice of coefficients (re-slice of Buff)
}

func NewMatrix(n, m int, r *ring.Ring) (mat *Matrix) {
	mat = new(Matrix)
	mat.Elements = make([][]*ring.Poly, n)
	for i := 0; i < n; i++ {
		mat.Elements[i] = make([]*ring.Poly, m)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			mat.Elements[i][j] = r.NewPoly()
		}
	}
	return
}
func (mat *Matrix) NTT(n, m int, r *ring.Ring) (mat2 *Matrix) {
	mat2 = NewMatrix(n, m, r)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			r.NTT(mat.Elements[i][j], mat2.Elements[i][j])
		}
	}
	return
}

func (mat *Matrix) InvNTT(n, m int, r *ring.Ring) (mat2 *Matrix) {
	mat2 = NewMatrix(n, m, r)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			r.InvNTT(mat.Elements[i][j], mat2.Elements[i][j])
		}
	}
	return
}

func NewRandomMatrix(n, m int, r *ring.Ring, prng utils.PRNG) (mat *Matrix) {
	mat = new(Matrix)
	mat.Elements = make([][]*ring.Poly, n)
	for i := 0; i < n; i++ {
		mat.Elements[i] = make([]*ring.Poly, m)
	}
	s := ring.NewUniformSampler(prng, r)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			mat.Elements[i][j] = s.ReadNew()
		}
	}
	return
}

func (mat *Matrix) Transpose() (mat2 *Matrix) {
	n := len(mat.Elements)
	m := len(mat.Elements[0])
	mat2 = new(Matrix)
	mat2.Elements = make([][]*ring.Poly, m)
	for i := 0; i < m; i++ {
		mat2.Elements[i] = make([]*ring.Poly, n)
	}
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			mat2.Elements[i][j] = mat.Elements[j][i]
		}
	}
	return
}

// multiplies Mv = w
func (mat *Matrix) MulVecRight(vec1 *Vector, r *ring.Ring) (vec2 *Vector) {
	n := len(mat.Elements)
	m := len(mat.Elements[0])
	vec2 = new(Vector)
	vec2.Elements = make([]*ring.Poly, n)
	for i := 0; i < n; i++ {
		p1 := r.NewPoly()
		for j := 0; j < m; j++ {
			p2 := r.NewPoly()
			r.MulCoeffsConstant(mat.Elements[i][j], vec1.Elements[j], p2)
			r.Add(p1, p2, p1)
		}
		vec2.Elements[i] = p1
	}
	return
}

// vM = w
func (mat *Matrix) MulVecLeft(vec1 *Vector, r *ring.Ring) (vec2 *Vector) {
	n := len(mat.Elements)
	m := len(mat.Elements[0])
	vec2 = new(Vector)
	vec2.Elements = make([]*ring.Poly, m)
	for i := 0; i < m; i++ {
		p1 := r.NewPoly()
		for j := 0; j < n; j++ {
			p2 := r.NewPoly()
			r.MulCoeffs(mat.Elements[j][i], vec1.Elements[j], p2)
			r.Add(p1, p2, p1)
		}
		vec2.Elements[i] = p1
	}
	return
}
