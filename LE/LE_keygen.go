package LE

import "lattice-LE-go/matrix"

/*
KeyGen this is the KGen algorithm
Inputs:

	le: public parameters of the system

Outputs:

	pk: public key in normal form
	sk: secret key in NTT form.
*/
func (le *LE) KeyGen() (pk, sk *matrix.Vector) {
	sk = matrix.NewRandomVecBinary(le.M2, le.R)
	sk = sk.NTT(le.R)
	pk = le.BNTT.MulVecRight(sk, le.R)
	pk = pk.InvNTT(le.R)
	return pk, sk
}
