package goheapix

import (
	"math"
)

//goheapix Healpix Go port from: chealpix(http://healpix.sourceforge.net/index.php)

//Vec2Ang (x, y, z) -> (θ, φ)
func Vec2Ang(v [3]float64) (θ, φ float64) {
	θ = math.Atan2((math.Sqrt(v[0]*v[0] + v[1]*v[1])), v[2])
	φ = math.Atan2(v[1], v[0])
	if φ < 0.0 {
		φ += π2
	}
	return
}

//Ang2Vec (θ, φ) -> (x, y, z)
func Ang2Vec(θ, φ float64) (v [3]float64) {
	sz := math.Sin(θ)
	v[0] = sz * math.Cos(φ)
	v[1] = sz * math.Sin(φ)
	v[2] = math.Cos(θ)
	return
}

//Npix2Nside npix to nsize
func Npix2Nside(npix int64) int64 {
	res := isqrt(npix / 12)
	if res*res*12 != npix {
		res = -1
	}
	return res
}

//Nside2Npix nside to npix
func Nside2Npix(nside int64) int64 {
	return 12 * nside * nside
}

//----------------

// Pix2AngNest pix -> ang @ nest
func Pix2AngNest(nside int64, ipix int64) (θ, φ float64) {
	var z, s float64
	z, s, φ = pix2angNestZPhi(nside, ipix)
	if s < -2. {
		θ = math.Acos(z)
	} else {
		θ = math.Atan2(s, z)
	}

	return
}

//Pix2VecNest pix -> vec @ nest
func Pix2VecNest(nside int64, ipix int64) (v [3]float64) {
	z, s, φ := pix2angNestZPhi(nside, ipix)
	if s < -2. {
		s = math.Sqrt((1. - z) * (1. + z))
	}
	v[0] = s * math.Cos(φ)
	v[1] = s * math.Sin(φ)
	v[2] = z

	return
}

//Ang2PixNest ang -> pix @ nest
func Ang2PixNest(nside int64, θ, φ float64) (ipnest int64) {
	if err := validTheta(θ); err != nil {
		ipnest = -1
	}
	ipnest = ang2pixNestZPhi(nside, math.Cos(θ), φ)

	return
}

//Vec2PixNest Vec -> pix @ nest
func Vec2PixNest(nside int64, v [3]float64) (ipnest int64) {
	vlen := math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
	ipnest = ang2pixNestZPhi(nside, v[2]/vlen, math.Atan2(v[1], v[0]))
	return
}

func ang2pixNestZPhi(nside int64, z, φ float64) int64 {
	za := math.Abs(z)
	tt := fmodulo(φ, π2) * invhπ // in [0,4)
	//tt := fmodulo(φ*invhπ, 4.0) // in [0,4)
	var faceNum, ix, iy int64

	if za <= twothird { //Equatorial region
		temp1 := float64(nside) * (0.5 + tt)
		temp2 := float64(nside) * (z * 0.75)
		jp := int64(temp1 - temp2)
		jm := int64(temp1 + temp2)
		ifp := jp / nside
		ifm := jm / nside
		if ifp == ifm { /* faces 4 to 7 */
			if ifp == 4 {
				faceNum = 4
			} else {
				faceNum = ifp + 4
			}
		} else if ifp < ifm { /* (half-)faces 0 to 3 */
			faceNum = ifp
		} else { /* (half-)faces 8 to 11 */
			faceNum = ifm + 8

		}
		ix = jm & (nside - 1)
		iy = nside - (jp & (nside - 1)) - 1
	} else { /* polar region, za > 2/3 */
		ntt := int64(tt)
		var jp, jm int64
		var tp, tmp float64
		if ntt >= 4 {
			ntt = 3
		}
		tp = tt - float64(ntt)

		if za > 0.99 { // 64 bits
			sth := math.Sqrt(1.0 - za*za) // sin(θ)
			tmp = float64(nside) * sth / math.Sqrt((1.0+za)/3.)
		} else {
			tmp = float64(nside) * math.Sqrt(3.0*(1.0-za))
		}

		jp = int64(tp * tmp)
		jm = int64((1.0 - tp) * tmp)
		if jp >= nside {
			jp = nside - 1
		}
		if jm >= nside {
			jm = nside - 1
		}
		if z >= 0 {
			faceNum = ntt /* in {0,3} */
			ix = nside - jm - 1
			iy = nside - jp - 1
		} else {
			faceNum = ntt + 8 /* in {8,11} */
			ix = jp
			iy = jm
		}
	}
	return xyf2nest(nside, ix, iy, faceNum)
}

func xyf2nest(nside int64, ix, iy, faceNum int64) int64 {
	return (faceNum * nside * nside) + spreadBits(ix) + (spreadBits(iy) << 1)
}

func pix2angNestZPhi(nside, ipix int64) (z, s, φ float64) {
	nl4 := nside * 4
	npix := 12 * nside * nside
	fact2 := 4. / float64(npix)
	s = -5.
	ix, iy, faceNum := nest2xyf(nside, ipix)

	jr := (jrll[faceNum] * nside) - ix - iy - 1
	var kshift, jp int64
	var nr float64
	if jr < nside {
		nr = float64(jr)
		tmp := nr * nr * fact2
		z = 1 - tmp
		if z > 0.99 {
			s = math.Sqrt(tmp * (2. - tmp))
		}
		kshift = 0
	} else if jr > 3*nside {
		nr = float64(nl4 - jr)
		tmp := nr * nr * fact2
		z = tmp - 1
		if z < 0.99 {
			s = math.Sqrt(tmp * (2. - tmp))
		}

	} else {
		fact1 := float64(nside<<1) * fact2
		nr = float64(nside)
		z = float64(2*nside-jr) * fact1
		kshift = (jr - nside) & 1
	}

	jp = (jpll[faceNum]*int64(nr) + ix - iy + 1 + kshift) / 2
	if jp > nl4 {
		jp -= nl4
	}
	if jp < 1 {
		jp += nl4
	}
	φ = (float64(jp) - float64(kshift+1)*0.5) * (hπ / nr)
	return
}

func nest2xyf(nside, pix int64) (ix, iy, faceNum int64) {
	npFace := nside * nside
	faceNum = pix / npFace
	pix &= (npFace - 1)
	ix = compressBits64(pix)
	iy = compressBits64(pix >> 1)
	return
}
