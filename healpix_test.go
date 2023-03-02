package goheapix

import (
	"math"
	"testing"
)

var (
	eps = 1e-8
)

func TestVec2Ang(t *testing.T) {
	v := []struct {
		vec  [3]float64
		t, f float64
	}{
		{[3]float64{1, 0, 0}, 1.57079633, 0},
		{[3]float64{0, 1, 0}, 1.57079633, 1.57079633},
		{[3]float64{0, 0, 1}, 0, 0},
	}

	for _, x := range v {
		t1, f := Vec2Ang(x.vec)
		if math.Abs(t1-x.t) > eps || math.Abs(f-x.f) > eps {
			t.Errorf("%v %v %v ", x, t1, f)
		}
	}
}

func TestAng2Vec(t *testing.T) {
	v := []struct {
		vec  [3]float64
		t, f float64
	}{
		{[3]float64{1, 0, 0}, 1.57079633, 0},
		{[3]float64{0, 1, 0}, 1.57079633, 1.57079633},
		{[3]float64{0, 0, 1}, 0, 0},
	}

	for _, x := range v {
		vec := Ang2Vec(x.t, x.f)
		if math.Abs(vec[0]-x.vec[0]) > eps || math.Abs(vec[1]-x.vec[1]) > eps || math.Abs(vec[2]-x.vec[2]) > eps {
			t.Errorf("%v %v", x, vec)
		}
	}
}

func TestNpix2Nside(t *testing.T) {
	v := [][]int64{
		{786432, 256},
	}
	for _, x := range v {
		if Npix2Nside(x[0]) != x[1] {
			t.Errorf("%v", x)
		}
	}
}

func TestNside2Npix(t *testing.T) {
	v := [][]int64{
		{786432, 256},
	}
	for _, x := range v {
		if Nside2Npix(x[1]) != x[0] {
			t.Errorf("%v", x)
		}
	}
}

func TestAng2PixNest(t *testing.T) {
	v := []struct {
		nside int64
		t, f  float64
		pix   int64
	}{
		{1, π / 2.0, 0, 4},
		{2, π / 2.0, 0, 19},
		{4, π / 2.0, 0, 76},
		{16, π / 2.0, 0, 1216},
		{64, π / 2.0, 0, 19456},
		{256, π / 2.0, 0, 311296},
		{1024, π / 2.0, 0, 4980736},
		{4096, π / 2.0, 0, 79691776},
		{32768, π / 2.0, 0, 5100273664},      //order = 15
		{1048576, π / 2.0, 0, 5222680231936}, //order = 20
	}

	for _, x := range v {
		t.Logf("%v", Ang2PixNest(x.nside, x.t, x.f))
		if Ang2PixNest(x.nside, x.t, x.f) != x.pix {
			t.Errorf("%v", x)
		}

	}
}

func TestPix2AngNest(t *testing.T) {
	v := []struct {
		nside int64
		t, f  float64
		pix   int64
	}{
		{1, 1.5707963267948966, 0.0, 4},
		{2, 1.2309594173407747, 0.0, 19},
		{4, 1.4033482475752073, 0.0, 76},
		{16, 1.5291175943723188, 0.0, 1216},
		{64, 1.5603794717389192, 0.0, 19456},
		{256, 1.5681921571847817, 0.0, 311296},
		{1024, 1.5701452850822386, 0.0, 4980736},
		{4096, 1.5706335663775113, 0.0, 79691776},
		{32768, 1.5707759817428117, 0.0, 5100273664},      //order = 15
		{1048576, 1.5707956910120189, 0.0, 5222680231936}, // order=20
	}

	for _, x := range v {
		if θ, φ := Pix2AngNest(x.nside, x.pix); θ != x.t || φ != x.f {
			t.Errorf("%v, %v, %v", x, θ, φ)
		}

	}

}

func TestVec2PixNest(t *testing.T) {
	v := []struct {
		nside int64
		v     [3]float64
		pix   int64
	}{
		{1, [3]float64{1, 0, 0}, 4},
		{2, [3]float64{1, 0, 0}, 17},
		{4, [3]float64{1, 0, 0}, 70},
		{16, [3]float64{1, 0, 0}, 1130},
		{64, [3]float64{1, 0, 0}, 18090},
		{256, [3]float64{1, 0, 0}, 289450},
		{1024, [3]float64{1, 0, 0}, 4631210},
		{4096, [3]float64{1, 0, 0}, 74099370},
		{32768, [3]float64{1, 0, 0}, 4742359722},
		{1048576, [3]float64{1, 0, 0}, 4856176356010},
	}

	for _, x := range v {
		t.Logf("%v", Vec2PixNest(x.nside, x.v))
		if Vec2PixNest(x.nside, x.v) != x.pix {
			t.Errorf("%v", x)
		}

	}
}

func TestPix2VecNest(t *testing.T) {
	v := []struct {
		nside int64
		v     [3]float64
		pix   int64
	}{
		{1, [3]float64{1.0, 0.0, 0.0}, 4},
		{2, [3]float64{0.9238795325112867, 0.3826834323650898, 0.0}, 17},
		{4, [3]float64{0.9807852804032304, 0.19509032201612825, 0.0}, 70},
		{16, [3]float64{0.9987954562051724, 0.049067674327418015, 0.0}, 1130},
		{64, [3]float64{0.9999247018391445, 0.012271538285719925, 0.0}, 18090},
		{256, [3]float64{0.9999952938095762, 0.003067956762965976, 0.0}, 289450},
		{1024, [3]float64{0.9999997058628822, 0.0007669903187427045, 0.0}, 4631210},
		{4096, [3]float64{0.9999999816164293, 0.0001917475973107033, 0.0}, 74099370},
		{32768, [3]float64{0.9999999997127567, 2.396844980841822e-05, 0.0}, 4742359722},
		{1048576, [3]float64{0.9999999999997194, 7.490140565847157e-07, 0.0}, 4856176356010},
	}

	for _, x := range v {
		//t.Logf("%v", Pix2VecNest(x.nside, x.pix))
		vec := Pix2VecNest(x.nside, x.pix)
		if math.Abs(vec[0]-x.v[0]) > eps || math.Abs(vec[1]-x.v[1]) > eps || math.Abs(vec[2]-x.v[2]) > eps {
			t.Errorf("%v ", vec)
		} else {
			t.Logf("%v", x)
		}

	}
}
