package goheapix

import (
	"fmt"
	"math"
)

func fmodulo(v1, v2 float64) float64 {
	if v1 >= 0 {
		if v1 < v2 {
			return v1
		}
		return math.Mod(v1, v2)
	}
	tmp := math.Mod(v1, v2) + v2
	if tmp == v2 {
		return 0.
	}
	return tmp
}

func imodulo(v1, v2 int64) int64 {
	v := v1 % v2
	if v >= 0 {
		return v
	}
	return v + v2
}

func isqrt(v int64) int64 {
	return int64(math.Sqrt(float64(v) + 0.5))
}

func spreadBits(v int64) int64 {
	return int64(utab[v&0xff]) | (int64(utab[(v>>8)&0xff]) << 16) | (int64(utab[(v>>16)&0xff]) << 32) | (int64(utab[(v>>24)&0xff]) << 48)
}

func compressBits64(v int64) int64 {
	var raw int64
	raw = v & 0x5555555555555555
	raw |= raw >> 15
	return ctab[raw&0xff] | (ctab[(raw>>8)&0xff] << 4) | (ctab[(raw>>32)&0xff] << 16) | (ctab[(raw>>40)&0xff] << 20)
}

func validTheta(θ float64) error {
	if θ < 0 || θ > π {
		return fmt.Errorf("θ = %v out of range, must 0 <= θ <= π", θ)
	}
	return nil
}

func r2d(r float64) float64 {
	return r * 180.0 / π
}

func d2r(d float64) float64 {
	return d * π / 180.0
}
