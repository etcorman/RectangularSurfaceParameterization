// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#pragma once

#if !defined(__SSE2__) && !defined(__x86_64__) && !defined(__amd64) && !defined(_M_AMD64) && !defined(_M_X64)
#include "vectorize_arm.h"
#else

#include <new>
#include <cstdint>
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__)
#include <x86intrin.h>
#else
#error Unknown compiler! Please contact author at yoann.coudert-osmont@univ-lorraine.fr
#endif

struct v2d;
struct v4d;

struct v4i {
friend v4d;
protected:
	__m128i xmm;
public:
	inline static constexpr std::align_val_t ALI{16};
	v4i(const int *x, const bool aligned=false):
		xmm(aligned ? _mm_load_si128((const __m128i*)x) : _mm_loadu_si128((const __m128i_u*)x)) {}
};

struct v2l {
friend v2d;
protected:
	__m128i xmm;
public:
	inline static constexpr std::align_val_t ALI{16};
	v2l(const int64_t *x, const bool aligned=false):
		xmm(aligned ? _mm_load_si128((const __m128i*)x) : _mm_loadu_si128((const __m128i_u*)x)) {}
};

struct v2d {
friend v4d;
protected:
	__m128d xmm;
	v2d(const __m128d &x): xmm(x) {}
public:
	inline static constexpr std::align_val_t ALI{16};
	v2d(double x): xmm(_mm_set1_pd(x)) {}
	v2d(double x0, double x1): xmm(_mm_setr_pd(x0, x1)) {}
	v2d(const double *x, const bool aligned=false): xmm(aligned ? _mm_load_pd(x) : _mm_loadu_pd(x)) {}
	#ifdef __AVX2__
	v2d(const double* array, const v2l &indices): xmm(_mm_i64gather_pd(array, indices.xmm, 8)) {}
	#else
	v2d(const double* array, const v2l &indices) {
		int64_t inds[2];
		_mm_store_si128((__m128i*)inds, indices.xmm);
		xmm = _mm_setr_pd(array[inds[0]], array[inds[1]]);
	}
	#endif

	inline void store(double *x, const bool aligned=false) const { if(aligned) _mm_store_pd(x, xmm); else _mm_storeu_pd(x, xmm); }

	inline double sum() const {
	#if defined(__AVX__) || TRUE
		return _mm_cvtsd_f64(_mm_add_pd(xmm, _mm_unpackhi_pd(xmm, xmm)));
	#else // In SEE Shuffle has an additional port at the cost of one more instruction byte
		return _mm_cvtsd_f64(_mm_add_pd(xmm, _mm_shuffle_pd(xmm, xmm, 0b01)));
	#endif
	}

	inline v2d operator+(const v2d &v) const { return { _mm_add_pd(xmm, v.xmm) }; }
	inline v2d& operator+=(const v2d &v) { return *this = *this + v; }
	inline v2d operator-(const v2d &v) const { return { _mm_sub_pd(xmm, v.xmm) }; }
	inline v2d& operator-=(const v2d &v) { return *this = *this - v; }
	inline v2d operator*(const v2d &v) const { return { _mm_mul_pd(xmm, v.xmm) }; }
	inline v2d& operator*=(const v2d &v) { return *this = *this * v; }

	friend v2d fma(const v2d &a, const v2d &b, const v2d &c) {
	#if defined(__FMA__)
		return _mm_fmadd_pd(a.xmm, b.xmm, c.xmm);
	#elif defined(__FMA4__)
		return _mm_macc_pd(a.xmm, b.xmm, c.xmm);
	#else
		return a * b + c;
	#endif
	}
};

#ifdef __AVX__

struct v4l {
friend v4d;
protected:
	__m256i ymm;
public:
	inline static constexpr std::align_val_t ALI{32};
	v4l(const int64_t *y, const bool aligned=false):
		ymm(aligned ? _mm256_load_si256((const __m256i*)y) : _mm256_loadu_si256((const __m256i_u*)y)) {}
};

struct v4d {
protected:
	__m256d ymm;
	v4d(const __m256d &y): ymm(y) {}
public:
	inline static constexpr std::align_val_t ALI{32};
	v4d(double y): ymm(_mm256_set1_pd(y)) {}
	v4d(double y0, double y1, double y2, double y3): ymm(_mm256_setr_pd(y0, y1, y2, y3)) {}
	v4d(const double *y, const bool aligned=false): ymm(aligned ? _mm256_load_pd(y) : _mm256_loadu_pd(y)) {}
	#ifdef __AVX2__
	v4d(const double* array, const v4i &indices): ymm(_mm256_i32gather_pd(array, indices.xmm, 8)) {}
	v4d(const double* array, const v4l &indices): ymm(_mm256_i64gather_pd(array, indices.ymm, 8)) {}
	#else
	v4d(const double* array, const v4i &indices) {
		int32_t inds[4];
		_mm_store_si128((__m128i*)inds, indices.xmm);
		ymm = _mm128_setr_pd(array[inds[0]], array[inds[1]], array[inds[2]], array[inds[3]]);
	}
	v4d(const double* array, const v4l &indices) {
		int64_t inds[4];
		_mm256_store_si256((__m256i*)inds, indices.ymm);
		ymm = _mm256_setr_pd(array[inds[0]], array[inds[1]], array[inds[2]], array[inds[3]]);
	}
	#endif

	inline void store(double *y, const bool aligned=false) const { if(aligned) _mm256_store_pd(y, ymm); else _mm256_storeu_pd(y, ymm); }

	inline v2d get_low() const { return _mm256_castpd256_pd128(ymm); }
	inline v2d get_high() const { return _mm256_extractf128_pd(ymm, 1); }

	inline double sum() const { return (get_low() + get_high()).sum(); }

	inline v4d operator+(const v4d &v) const { return _mm256_add_pd(ymm, v.ymm); }
	inline v4d& operator+=(const v4d &v) { return *this = *this + v; }
	inline v4d operator-(const v4d &v) const { return _mm256_sub_pd(ymm, v.ymm); }
	inline v4d& operator-=(const v4d &v) { return *this = *this - v; }
	inline v4d operator*(const v4d &v) const { return _mm256_mul_pd(ymm, v.ymm); }
	inline v4d& operator*=(const v4d &v) { return *this = *this * v; }

	#if defined(__FMA__)
	inline friend v4d fma(const v4d &a, const v4d &b, const v4d &c) { return _mm256_fmadd_pd(a.ymm, b.ymm, c.ymm); }
	inline friend v4d fnma(const v4d &a, const v4d &b, const v4d &c) { return _mm256_fnmadd_pd(a.ymm, b.ymm, c.ymm); }
	#elif defined(__FMA4__)
	inline friend v4d fma(const v4d &a, const v4d &b, const v4d &c) { return _mm256_macc_pd(a.ymm, b.ymm, c.ymm); }
	inline friend v4d fnma(const v4d &a, const v4d &b, const v4d &c) { return _mm256_nmacc_pd(a.ymm, b.ymm, c.ymm); }
	#else
	inline friend v4d fma(const v4d &a, const v4d &b, const v4d &c) { return a * b + c; }
	inline friend v4d fnma(const v4d &a, const v4d &b, const v4d &c) { return c - a * b; }
	#endif
};


#else

struct v4l {
friend v4d;
protected:
	__m128i xmm0, xmm1;
public:
	inline static constexpr std::align_val_t ALI{16};
	v4l(const int64_t *y, const bool aligned=false):
		xmm0(aligned ? _mm_load_si128((const __m128i*)y) : _mm_loadu_si128((const __m128i_u*)y)),
		xmm1(aligned ? _mm_load_si128((const __m128i*)(y+2)) : _mm_loadu_si128((const __m128i_u*)(y+2))) {}
};

struct v4d {
protected:
	__m128d xmm0, xmm1;
	v4d(const v2d x0, const v2d x1): xmm0(x0.xmm), xmm1(x1.xmm) {}
public:
	inline static constexpr std::align_val_t ALI{16};
	v4d(double y):
		xmm0(_mm_set1_pd(y)),
		xmm1(_mm_set1_pd(y)) {}
	v4d(double y0, double y1, double y2, double y3):
		xmm0(_mm_setr_pd(y0, y1)),
		xmm1(_mm_setr_pd(y2, y3)) {}
	v4d(const double *y, const bool aligned=false):
		xmm0(aligned ? _mm_load_pd(y) : _mm_loadu_pd(y)),
		xmm1(aligned ? _mm_load_pd(y+2) : _mm_loadu_pd(y+2)) {}
	v4d(const double* array, const v4i &indices) {
		int32_t inds[4];
		_mm_store_si128((__m128i*)inds, indices.xmm);
		xmm0 = _mm_setr_pd(array[inds[0]], array[inds[1]]);
		xmm1 = _mm_setr_pd(array[inds[2]], array[inds[3]]);
	}
	v4d(const double* array, const v4l &indices) {
		int64_t inds[4];
		_mm_store_si128((__m128i*)inds, indices.xmm0);
		_mm_store_si128((__m128i*)(inds+2), indices.xmm1);
		xmm0 = _mm_setr_pd(array[inds[0]], array[inds[1]]);
		xmm1 = _mm_setr_pd(array[inds[2]], array[inds[3]]);
	}

	inline void store(double *y, const bool aligned=false) const {
		if(aligned) { _mm_store_pd(y, xmm0); _mm_store_pd(y+2, xmm1); }
		else { _mm_storeu_pd(y, xmm0); _mm_storeu_pd(y+2, xmm1); }
	}

	inline v2d get_low() const { return xmm0; }
	inline v2d get_high() const { return xmm1; }

	inline double sum() const { return (get_low() + get_high()).sum(); }

	inline v4d operator+(const v4d &v) const { return v4d(_mm_add_pd(xmm0, v.xmm0), _mm_add_pd(xmm1, v.xmm1)); }
	inline v4d& operator+=(const v4d &v) { return *this = *this + v; }
	inline v4d operator-(const v4d &v) const { return v4d(_mm_sub_pd(xmm0, v.xmm0), _mm_sub_pd(xmm1, v.xmm1)); }
	inline v4d& operator-=(const v4d &v) { return *this = *this - v; }
	inline v4d operator*(const v4d &v) const { return v4d(_mm_mul_pd(xmm0, v.xmm0), _mm_mul_pd(xmm1, v.xmm1)); }
	inline v4d& operator*=(const v4d &v) { return *this = *this * v; }

	#if defined(__FMA__)
	inline friend v4d fma(const v4d &a, const v4d &b, const v4d &c) { return v4d(_mm_fmadd_pd(a.xmm0, b.xmm0, c.xmm0), _mm_fmadd_pd(a.xmm1, b.xmm1, c.xmm1)); }
	inline friend v4d fnma(const v4d &a, const v4d &b, const v4d &c) { return v4d(_mm_fnmadd_pd(a.xmm0, b.xmm0, c.xmm0), _mm_fnmadd_pd(a.xmm1, b.xmm1, c.xmm1)); }
	#elif defined(__FMA4__)
	inline friend v4d fma(const v4d &a, const v4d &b, const v4d &c) { return v4d(_mm_macc_pd(a.xmm0, b.xmm0, c.xmm0), _mm_macc_pd(a.xmm1, b.xmm1, c.xmm1)); }
	inline friend v4d fnma(const v4d &a, const v4d &b, const v4d &c) { return v4d(_mm_nmacc_pd(a.xmm0, b.xmm0, c.xmm0), _mm_nmacc_pd(a.xmm1, b.xmm1, c.xmm1)); }
	#else
	inline friend v4d fma(const v4d &a, const v4d &b, const v4d &c) { return a * b + c; }
	inline friend v4d fnma(const v4d &a, const v4d &b, const v4d &c) { return c - a * b; }
	#endif
};

#endif

#endif