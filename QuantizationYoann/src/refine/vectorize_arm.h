// Copyright (C) 2024, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#pragma once

#include <cstdint>
#include <algorithm>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

struct v2d;
struct v4d;

struct v4i {
friend v4d;
protected:
	int32_t m[4];
public:
	inline static constexpr std::align_val_t ALI{alignof(int32_t)};
	v4i(const int *x, const bool aligned=false) { std::copy_n(x, 4, m); }
};

struct v2l {
friend v2d;
protected:
	int64_t mm[2];
public:
	inline static constexpr std::align_val_t ALI{alignof(int64_t)};
	v2l(const int64_t *x, const bool aligned=false) { std::copy_n(x, 2, mm); }
};

struct v2d {
friend v4d;
protected:
	double mm[2];
public:
	inline static constexpr std::align_val_t ALI{alignof(double)};
	v2d(double x): mm{x,x} {}
	v2d(double x0, double x1): mm{x0,x1} {}
	v2d(const double *x, const bool aligned=false) { std::copy_n(x, 2, mm); }
	v2d(const double* array, const v2l &indices): mm{array[indices.mm[0]], array[indices.mm[1]]} {}

	inline void store(double *x, const bool aligned=false) const { std::copy_n(mm, 2, x); }

	inline double sum() const { return mm[0]+mm[1]; }

	inline v2d operator+(const v2d &v) const { return { mm[0]+v.mm[0], mm[1]+v.mm[1] }; }
	inline v2d& operator+=(const v2d &v) { return *this = *this + v; }
	inline v2d operator-(const v2d &v) const { return { mm[0]-v.mm[0], mm[1]-v.mm[1] }; }
	inline v2d& operator-=(const v2d &v) { return *this = *this - v; }
	inline v2d operator*(const v2d &v) const { return { mm[0]*v.mm[0], mm[1]*v.mm[1] }; }
	inline v2d& operator*=(const v2d &v) { return *this = *this * v; }

	friend v2d fma(const v2d &a, const v2d &b, const v2d &c) { return a * b + c; }
};

struct v4l {
friend v4d;
protected:
	int64_t mm[4];
public:
	inline static constexpr std::align_val_t ALI{alignof(int64_t)};
	v4l(const int64_t *y, const bool aligned=false) { std::copy_n(y, 4, mm); }
};

struct v4d {
protected:
	double mm[4];
public:
	inline static constexpr std::align_val_t ALI{alignof(double)};
	v4d(double y): mm{y,y,y,y} {}
	v4d(double y0, double y1, double y2, double y3): mm{y0,y1,y2,y3} {}
	v4d(const double *y, const bool aligned=false) { std::copy_n(y, 4, mm); }
	v4d(const double* array, const v4i &indices):
		mm{array[indices.m[0]], array[indices.m[1]], array[indices.m[2]], array[indices.m[3]]} {}
	v4d(const double* array, const v4l &indices):
		mm{array[indices.mm[0]], array[indices.mm[1]], array[indices.mm[2]], array[indices.mm[3]]} {}

	inline void store(double *y, const bool aligned=false) const { std::copy_n(mm, 4, y); }

	inline v2d get_low() const { return v2d(mm[0], mm[1]); }
	inline v2d get_high() const { return v2d(mm[2], mm[3]); }

	inline double sum() const { return mm[0]+mm[1]+mm[2]+mm[3]; }

	inline v4d operator+(const v4d &v) const { return { mm[0]+v.mm[0], mm[1]+v.mm[1], mm[2]+v.mm[2], mm[3]+v.mm[3] }; }
	inline v4d& operator+=(const v4d &v) { return *this = *this + v; }
	inline v4d operator-(const v4d &v) const { return { mm[0]-v.mm[0], mm[1]-v.mm[1], mm[2]-v.mm[2], mm[3]-v.mm[3] }; }
	inline v4d& operator-=(const v4d &v) { return *this = *this - v; }
	inline v4d operator*(const v4d &v) const { return { mm[0]*v.mm[0], mm[1]*v.mm[1], mm[2]*v.mm[2], mm[3]*v.mm[3] }; }
	inline v4d& operator*=(const v4d &v) { return *this = *this * v; }

	inline friend v4d fma(const v4d &a, const v4d &b, const v4d &c) { return a * b + c; }
	inline friend v4d fnma(const v4d &a, const v4d &b, const v4d &c) { return c - a * b; }
};

#pragma GCC diagnostic pop