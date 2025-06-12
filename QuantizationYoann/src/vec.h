// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#pragma once

#include <iostream>
#include <cmath>

#define M_PI_2 1.57079632679489661923
#define M_PI   3.14159265358979323846

template<typename T, typename Enable=void>
struct get_scalar { using value = typename T::Scalar; };
template<typename T>
struct get_scalar<T, std::enable_if_t<std::is_arithmetic_v<T>>> { using value = T; };
template<typename T>
using get_scalar_v = typename get_scalar<T>::value;

template<int N, typename T, typename V>
struct vec_base {
using Scalar = get_scalar_v<T>;
#define _vecFOR for(int i = 0; i < N; ++i)
private:
	// Casting to V
	inline V& toV() { return *reinterpret_cast<V*>(this); }
	inline const V& toV() const { return *reinterpret_cast<const V*>(this); }

public:
	// Constructor
	vec_base() = default;
	template<typename U, typename W>
	vec_base(const vec_base<N, U, W> &other) { _vecFOR (*this)[i] = other[i]; }

	// Coordinate access
	inline const T& operator[](int i) const { return reinterpret_cast<const T*>(this)[i]; }
	inline T& operator[](int i) { return reinterpret_cast<T*>(this)[i]; }

	// Unary minus
	inline V operator-() const { V v; _vecFOR v[i] = -(*this)[i]; return v; }

	// Operators with scalar
	inline V& operator*=(Scalar x) { _vecFOR (*this)[i] *= x; return toV(); }
	inline V operator*(Scalar x) const { return V(toV()) *= x; }
	inline friend V operator*(Scalar x, const vec_base &v) { return v * x; }
	inline V& operator/=(Scalar x) { _vecFOR (*this)[i] /= x; return toV(); }
	inline V operator/(Scalar x) const { return V(toV()) /= x; }

	// Operators with vectors
	template<typename W>
	inline V& operator+=(const W& other) { _vecFOR (*this)[i] += other[i]; return toV(); }
	template<typename U, typename W>
	inline auto operator+(const vec_base<N, U, W>& other) const {
		using TU = std::common_type_t<T, U>;
		if constexpr (std::is_same_v<TU, T>) return V(toV()) += other;
		if constexpr (std::is_same_v<TU, U>) return W(*this) += other;
	}
	template<typename W>
	inline V& operator-=(const W& other) { _vecFOR (*this)[i] -= other[i]; return toV(); }
	template<typename U, typename W>
	inline auto operator-(const vec_base<N, U, W>& other) const {
		using TU = std::common_type_t<T, U>;
		if constexpr (std::is_same_v<TU, T>) return V(toV()) -= other;
		if constexpr (std::is_same_v<TU, U>) return W(*this) -= other;
	}

	// Dot product
	template<typename U, typename W>
	inline auto operator*(const vec_base<N, U, W>& other) const { std::common_type_t<T, U> d = 0.; _vecFOR d += (*this)[i]*other[i]; return d; }

	// Equality
	inline bool operator==(const V& other) { _vecFOR if((*this)[i] != other[i]) return false; return true; }
	inline bool operator!=(const V& other) { _vecFOR if((*this)[i] != other[i]) return true; return false; }

	// Norm
	inline Scalar norm2() const { Scalar n2 = 0.; _vecFOR n2 += (*this)[i]*(*this)[i]; return n2; }
	inline Scalar norm() const { return std::sqrt(norm2()); }
	inline V& normalize() { return *this *= 1. / norm(); }

	// Print
	friend std::ostream& operator<<(std::ostream &stream, const vec_base &v) {
		stream << '(' << v[0];
		for(int i = 1; i < N; ++i) stream << ", " << v[i];
		return stream << ')';
	}
#undef _vecFOR
};

template <typename T>
struct vec2T : vec_base<2, T, vec2T<T>> {
	T x, y;
	vec2T() = default;
	vec2T(T x, T y=0.): x(x), y(y) {}
	template<typename U, typename W>
	vec2T(const vec_base<2, U, W> &other): vec_base<2, T, vec2T>(other) {}
	vec2T& rotate90(int r) {
		if(r&1) std::swap(x, y);
		if((r+1)&2) x = -x;
		if(r&2) y = -y;
		return *this;
	}
};
typedef vec2T<double> vec2;
typedef vec2T<int> vec2i;
typedef vec2T<bool> vec2b;

struct vec3 : vec_base<3, double, vec3> {
	double x, y, z;
	vec3() = default;
	vec3(double x, double y=0., double z=0.): x(x), y(y), z(z) {}
	template<typename U, typename W>
	vec3(const vec_base<3, U, W> &other): vec_base<3, double, vec3>(other) {}
};

struct mat2 : vec_base<2, vec2, mat2> {
	vec2 rows[2];
	mat2() = default;
	mat2(const vec2 &row0, const vec2 &row1): rows{row0, row1} {}
	template<typename U, typename W>
	mat2(const vec_base<2, U, W> &other): vec_base<2, vec2, mat2>(other) {}

	inline double operator()(int i, int j) const { return rows[i][j]; }
	inline double& operator()(int i, int j) { return rows[i][j]; }

	inline double det() const { return rows[0][0]*rows[1][1] - rows[0][1]*rows[1][0]; }
	inline mat2 com() const { return mat2(vec2(rows[1][1], -rows[1][0]), vec2(-rows[0][1], rows[0][0])); }
	inline mat2 inv() const { const double d = 1./det(); return mat2(vec2(d*rows[1][1], -d*rows[0][1]), vec2(-d*rows[1][0], d*rows[0][0])); }
};