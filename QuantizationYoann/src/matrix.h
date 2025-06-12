// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#pragma once

#include <algorithm>
#include <ostream>
#include <vector>

template<typename T>
struct SparseVector : std::vector<std::pair<int, T>> {
	using Base = std::vector<std::pair<int, T>>;
	using Scalar = T;
	using value_type = std::pair<int, T>;
	using Base::size; using Base::begin; using Base::end;
	inline SparseVector(std::size_t size=0): Base(size) {}

	SparseVector operator+(const SparseVector &v) const;
	SparseVector operator-(const SparseVector &v) const;
	inline SparseVector& operator+=(const SparseVector &v) { return *this = move(*this+v); }
	inline SparseVector& operator-=(const SparseVector &v) { return *this = move(*this-v); }
	SparseVector operator-() const;
	SparseVector operator*(T a) const;
	inline friend SparseVector operator*(T a, const SparseVector &v) { return v * a; } 
		
	inline friend std::ostream& operator<<(std::ostream &stream, const SparseVector &v) {
		if(!v.empty()) for(const auto &[i, a] : v) {
			if(i != v.begin()->first) stream << " + ";
			stream << a << " X" << i;
		} else stream << '0';
		return stream;
	}

	static inline bool isNull(T a) {
		if constexpr (std::is_integral_v<T>) return !a;
		if constexpr (std::is_floating_point_v<T>) return std::abs(a) < 1e-6;
		static_assert(true, "Bad template type for SparseVector (must be integer or floating point)");
	}

private:
	static inline value_type neg(const value_type &a) { return {a.first, -a.second}; };
};

struct SparseMatrix {
	using Row = SparseVector<int>;

	inline void setIdentity(int n) {
		_rows.resize(n);
		for(int i = 0; i < n; ++i) {
			_rows[i].clear();
			_rows[i].emplace_back(i, 1);
		}
	}

	inline bool empty() const { return _rows.empty(); }
	inline std::size_t size() const { return _rows.size(); }
	inline void resize(int n) { _rows.resize(n); }

	inline const std::vector<Row>& rows() const { return _rows; }
	inline const Row& getRow(int i) const { return _rows[i]; }
	inline Row& getRow(int i) { return _rows[i]; }

protected:
	std::vector<Row> _rows;
};


//===========================================
//== IMPLEMENTATION OF TEMPLATED FUNCTIONS ==
//===========================================


template<typename T>
SparseVector<T> SparseVector<T>::operator+(const SparseVector &r) const {
	SparseVector s; s.reserve(size() + r.size());
	auto i = begin(), j = r.begin();
	while(i != end() && j != r.end()) {
		if(i->first < j->first) s.push_back(*(i++));
		else if(i->first > j->first) s.push_back(*(j++));
		else {
			const T v = i->second + j->second;
			if(!isNull(v)) s.emplace_back(i->first, v);
			++i; ++j;
		}
	}
	s.insert(s.end(), i, end());
	s.insert(s.end(), j, r.end());
	return s;
}

template<typename T>
SparseVector<T> SparseVector<T>::operator-(const SparseVector &r) const {
	SparseVector s; s.reserve(size() + r.size());
	auto i = begin(), j = r.begin();
	while(i != end() && j != r.end()) {
		if(i->first < j->first) s.push_back(*(i++));
		else if(i->first > j->first) s.emplace_back(move(neg(*(j++))));
		else {
			const T v = i->second - j->second;
			if(!isNull(v)) s.emplace_back(i->first, v);
			++i; ++j;
		}
	}
	s.insert(s.end(), i, end());
	std::transform(j, r.end(), std::back_inserter(s), neg);
	return s;
}

template<typename T>
SparseVector<T> SparseVector<T>::operator-() const {
	SparseVector r(size());
	std::transform(begin(), end(), r.begin(), neg);
	return r;
}

template<typename T>
SparseVector<T> SparseVector<T>::operator*(T a) const {
	SparseVector r(*this);
	for(auto &[i, b] : r) b *= a; 
	return r;
}