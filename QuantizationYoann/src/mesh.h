// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#pragma once

#include "vec.h"

#include <vector>
#include <string>
#include <tuple>

struct Range {
	const int a, b;
	inline Range(int a, int b): a(a), b(b) {}
	inline Range(int b): a(0), b(b) {}
	struct iterator {
		int i;
		inline iterator(int i): i(i) {}
		inline iterator& operator++() { ++i; return *this; }
		inline bool operator!=(const iterator& rhs) const { return i != rhs.i; }
		inline const int &operator*() const { return i; }
	};
	inline iterator begin() const { return {a}; }
	inline iterator end  () const { return {b}; }
};
template<> struct std::iterator_traits<Range::iterator> {
	using iterator_category = std::forward_iterator_tag;
};

struct Mesh {
	std::vector<vec3> points;
	std::vector<int> h2v;

	inline int nverts() const { return points.size(); }
	inline int ncorners() const { return h2v.size(); }
	inline int nfacets() const { return h2v.size()/3; }
	inline Range verts() const { return Range(nverts()); }
	inline Range corners() const { return Range(ncorners()); }
	inline Range facets() const { return Range(nfacets()); }

	// Half-edge navigation
	inline int v2h(int v) const { return _v2h[v]; }
	inline int next(int h) const { return h%3==2 ? h-2 : h+1; }
	inline int prev(int h) const { return h%3==0 ? h+2 : h-1; }
	inline int opp(int h) const { return _opp[h]; }
	inline int next_around_vertex(int h) const { return opp(prev(h)); }
	inline int prev_around_vertex(int h) const { return opp(h) == -1 ? -1 : next(opp(h)); }

	// Half-edge geometry
	inline int from(int h) const { return h2v[h]; }
	inline int to(int h) const { return h2v[next(h)]; }
	inline vec3 geom(int h) const { return points[to(h)] - points[from(h)]; } 

	// One ring iterator
	struct OneRing {
		const Mesh *m;
		const int h0;
		OneRing(const Mesh *m, int v): m(m), h0(m->_v2h[v]) {}
		struct iterator {
			const Mesh *m;
			int h;
			inline iterator(): m(nullptr), h(-1) {}
			inline iterator(const Mesh *m, int h): m(m), h(h) {}
			inline iterator& operator++() { h = m->next_around_vertex(h); if(h != -1 && h == m->v2h(m->h2v[h])) h = -1; return *this; }
			inline bool operator!=(const iterator& rhs) const { return h != rhs.h; }
			inline const int &operator*() const { return h; }
		};
		inline iterator begin() const { return {m, h0}; }
		inline iterator end  () const { return {}; }
	};
	inline OneRing one_ring(int v) const { return OneRing(this, v); }
	inline bool is_boundary_vert(int v) const { return opp(v2h(v)) == -1; }

	// Convertion
	inline int f(int h) const { return h/3; }
	inline int h(int f, int i) const { return 3*f+i; }
	inline int v(int f, int i) const { return h2v[h(f, i)]; }
	inline Range hs(int f) const { return Range(h(f, 0), h(f, 3)); }

	void compute_opp();

	// DANGER!!! BE CAREFULL!!!
	inline std::vector<int>& getOpp() { return _opp; }
	inline std::vector<int>& getV2H() { return _v2h; }

private:
	std::vector<int> _opp, _v2h;
};
template<> struct std::iterator_traits<Mesh::OneRing::iterator> {
	using iterator_category = std::forward_iterator_tag;
};

std::tuple<Mesh, std::vector<vec2>, std::vector<bool>> readObj(const std::string filename);

void writeObj(const std::string filename, const Mesh &m, const std::vector<vec2> &uv, const std::vector<bool> &feature={});

inline double uvArea(const Mesh &m, const std::vector<vec2> &uv, int f) {
	const vec2 u = uv[m.h(f, 1)] - uv[m.h(f, 0)];
	const vec2 v = uv[m.h(f, 2)] - uv[m.h(f, 0)];
	return 0.5 * (u.x*v.y - u.y*v.x);
}
