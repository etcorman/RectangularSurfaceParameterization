// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#include "refine/imprint.h"

#include <array>
#include <numeric>
#include <map>
#include <algorithm>

using namespace std;

using Triangle = array<vec2, 3>;

void clip_triangle_by_half_space(const Triangle &tri, const vec2 &O, const vec2 &n, vector<Triangle> &clipped) {
	constexpr double EPS = 1e-6, EPS_b = 1e-10;
	bool in[3];
	for(const int i : Range(3)) in[i] = ((O - tri[i]) * n > EPS);

	if(!in[0] && !in[1] && !in[2]) return;
	if(in[0] && in[1] && in[2])	return clipped.push_back(tri);

	const int i = in[0] == in[1] ? 2 : (in[0] == in[2] ? 1 : 0);

	const double dO0 = (O - tri[i]) * n;
	const double c1 = clamp(dO0 / ((tri[(i+1)%3] - tri[i]) * n), 0., 1.);
	const double c2 = clamp(dO0 / ((tri[(i+2)%3] - tri[i]) * n), 0., 1.);
	const vec2 P1 = (1. - c1) * tri[i] + c1 * tri[(i+1)%3];
	const vec2 P2 = (1. - c2) * tri[i] + c2 * tri[(i+2)%3];

	if(in[i]) {
		if(c1 > EPS_b && c2 > EPS_b) clipped.push_back({tri[i],P1,P2});
		return;
	}
	if(c1 < 1.-EPS_b) clipped.push_back({tri[(i+2)%3],P1,tri[(i+1)%3]});
	if(c2 < 1.-EPS_b && (c1 > EPS_b || c2 > EPS_b)) clipped.push_back({P1,tri[(i+2)%3],P2});
}

vector<Triangle> clip_triangle_by_triangle(const Triangle& tri0, const Triangle& tri1) {
	vector<Triangle> clipped[2];
	clipped[1] = { tri0 };
	for(const int i : Range(3)) {
		const vec2 P = tri1[i];
		const vec2 n = (tri1[(i + 1) % 3] - P).rotate90(3);
		swap(clipped[0], clipped[1]);
		clipped[1].clear();
		for(const Triangle &tr : clipped[0]) clip_triangle_by_half_space(tr, P, n, clipped[1]);
	}
	return clipped[1];
}

vec3 bary_coord(const Triangle &tri, const vec2 &P) {
	vec3 C;
	for(const int e : Range(3)) {
		const vec2 u = tri[(e+1)%3] - P;
		const vec2 v = tri[(e+2)%3] - P;
		C[e] = max(1e-30, u.x*v.y - u.y*v.x);
	}
	return C /= (C[0]+C[1]+C[2]);
}

vector<int> gather_all_overlapping_triangles(const Mesh &m, int f0, vector<vec2>& uv, const Triangle& tri) {
	static vector<bool> visited;
	if((int) visited.size() < m.nfacets()) visited.resize(m.nfacets(), false);
	vector<int> faces;
	const auto add = [&](int f) {
		faces.push_back(f);
		visited[f] = true;
	};
	add(f0);
	for(int ind = 0; ind < (int) faces.size(); ++ind) {
		const int f = faces[ind];
		const vector<Triangle> ts = clip_triangle_by_triangle({ uv[m.h(f, 0)], uv[m.h(f, 1)], uv[m.h(f, 2)] }, tri);
		for(const int h : m.hs(f)) {
			const int opp = m.opp(h);
			if(opp == -1) continue;
			if(visited[m.f(opp)]) continue;
			const vec2 ref = uv[m.next(h)] - uv[h];
			const vec2 nref = vec2(ref).normalize().rotate90(1);
			bool parallel = false;
			for(const Triangle &t : ts) {
				int on_edge = 0;
				for(const vec2 &p : t) if(abs((p-uv[h])*nref) < 1e-6) ++on_edge;
				if(on_edge >= 2) {
					parallel = true;
					break;
				}
			}
			if(!parallel) continue;
			vec2 ref2 = uv[m.next(opp)] - uv[opp];
			int best_rot = 0;
			double best_dist = numeric_limits<double>::max();
			for(const int rot : Range(4)) {
				const double dist = (ref + ref2).norm2();
				if(dist < best_dist) {
					best_rot = rot;
					best_dist = dist;
				}
				ref2.rotate90(1);
			}
			for(const int e : {opp, m.next(opp), m.prev(opp)}) uv[e].rotate90(best_rot);
			const vec2 t = .5 * ((uv[h]+uv[m.next(h)]) - (uv[opp]+uv[m.next(opp)]));
			for(const int e : {opp, m.next(opp), m.prev(opp)}) uv[e] += t;
			add(m.f(opp));
		}
	}
	for(const int f : faces) visited[f] = false;
	return faces;
}

void reconnect(Mesh &m, vector<vec2> &uv) {
	const double EPS = 5.e-6, EPS2 = EPS*EPS;

	//====================
	//== MERGE VERTICES ==
	//====================

	vector<int> g(m.nverts(), -1);
	const auto comp = [](const vec3 &a, const vec3 &b) {
		return make_tuple(a.x, a.y, a.z) < make_tuple(b.x, b.y, b.z);
	};
	map<vec3, int, decltype(comp)> ps(comp);
	const auto merge = [&](int v, auto it) {
		const int u = it->second;
		auto node = ps.extract(it);
		node.key() = (m.points[v] - g[u] * node.key()) / (1 - g[u]);
		auto nit = it; ++nit;
		ps.insert(nit, move(node));
		g[v] = u;
		-- g[u];
	};
	for(const int v : m.verts()) {
		const vec3 &p = m.points[v];
		const auto it0 = ps.lower_bound(p);
		for(auto it = it0; it != ps.end() && it->first.x < p.x+EPS; ++it)
			if((p - it->first).norm2() < EPS2) {
				merge(v, it);
				break;
			}
		if(g[v] != -1) continue;
		for(auto it = it0; it != ps.begin() && (--it)->first.x > p.x-EPS;)
			if((p - it->first).norm2() < EPS2) {
				merge(v, it);
				break;
			}
		if(g[v] != -1) continue;
		ps.emplace_hint(it0, p, v);
	}
	int nv = 0;
	for(const int v : m.verts()) {
		if(g[v] < 0) g[v] = nv++;
		else g[v] = g[g[v]];
		m.h2v[v] = g[v];
	}
	for(const auto [p, v] : ps) m.points[g[v]] = p;
	m.points.resize(nv);

	//===============================
	//== REMOVE DEGENERATED FACETS ==
	//===============================

	for(int f = 0; f < m.nfacets(); ++f) {
		if(m.v(f, 0) == m.v(f, 1) || m.v(f, 0) == m.v(f, 2) || m.v(f, 1) == m.v(f, 2)) {
			for(const int i : Range(3)) {
				m.h2v[3*f+i] = m.h2v[m.ncorners()-3+i];
				uv[3*f+i] = uv[m.ncorners()-3+i];
			}
			m.h2v.resize(m.h2v.size()-3);
			uv.resize(uv.size()-3);
			--f;
		}
	}

	//=================================
	//== ADD VERTICES FOR CONNECTION ==
	//=================================

	vector<int> hs;
	for(int step = 0; step < 4; ++step) {
		hs.resize(m.ncorners());
		iota(hs.begin(), hs.end(), 0);
		sort(hs.begin(), hs.end(), [&](int h, int h2) { return m.h2v[h] < m.h2v[h2]; });
		for(int i = 0; i < (int)hs.size();) {
			const int i0 = i;
			while(i < (int)hs.size() && m.h2v[hs[i]] == m.h2v[hs[i0]]) ++i;
			const int i1 = i;
			for(i = i0; i < i1; ++i) {
				const int h = hs[i];
				const int v = m.from(h);
				const int v2 = m.to(h);
				bool has_opposite = false;
				for(int j = i0; j < i1; ++j)
					has_opposite |= m.from(m.prev(hs[j])) == v2;
				if(has_opposite) continue;
				vec3 eh = m.geom(h);
				const double neh = eh.norm();
				eh /= neh;
				for(int j = i0; j < i1; ++j) if(j != i) {
					const int h2 = m.prev(hs[j]);
					vec3 eh2 = m.geom(h2);
					const double neh2 = eh2.norm();
					eh2 /= neh2;
					if(eh*eh2 > EPS2 - 1.) continue;
					if(neh < neh2) {
						m.points[v2] = .5*(m.points[v2] + m.points[v] - neh*eh2);
						const double t = neh/neh2;
						m.h2v.push_back(m.h2v[h2]); uv.push_back(uv[h2]);
						m.h2v[h2] = v2;
						uv[h2] = (1.-t)*uv[m.next(h2)] + t*uv[h2];
						m.h2v.push_back(v2); uv.push_back(uv[h2]);
						m.h2v.push_back(m.h2v[m.prev(h2)]); uv.push_back(uv[m.prev(h2)]);
					} else {
						const int v3 = m.h2v[h2];
						m.points[v3] = .5*(m.points[v3] + m.points[v] + neh2*eh);
						const double t = neh2/neh;
						m.h2v.push_back(v); uv.push_back(uv[h]);
						m.h2v[h] = v3;
						uv[h] = (1.-t)*uv[h] + t*uv[m.next(h)];
						m.h2v.push_back(v3); uv.push_back(uv[h]);
						m.h2v.push_back(m.h2v[m.prev(h)]); uv.push_back(uv[m.prev(h)]);
					}
					break;
				}
			}
		}
	}
}

pair<Mesh, vector<vec2>> transfert_coarse_to_fine(const Mesh &fine, const vector<vec2> &f_uv,
													const Mesh &coarse, const vector<vec2> &c_uv, const vector<vec2> &c_uv2) {
	vector<int> c2f(coarse.nverts(), 0);
	for(const int vc : coarse.verts()) 
		for(const int vf : fine.verts()) 
			if((fine.points[vf] - coarse.points[vc]).norm2() < (fine.points[c2f[vc]] - coarse.points[vc]).norm2())
				c2f[vc] = vf;

	Mesh m;
	vector<vec2> uv, fuv = f_uv;
	for(const int fc : coarse.facets()) {
		const Triangle tri = { c_uv[coarse.h(fc, 0)], c_uv[coarse.h(fc, 1)], c_uv[coarse.h(fc, 2)] };
		bool done_with_fc = false;
		for(const int h : fine.one_ring(c2f[coarse.v(fc, 0)])) {
			for([[maybe_unused]] const int rot : Range(4)) {
				// apply rot and translation
				for(const int e : {h, fine.next(h), fine.prev(h)}) fuv[e].rotate90(1); 
				const vec2 t = tri[0] - fuv[h];
				for(const int e : {h, fine.next(h), fine.prev(h)}) fuv[e] += t;
			
				vector<int> subset = gather_all_overlapping_triangles(fine, fine.f(h), fuv, tri);

				if(find_if_not(Range(3).begin(), Range(3).end(), [&](const int i) {
					return find_if(subset.begin(), subset.end(), [&](const int f) {
						return find_if(m.hs(f).begin(), m.hs(f).end(), [&](const int e) {
							return c2f[coarse.v(fc, i)] == fine.from(e) && (tri[i] - fuv[e]).norm2() < 1e-10;
						}) != m.hs(f).end();
					}) != subset.end();
				}) != Range(3).end()) continue;

				for(const int f : subset) {
					const Triangle local_tri = { fuv[fine.h(f, 0)], fuv[fine.h(f, 1)], fuv[fine.h(f, 2)] };
					for(const Triangle t : clip_triangle_by_triangle(local_tri, tri))
						for(const int i : Range(3)) {
							const vec3 c = bary_coord(local_tri, t[i]);
							m.points.push_back(
								  c[0] * fine.points[fine.from(fine.h(f, 0))]
								+ c[1] * fine.points[fine.from(fine.h(f, 1))]
								+ c[2] * fine.points[fine.from(fine.h(f, 2))]
							);
							const vec3 c2 = bary_coord(tri, t[i]);
							uv.push_back(
								  c2[0] * c_uv2[coarse.h(fc, 0)]
								+ c2[1] * c_uv2[coarse.h(fc, 1)]
								+ c2[2] * c_uv2[coarse.h(fc, 2)]
							);
						}
				}
				done_with_fc = true;
				break;
			}
			if(done_with_fc) break;
		}
		if(!done_with_fc) cerr << "[Inprinting]: Failed to inprint coarse facet " << fc << endl;
	}

	m.h2v.resize(m.nverts());
	iota(m.h2v.begin(), m.h2v.end(), 0);
	reconnect(m, uv);
	return make_pair(move(m), move(uv));
}
