// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#include "cut_graph.h"

#include <queue>

using namespace std;

CutGraph compute_cut_graph(const Mesh &m, vector<vec2> &uv) {
	CutGraph cg(m);
	vector<bool> seen(m.nfacets(), false);
	vector<int> Qring;
	vector<int> nz(m.nverts(), 0); // nz[v] is the number of non cut graph half-edges starting at v
	for(int v : m.h2v) ++ nz[v];
	const auto setZero = [&](int h)->void {
		cg.isCut[h] = false;
		cg.jump[h] = 0;
		const int v = m.from(h);
		if(--nz[v] == 1) Qring.push_back(v);
	};
	for(const int f_seed : m.facets()) if(!seen[f_seed]) {
		queue<int> Qtree;
		seen[f_seed] = true;
		Qtree.push(f_seed);
		while(!Qtree.empty()) {
			const int f = Qtree.front(); Qtree.pop();
			for(const int h : m.hs(f)) {
				const int h2 = m.opp(h);
				if(h2 == -1) {
					setZero(h);
					continue;
				}
				const int f2 = m.f(h2);
				const vec2 a = uv[m.next(h )] - uv[h ];
				const vec2 b = uv[m.next(h2)] - uv[h2];
				double best_dist = numeric_limits<double>::max();
				for(const int j : Range(4)) {
					const double dist = (a + vec2(b).rotate90(j)).norm2();
					if(dist < best_dist) {
						best_dist = dist;
						cg.jump[h] = j;
					}
				}
				if(!seen[f2]) {
					for(const int i : Range(3)) uv[m.h(f2, i)].rotate90(cg.jump[h]);
					setZero(h);
					setZero(h2);
					seen[f2] = true;
					Qtree.push(f2);
				}
			}
		}
	}
	// Compute valences
	for(const int v : m.verts()) {
		double angle = 0.;
		for(const int h : m.one_ring(v)) {
			const vec2 a = uv[m.next(h)] - uv[h];
			const vec2 b = uv[m.prev(h)] - uv[h];
			angle += acos((a*b) / (a.norm()*b.norm()));
		}
		cg.valence[v] = round(angle / M_PI_2);
	}
	// Remove some edges from cut graph
	while(!Qring.empty()) {
		const int v = Qring.back(); Qring.pop_back();
		if(cg.valence[v] != 4 || m.is_boundary_vert(v)) continue;
		int h = m.v2h(v);
		while(!cg.isCut[h]) h = m.next_around_vertex(h);
		setZero(h);
		setZero(m.opp(h));
	}
	return cg;
}

vec2b CutGraph::getHas2BeInt(const vector<vec2> &uv, const std::vector<bool> &feature, int v) const {
	if(isSingularity(v)) return {true, true};
	vec2b c(false, false);
	int dim =  0;
	for(const int h : m.one_ring(v)) {
		for(const int h : {h, m.prev(h)}) if(feature[h]) {
			const vec2 a = uv[m.next(h)] - uv[h];
			c[dim ^ int(abs(a.y) < abs(a.x))] = true;
		}
		dim ^= jump[m.prev(h)]&1;
	}
	return c;
}
