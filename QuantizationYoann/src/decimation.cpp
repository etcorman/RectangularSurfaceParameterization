// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#include "decimation.h"

#include <algorithm>
#include <numeric>
#include <queue>
#include <set>
#include <cassert>

using namespace std;

void makeDelaunay(Mesh &m, CutGraph &cg, vector<vec2> &uv, vector<bool> &feature, SparseMatrix &D) {
	// Create a queue containing all haf-edges
	vector<int> Q(m.ncorners()); iota(Q.begin(), Q.end(), 0);
	vector<bool> inQ(m.ncorners(), true);
	// Process the queue
	while(!Q.empty()) {
		const int h = Q.back(); Q.pop_back();
		inQ[h] = false;

		// We need to keep feature edges
		if(feature[h]) continue; // can't flip feature

		/*
				b					 		     b
		 nh/       |\ph2                 ph2/        |\nh2
		  /_  /||    \                     /_   h2 \   \
		c	  h||h2  _ d       ===>      c     =====   _ d
		 \	   ||/   /                     \   \ h     /
		ph\|        /nh2			      nh\|        /ph
			   a                                 a
		*/

		// Early check
		const int h2 = m.opp(h);
		if(h2 == -1) continue; // Can't flip border
		const int ph = m.prev(h), ph2 = m.prev(h2);
		const int c = m.from(ph), d = m.from(ph2);
		if(c == d) continue; // Can't create flipped edge
		const auto existsHalfEdge = [&](int a, int b)->bool {
			return find_if(m.one_ring(a).begin(), m.one_ring(a).end(), [&](int h) { return m.to(h) == b; }) != m.one_ring(a).end();
		};
		if(existsHalfEdge(c, d) || existsHalfEdge(d, c)) continue; // Flipped edge already exist

		// Check Delaunay criterion
		const auto getOppAngle = [&](int h)->double {
			vec2 u = uv[h];
			h = m.next(h);
			vec2 v = uv[h];
			h = m.next(h);
			u -= uv[h];
			v -= uv[h];
			return acos((u * v) / (u.norm() * v.norm()));
		};
		if(getOppAngle(h) + getOppAngle(h2) <= M_PI) continue;

		const int nh = m.next(h);
		const int nh2 = m.next(h2);

		// create a 0 jump between the two facets
		const int j = cg.jump[h];
		cg.jump[h] = cg.jump[h2] = 0;
		for(const int e : {nh2, ph2}) {
			const int oe = m.opp(e);
			if(oe == -1) continue;
			cg.jump[e ] = (cg.jump[e] + j) % 4;
			cg.jump[oe] = (4 - cg.jump[e]) % 4;
		}

		// Store new uvs
		const int a = m.from(h), b = m.from(h2);
		const auto getVec2Row = [&](int h) {
			if(D.empty()) return vec2T<SparseMatrix::Row> {{},{}}; 
			return vec2T<SparseMatrix::Row> { D.getRow(2*h), D.getRow(2*h+1) };
		};
		const vec2 aUV = uv[h];  const auto aMat = getVec2Row(h);
		const vec2 bUV = uv[nh]; const auto bMat = getVec2Row(nh);
		const vec2 cUV = uv[ph]; const auto cMat = getVec2Row(ph);
		const vec2 dUV = (uv[ph2] - uv[h2]).rotate90(j) + bUV;
		const auto dMat = (getVec2Row(ph2) - getVec2Row(h2)).rotate90(j) + bMat;

		// Flip and update uvs and jumps
		#define FLIP_UPDATE(_h_, _v_) m.h2v[_h_] = _v_; uv[_h_] = _v_##UV; if(!D.empty()) { D.getRow(2*_h_) = _v_##Mat.x; D.getRow(2*_h_+1) = _v_##Mat.y; }
		FLIP_UPDATE(h, d);
		FLIP_UPDATE(nh, c);
		FLIP_UPDATE(ph, a);
		FLIP_UPDATE(h2, c);
		FLIP_UPDATE(nh2, d);
		FLIP_UPDATE(ph2, b);
		#define FLIP_UPDATE_2(_h1_, _h2_) swap(cg.jump[_h1_],  cg.jump[_h2_]); swap(m.getOpp()[_h1_], m.getOpp()[_h2_]); vector<bool>::swap(feature[_h1_], feature[_h2_])
		FLIP_UPDATE_2(ph, nh);
		FLIP_UPDATE_2(ph, ph2);
		FLIP_UPDATE_2(ph, nh2);
		for(const int e : {ph, ph2, nh, nh2}) if(m.opp(e) != -1) m.getOpp()[m.opp(e)] = e;
		if(m.v2h(a) == h || m.v2h(a) == nh2) m.getV2H()[a] = ph;
		if(m.v2h(b) == nh || m.v2h(b) == h2) m.getV2H()[b] = ph2;
		if(m.v2h(c) == ph)  m.getV2H()[c] = nh;
		if(m.v2h(d) == ph2) m.getV2H()[d] = nh2;
	
		// Fill Q
		for(const int e : {ph, ph2, nh, nh2}) if(!inQ[e]) {
			inQ[e] = true;
			Q.push_back(e);
		}
	}
}

bool collapse(Mesh &m, CutGraph &cg, vector<vec2> &uv, vector<bool> &feature, SparseMatrix &D) {
	vector<bool> to_kill_v(m.nverts(), false);
	vector<bool> to_kill_f(m.nfacets(), false);
	vector<double> score(m.ncorners(), 0);
	vector<double> minArea(m.nverts(), 0);

	const auto getVec2Row = [&](int h) { return vec2T<SparseMatrix::Row> { D.getRow(2*h), D.getRow(2*h+1) }; };
	const auto setDrow = [&](int h, const vec2T<SparseMatrix::Row> &r) { D.getRow(2*h) = r.x; D.getRow(2*h+1) = r.y; };
	const auto cpyDrow = [&](int h, int h2) { if(D.empty()) return; D.getRow(2*h) = D.getRow(2*h2); D.getRow(2*h+1) = D.getRow(2*h2+1); };

	// Compute min area per vertex
	constexpr double INF = std::numeric_limits<double>::infinity();
	const auto getMinArea = [&](int v)->double {
		if(cg.has2BeInt(uv, feature, v)) return INF;
		double max_a = std::numeric_limits<double>::max();
		for(const int h : m.one_ring(v)) {
			const int nh = m.next(h), ph = m.prev(h);
			const double a = (uv[nh].x - uv[h].x) * (uv[ph].y - uv[h].y) - (uv[nh].y - uv[h].y) * (uv[ph].x - uv[h].x);
			max_a = min(max_a, a);
		}
		return max_a;
	};

	// Make zero jump func (assume a is not a singularity)
	const auto makeZJ = [&](int v, bool updateD)->void {
		for(const int h : m.one_ring(v)) if(h != m.v2h(v)) {
			const int h2 = m.opp(h);
			assert(h2 != -1);
			const int j = cg.jump[h2];
			cg.jump[h] = cg.jump[h2] = 0;
			const int nh = m.next(h), ph = m.prev(h);
			for(const int e : {nh, ph}) {
				const int oe = m.opp(e);
				if(oe == -1) continue;
				cg.jump[e ] = (cg.jump[e] + j) % 4;
				cg.jump[oe] = (4 - cg.jump[e]) % 4;
			}
			const int nh2 = m.next(h2);
			uv[ph] = (uv[ph] -= uv[nh]).rotate90(j) += uv[h2];
			uv[nh] = uv[h2];
			uv[h] = uv[nh2];
			if(updateD) {
				setDrow(ph, (getVec2Row(ph) -= getVec2Row(nh)).rotate90(j) += getVec2Row(h2));
				cpyDrow(nh, h2);
				cpyDrow(h, nh2);
			}
		}
	};
	
	// Score func (assume from(h) is not a singularity)
	const auto getScore = [&](int h)->double {
		const int v = m.from(h);
		if(m.is_boundary_vert(v) && m.opp(h) != -1) return INF; // Can't collapse a border vertex to an interior vertex
		const auto willBeIsolated = [&](int h) { return m.next_around_vertex(h) == -1 && m.prev_around_vertex(h) == -1; };
		if(willBeIsolated(m.prev(h))) return INF; // This collapse will leave an isolated vertex
		const int ph = m.prev_around_vertex(h);
		if(ph != -1 && willBeIsolated(m.next(ph))) return INF; // This collapse will leave an isolated vertex
		// check feature
		int countFeature = m.is_boundary_vert(v);
		for(const int e : m.one_ring(v)) countFeature += feature[e];
		if(countFeature > 0 && (countFeature != 2 || !feature[h])) return INF;
		// =============
		vector<tuple<int, vec2, int>> old;
		for(const int e0 : m.one_ring(v)) for(const int e : {e0, m.next(e0), m.prev(e0)}) old.emplace_back(e, uv[e], cg.jump[e]);
		makeZJ(v, false);
		const vec2 uvb = uv[m.next(h)];
		double nma = std::numeric_limits<double>::max();
		for(const int e : m.one_ring(v)) if(e != h && e != ph) {
			const int ne = m.next(e), pe = m.prev(e);
			nma = std::min(nma, (uv[ne].x - uvb.x) * (uv[pe].y - uvb.y) - (uv[ne].y - uvb.y) * (uv[pe].x - uvb.x));
		}
		for(const auto &[e, u, j] : old) { uv[e] = u; cg.jump[e] = j; if(const int oe = m.opp(e); oe != -1) cg.jump[oe] = (4-j)%4; }
		// The new ring can't be empty and new triangles can't be inverted
		if(nma == std::numeric_limits<double>::max() || (nma < 1.e-2 && minArea[v] > nma)) return INF;
		// Highest priotity to small edges
		const double n2 = (uv[m.next(h)] - uv[h]).norm2();
		if(n2 < 1.5) return n2;
		// Then to small triangles that can be improved
		if(nma > minArea[v]) return 2. + minArea[v] * (1. + minArea[v] / nma); // We can improve the min area
		// And finally to least bad
		return 1.e5 + minArea[v] / nma;
	};
	const double SCORE_EPS = 1.e-6;

	// Init queue
	typedef pair<double, int> QEl;
	vector<QEl> initQ(0);
	for(const int v : m.verts()) if(!isinf(minArea[v] = getMinArea(v)))
		for(const int h : m.one_ring(v)) if(!isinf(score[h] = getScore(h)))
			initQ.emplace_back(score[h], h);
	priority_queue<QEl, std::vector<QEl>, greater<QEl>> Q(greater<QEl>(), move(initQ));

	// Did any collapse occur
	bool answer = false;

	while(!Q.empty()) {
		const auto [sco, h] = Q.top(); Q.pop();
		if(to_kill_f[m.f(h)] || sco != score[h]) continue;
		if(const double sco2 = getScore(h); std::abs(sco - sco2) > SCORE_EPS) {
			if(!isinf(sco2)) Q.emplace(score[h] = sco2, h);
			continue;
		}
		score[h] = INF;

		// Check the possibility to collapse
		const int nh = m.next(h);
		const int h2 = m.opp(h);
		const int v = m.from(h), v2 = m.from(nh);

		// construct one ring
		set<int> ring;
		vector<int> incident;
		for(const int e : m.one_ring(v)) {
			ring.insert(m.to(e));
			ring.insert(m.from(m.prev(e)));
			incident.push_back(e);
		}
		ring.erase(v2);
		ring.erase(m.from(m.prev(h)));
		if(h2 != -1) ring.erase(m.from(m.prev(h2)));
		bool edgeAlreadyExists = false;
		for(const int e : m.one_ring(v2)) if(ring.count(m.to(e)) || ring.count(m.from(m.prev(e))))
			edgeAlreadyExists = true;
		if(edgeAlreadyExists) continue;
		answer = true;

		// Make zero jumps around vertex a
		makeZJ(v, !D.empty());

		// update jumps through deleted faces
		const int ph = m.prev(h);
		const int oph = m.opp(ph), onh = m.opp(nh);
		const auto reconnect = [&](const int in, const int out)->bool {
			if(in == -1) {
				if(out == -1) return true;
				cg.jump[out] = 0;
				m.getOpp()[out] = -1;
				feature[out] = true;
			} else if(out == -1) {
				cg.jump[in] = 0;
				m.getOpp()[in] = -1;
				feature[in] = true;
			} else {
				cg.jump[in] = (4 - cg.jump[out]) % 4;
				m.getOpp()[in]  = out;
				m.getOpp()[out] = in;
				feature[in] = feature[out] = feature[in] || feature[out];
			}
			return false;
		};
		if(reconnect(oph, onh)) to_kill_v[m.from(ph)] = true;
		if(h2 != -1) {
			const int ph2 = m.prev(h2), nh2 = m.next(h2);
			const int oph2 = m.opp(ph2), onh2 = m.opp(nh2);
			if(reconnect(onh2, oph2)) to_kill_v[m.from(ph2)] = true;
		}

		// Make uv[a] = uv[b]
		for(const int e : incident) {
			uv[e] = uv[nh];
			cpyDrow(e, nh);
			m.h2v[e] = v2;
		}

		// Collapse
		if(m.v2h(v2) == h2 || m.v2h(v2) == nh) {
			assert(oph != -1);
			m.getV2H()[v2] = oph;
		}
		if(m.v2h(m.from(ph)) == ph) {
			if(oph == -1) m.getV2H()[m.from(ph)] = onh;
			else m.getV2H()[m.from(ph)] = m.next(oph);
		}
		if(h2 != -1) {
			const int ph2 = m.prev(h2), nh2 = m.next(h2);
			const int oph2 = m.opp(ph2), onh2 = m.opp(nh2);
			if(m.v2h(m.from(ph2)) == ph2) {
				if(oph2 == -1) m.getV2H()[m.from(ph2)] = onh2;
				else m.getV2H()[m.from(ph2)] = m.next(oph2);
			}
		}
		to_kill_v[v] = true;
		to_kill_f[m.f(h)] = true;
		if(h2 != -1) to_kill_f[m.f(h2)] = true;

		// Update Q
		ring.insert(v2);
		ring.insert(m.from(m.prev(h)));
		if(h2 != -1) ring.insert(m.from(m.prev(h2)));
		for(const int u : ring) if(isinf(minArea[u]) || to_kill_v[u]) {
			for(const int e : m.one_ring(u)) score[e] = INF;
		} else {
			minArea[u] = getMinArea(u);
			for(const int e : m.one_ring(u)) {
				const double sco = getScore(e);
				if(abs(sco - score[e]) > SCORE_EPS && !isinf(score[e] = sco))
					Q.emplace(sco, e);
			}
		}
	}

	if(!answer) return false;

	// Delete vertices and faces
	int nv = 0, nf = 0;
	vector<int> old2new(m.nverts());
	for(const int v : m.verts()) if(!to_kill_v[v]) {
		m.points[nv] = m.points[v];
		cg.valence[nv] = cg.valence[v];
		old2new[v] = nv++; 
	}
	for(const int f : m.facets()) if(!to_kill_f[f]) {
		for(const int i : Range(3)) {
			const int h = m.h(f, i), nh = m.h(nf, i);
			m.h2v[nh] = old2new[m.h2v[h]];
			cg.jump[nh] = cg.jump[h];
			uv[nh] = uv[h];
			cpyDrow(nh, h);
			feature[nh] = feature[h];
		}
		++ nf;
	}
	const int nh = 3*nf;
	m.points.resize(nv);
	m.h2v.resize(nh);
	cg.jump.resize(nh);
	cg.valence.resize(nv);
	uv.resize(nh);
	D.resize(2*nh);
	feature.resize(nh);
	m.compute_opp();

	return true;
}
