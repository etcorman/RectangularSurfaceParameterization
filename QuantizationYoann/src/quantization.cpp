// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#include "quantization.h"

#include <algorithm>
#include <random>
#include <queue>
#include <cassert>
#include <numeric>
#include <functional>

using namespace std;
typedef long long ll;

vector<vec2i> quantize(const Mesh &m, const CutGraph &cg, const vector<vec2> &uv, const vector<bool> &feature, bool coarsest) {
	vector<vec2i> quv(m.ncorners());
	mt19937 re;
	
	// Compute delta uv
	vector<vec2> duv(m.ncorners());
	for(const int h : m.corners()) duv[h] = uv[m.next(h)] - uv[h];
	
	// Has2BeNull
	vector<vec2b> has2BeNull(m.ncorners(), {false, false});
	for(const int h : m.corners()) if(feature[h])
		has2BeNull[h][int(abs(duv[h].y) < abs(duv[h].x))] = true;

	// DATA
	vector<int> bfsIdx(m.ncorners() << 2, 0);
	vector<int> pred(bfsIdx.size());
	vector<vec3> weight(8 * m.ncorners());
	vector<vec3> bfsDist(4 * m.ncorners());
	int bfsID = 0u;

	// *****************
	// Usefull functions
	// *****************

	const auto getDet = [&](int h)->ll {
		const int nh = m.next(h);
		return (ll) quv[h].x * quv[nh].y - (ll) quv[h].y * quv[nh].x;
	};

	const auto getOp = [&](int oe, int d)->int {
		return (oe << 2) + (((d + cg.jump[oe]) & 0b11) ^ 0b10);
	};

	const vec2i addV[4] {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
	const auto add = [&](int h, int d)->void { quv[h][d&1] += 1 - (d&2); };

	// const auto getEfacet = [&](int f) {
	// 	const int h1 = m.h(f, 0), h2 = m.h(f, 1);
	// 	const double detM = duv[h1].x * duv[h2].y - duv[h1].y * duv[h2].x;
	// 	vec2 Jx = duv[h2].y * (vec2)quv[h1] - duv[h1].y * (vec2)quv[h2];
	// 	vec2 Jy = duv[h1].x * (vec2)quv[h2] - duv[h2].x * (vec2)quv[h1];
	// 	if(!coarsest) {
	// 		Jx.x -= detM;
	// 		Jy.y -= detM;
	// 	}
	// 	return (Jx.norm2() + Jy.norm2()) / detM;
	// };

	const auto getE = [&]()->double {
		double E = 0.;
		for(const int h : m.corners()) {
			const double dx = coarsest ? quv[h].x : quv[h].x - duv[h].x;
			const double dy = coarsest ? quv[h].y : quv[h].y - duv[h].y;
			const double duv_norm = max(2.e-4, duv[h].x*duv[h].x + duv[h].y*duv[h].y);
			E += (dx*dx + dy*dy) / duv_norm;
		}
		// for(const int f : m.facets()) E += getEfacet(f);
		return E;
	};

	const auto computeW = [&](int h0)->void {
		// const double E0 = getEfacet(m.f(h0));
		for(const int h : {h0, m.next(h0), m.prev(h0)}) {
			// Check the update of det is valid
			for(const auto [h2, i] : { vec2i{m.next(h), 0}, vec2i{m.prev(h), 1} }) for(const int d : Range(4)) {
				const vec2 diff = coarsest ? vec2(quv[h2]) : quv[h2] - duv[h2];
				const double dd = diff*addV[d];
				weight[((4*h+d)<<1)+i] = dd > .5 ? vec3{0., 0., 1. / (.5 + dd)} : vec3{0., (1.5 - dd), 0.};
				add(h, d);
				add(h2, d^2);
				// double E = getEfacet(m.f(h));
				// if(E <= E0) weight[(4*h+d)<<1] = {0., 0., E/E0};
				// else weight[(4*h+d)<<1] = {0., E/E0, 0.};
				if(getDet(h) <= 0) ++ weight[((4*h+d)<<1)+i].x;
				add(h2, d);
				add(h, d^2);
			}
		}
		// POSSIBLE IMPROVEMENT: - for non-zero edge
	};

	const auto bfs = [&](int ed0)->int {
		++ bfsID;
		typedef pair<vec3, int> QEl;
		const auto Comp = [](const vec3 &a, const vec3 &b)->bool {
			return make_tuple(a.x, a.y, a.z) < make_tuple(b.x, b.y, b.z);
		};
		const auto Comp2 = [&](const QEl &a, const QEl &b)->bool { return !Comp(a.first, b.first); };
		priority_queue<QEl, vector<QEl>, decltype(Comp2)> Q(Comp2);
		vec3 border_dist = { numeric_limits<double>::max(), 0., 0. };
		int border_opp = -1; // half-edge across which we leave the mesh
		Q.emplace(bfsDist[ed0] = {0., 0., 0.}, ed0);
		while(!Q.empty()) {
			auto [dist, ed] = Q.top(); Q.pop();
			if(ed == -1) return border_opp;
			if(dist != bfsDist[ed]) continue;
			if(ed == ed0 && bfsIdx[ed] == bfsID) return getOp(m.opp(ed>>2), ed&0b11);
			const int e = ed >> 2, d2 = (ed & 0b11)^0b10;
			for(const auto [e2, i] : { vec2i{m.next(e), 0}, vec2i{m.prev(e), 1} }) {
				if(has2BeNull[e2][d2&0b01]) continue;
				const int ed2 = (e2 << 2) + d2;
				const vec3 new_dist = dist + weight[(ed<<1)+i];
				const int oe2 = m.opp(e2);
				if(oe2 == -1) {
					if(!Comp(new_dist, border_dist)) continue; // dist not lower
					pred[ed2] = ed;
					border_dist = new_dist;
					border_opp = ed2;
					Q.emplace(new_dist, -1);
				} else {
					const int oed2 = getOp(oe2, d2);
					if(bfsIdx[oed2] == bfsID && !Comp(new_dist, bfsDist[oed2])) continue; // dist not lower
					pred[ed2] = ed;
					bfsIdx[oed2] = bfsID;
					bfsDist[oed2] = new_dist;
					Q.emplace(new_dist, oed2);
				}
			}
		}
		return -1;
	};

	const auto getPath = [&](int ed0, int ed1, vector<int> &path)->void {
		while(true) {
			path.push_back(ed1);
			path.push_back(ed1 = pred[ed1]);
			if(ed1 == ed0) return;
			ed1 = getOp(m.opp(ed1>>2), ed1&0b11);
		}
	};

	const auto find_cycle = [&](int h, int d)->vector<int> {
		int ed0 = (h << 2) + d;
		const int ed0b = m.opp(h) == -1 ? -1 : getOp(m.opp(h), d);
		vector<int> path;
		int ed1 = bfs(ed0);
		if(ed1 == -1) return {};
		if(ed0b != -1 && ed1 != ed0b) { // End on border
			getPath(ed0, ed1, path);
			for(const int ed : path) {
				add(ed>>2, ed&0b11);
				computeW(ed>>2);
			}
			ed1 = bfs(ed0b);
			for(const int ed : path) {
				add(ed>>2, (ed&0b11)^0b10);
				computeW(ed>>2);
			}
			if(ed1 == -1) return {};
			if(ed1 == ed0) path.clear();
			ed0 = ed0b;
		}
		getPath(ed0, ed1, path);
		return path;
	};

	const auto applyPath = [&](const vector<int> &path, double &E)->bool {
		if(path.empty()) return false;
		for(const int ed : path) add(ed>>2, ed&3);
		for(const int ed0 : path) if(getDet(ed0>>2) <= 0) {
			for(const int ed : path) add(ed>>2, (ed&3)^2);
			return false;
		}
		const double E2 = getE();
		if(E2 < E) {
			for(const int ed : path) computeW(ed>>2);
			E = E2;
			return true;
		}
		for(const int ed : path) add(ed>>2, (ed&3)^2);
		return false;
	};

	// Init quantization
	double scale0 = 1.;
	while(true) {
		for(const int h : m.corners()) {
			const int oh = m.opp(h);
			if(oh > h) continue;
			for(const int d : {0, 1}) quv[h][d] = round(scale0 * duv[h][d]);
			if(quv[h] == vec2i{0, 0}) {
				if(abs(duv[h][0]) > abs(duv[h][1])) quv[h][0] = duv[h][0] < 0. ? -1 : 1; 
				else quv[h][1] = duv[h][1] < 0. ? -1 : 1; 
			}
			if(oh != -1) quv[oh] = - vec2i(quv[h]).rotate90(cg.jump[oh]);
		}
		vector<vec2i> sumF(m.nfacets(), {0, 0});
		for(const int h : m.corners()) sumF[m.f(h)] += quv[h];
		for(const int d0 : {0, 1}) for(const int f0 : m.facets()) {
			while(sumF[f0][d0] != 0) {
				const int fd0 = (f0<<2) | (sumF[f0][d0] < 0 ? 2 : 0) | d0;
				using QEl = pair<int, int>;
				priority_queue<QEl, vector<QEl>, greater<QEl>> Q;
				vector<int> pred(4*m.nfacets()+1, -1);
				Q.emplace(0, fd0);
				pred[fd0] = -2;
				while(!Q.empty()) {
					const auto [dist, fd] = Q.top(); Q.pop();
					const int f = fd >> 2, d = fd & 0b11;
					if(f == m.nfacets() || (((fd0^fd) != 0b10 || abs(sumF[f0][d0]) > 1) && ((d&0b10) ? sumF[f][d&0b01] > 0 : sumF[f][d&0b01] < 0))) {
						for(int hd = pred[fd]; hd != -2;) {
							const int h = hd >> 2, d = hd & 0b11;
							quv[h] += addV[d];
							sumF[m.f(h)] += addV[d];
							const int oh = m.opp(h);
							if(oh != -1) {
								const int d2 = (d + cg.jump[oh]) & 0b11;
								quv[oh] -= addV[d2];
								sumF[m.f(oh)] -= addV[d2];
							}
							hd = pred[(m.f(h)<<2) | (d^0b10)];
						}
						break;
					}
					for(const int h : m.hs(f)) {
						if(has2BeNull[h][d&1]) continue;
						if((quv[h] + addV[d^2]) == vec2i{0, 0}) continue;
						const int h2 = m.opp(h);
						const int f2 = h2 == -1 ? m.nfacets() : m.f(h2);
						const int fd2 = (f2<<2) | (h2 == -1 ? 0 : ((d+cg.jump[h2]) & 0b11));
						if(pred[fd2] != -1) continue;
						pred[fd2] = (h<<2) | (d^0b10);
						if(pred[fd] != -2) quv[m.opp(pred[fd]>>2)] += addV[d];
						quv[h] -= addV[d];
						const int det2 = getDet(h);
						if(pred[fd] != -2) quv[m.opp(pred[fd]>>2)] -= addV[d];
						quv[h] += addV[d];
						const int det = getDet(h);
						if(det2 <= 0 && det2 < det) Q.emplace(dist + (1-det2)*1000, fd2);
						else if(det <= 0 && det < det2) Q.emplace(dist-100*(1-det), fd2);
						else Q.emplace(dist+1, fd2);
					}
				}
			}
		}
		bool ok = true;
		for(const int f : m.facets()) if(getDet(m.h(f, 0)) <= 0) {
			ok = false;
			break;
		}
		if(ok) break;
		scale0 *= sqrt(2.);
		if(scale0 > 1e3) {
			cerr << "Scaled snapping failed!!" << endl;
			exit(1);
		}
	}
	cerr << "Scale: " << scale0 << endl;

	// Check features
	for(const int h : m.corners()) for(const int d : Range(2)) if(has2BeNull[h][d]) assert(quv[h][d] == 0);

	// Compute weights and E
	for(const int f : m.facets()) computeW(m.h(f, 0));
	double E = getE();
	cerr << "Energy: " << E << endl;

	int countBetter0 = 0, countBetter = 1;
	const auto update = [&](int h) {
		for(const int d : Range(4)) if(!has2BeNull[h][d&1])
			if(applyPath(find_cycle(h, d), E)) {
				cerr << '\r' << E << "        ";
				if(((++ countBetter)%20) == 0) cerr << '\n';
				#ifndef NDEBUG
				cerr << flush;
				#endif
				return;
			}
	};
	// random order for edges
	vector<int> order(m.ncorners()); iota(order.begin(), order.end(), 0);
	shuffle(order.begin(), order.end(), re);
	while(countBetter0 < countBetter) {
		countBetter0 = countBetter;
		shuffle(order.begin(), order.end(), re);
		for(const int h : order) update(h);
	}
	if(countBetter%20 != 0) cerr << endl;

	// integrate quv
	vector<bool> seen(m.ncorners(), false);
	function<void(int, vec2i)> dfs;
	dfs = [&](int h, vec2i val)->void {
		seen[h] = true;
		vec2i add = quv[h];
		quv[h] = val;
		const int nh = m.next(h);
		if(!seen[nh]) dfs(nh, val + add);
		const int ph = m.prev_around_vertex(h);
		if(ph != -1 && !seen[ph]) dfs(ph, val.rotate90(4-cg.jump[h]));
	};
	for(const int h : m.corners()) if(!seen[h]) dfs(h, {0, 0});

	return quv;
}