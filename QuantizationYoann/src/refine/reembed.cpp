// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#include "refine/reembed.h"

#include "cgls.h"

#include <array>

using namespace std;

//======================
//== Reduction Matrix ==
//======================

struct ReductionMatrix {
	using Row = SparseVector<double>;
	inline static const int ONE = -1;

	ReductionMatrix() = default;
	ReductionMatrix(int n): _rows(n), refSize(n, 1) { for(const int i : Range(n)) _rows[i].emplace_back(i, 1.); }

	inline void addNewVar(int n) {
		_rows.resize(size()+n);
		refSize.resize(size(), 1);
		for(const int i : Range(n)) _rows[size()-1-i].emplace_back(size()-1-i, 1.);
	}
	inline void simplify(int v) { if(!isVar(v)) sub(_rows[v], refSize[v]); }
	void addEquality(Row &e);
	inline void addEquality(Row &&e) { addEquality(e); }
	inline bool isFixed(Row &e) { sub(e); return e.empty() || (e.size() == 1 && e[0].first == ONE); }
	inline bool isVar(int v) { return v == ONE || (_rows[v].size() == 1 && _rows[v][0].first == v); } 
	inline void reIndex() {
		_m = 0;
		vector<int> map_var(size());
		for(const int v : Range(size())) if(simplify(v), isVar(v)) map_var[v] = _m++;
		for(Row &r : _rows) for(auto &[i, a] : r) if(i != ONE) i = map_var[i];
	}

	inline std::size_t size() const { return _rows.size(); }
	inline int m() const { return _m; }
	inline const Row& getRow(int i) const { return _rows[i]; }
private:
	vector<Row> _rows;
	vector<uint32_t> refSize;
	int _m;

	void sub(Row &e, uint32_t ref=0);
};

//==========================
//== Re-embeding function ==
//==========================

void ES(const Mesh &m, const vector<vec2> &U, const ReductionMatrix &M, vector<double> &x, const double lambda=.05, const double delta=1.e-9) {
	vector<array<ReductionMatrix::Row, 4>> Jreduc(m.nfacets());
	vector<double> areas(m.nfacets());
	SymetricSystem SS(x.size(), 1e-6);
	vector<double> dx, x0, g;
	double eps2 = 1e9;

	// precompute
	for(const int f : m.facets()) {
		const mat2 dU(U[m.h(f,1)] - U[m.h(f,0)], U[m.h(f,2)] - U[m.h(f,0)]);
		const mat2 inv_dU = dU.inv();
		areas[f] = uvArea(m, U, f);
		for(const int j : {0,1}) for(const int i : {0,1})
			Jreduc[f][(j<<1)|i] =
				  inv_dU(i,0) * (M.getRow(2*m.h(f,1)+j) - M.getRow(2*m.h(f,0)+j))
				+ inv_dU(i,1) * (M.getRow(2*m.h(f,2)+j) - M.getRow(2*m.h(f,0)+j));
		for(const int i : Range(4)) for(const auto &[ind, a] : Jreduc[f][i]) if(ind != ReductionMatrix::ONE)
			for(const int j : Range(4))for(const auto &[ind2, a2] : Jreduc[f][j]) if(ind2 != ReductionMatrix::ONE)
				SS.addCoeff(ind, ind2);
	}
	const auto getJ = [&](int f) {
		mat2 J;
		for(const int i : {0,1}) for(const int j : {0,1}) {
			double val = 0.;
			for(const auto &[ind, a] : Jreduc[f][(i<<1)|j])
				if(ind == ReductionMatrix::ONE) val += a;
				else val += a * x[ind];
			J(i,j) = val;
		}
		return J;
	};
	const auto getF = [&]() {
		double F = 0.;
		const double C = lambda*lambda / ((1. + lambda) * sqrt(eps2 + 1.) - 1.);
		const double B = (2. * (1. + C) * sqrt(eps2 + 1.) - 2.) / C - 1.;
		for(const int f : m.facets()) {
			const mat2 J = getJ(f);
			const double detJ = J.det();
			const double detJ2 = detJ*detJ;
			const double eabs_detJ = sqrt(eps2 + detJ2);
			const double khi = .5 * (detJ > 0. ? detJ + eabs_detJ : eps2 / (eabs_detJ - detJ));
			F += areas[f] * (J.norm2() + C * (detJ2 + B)) / khi;
		}
		return F;
	};
	SS.setPattern();
	
	for(int step = 0; step < 200; ++step) {
		double minDet = numeric_limits<double>::max();
		for(const int f : m.facets()) minDet = min(minDet, getJ(f).det());
		const double detEPS = 1.e-14 + 5.e-2 * pow(min(0., minDet), 2);
		eps2 = min(max(.01*detEPS, eps2*.9), detEPS);
		if(!step && eps2 > 1.) eps2 = max(1., .1*detEPS);
		double F = 0.;
		SS.resetMb();
		const double C = lambda*lambda / ((1. + lambda) * sqrt(eps2 + 1.) - 1.);
		const double B = (2. * (1. + C) * sqrt(eps2 + 1.) - 2.) / C - 1.;
		for(const int f : m.facets()) {
			const mat2 J = getJ(f);
			const mat2 Jcom = J.com();
			const double detJ = J.det();
			const double detJ2 = detJ*detJ;
			const double eabs_detJ = sqrt(eps2 + detJ2);
			const double khi = .5 * (detJ > 0. ? detJ + eabs_detJ : eps2 / (eabs_detJ - detJ));
			const double F_num = J.norm2() + C * (detJ2 + B);
			const double area = areas[f];
			F += area * F_num / khi;
			const mat2 dJ = (area / khi) * (2.*J + (2.*C*detJ - F_num / eabs_detJ) * Jcom);
			const double daaJ = 2. * area / khi;
			const double daDJ = - 2. * area / (eabs_detJ * khi);
			const double dDDJ = 2. * area * (F_num / (eabs_detJ*eabs_detJ) + C * (1 - 2.*detJ/eabs_detJ)) / khi;
			for(const int i : Range(4)) {
				for(const auto &[ind, a] : Jreduc[f][i]) if(ind != ReductionMatrix::ONE) SS.addB(ind, a * dJ(i>>1,i&1));
				for(const int j : Range(4)) {
					double d2J = i==j ? daaJ : 0.;
					d2J += daDJ * (J(i>>1,i&1)*Jcom(j>>1,j&1) + Jcom(i>>1,i&1)*J(j>>1,j&1));
					d2J += dDDJ * (Jcom(i>>1,i&1)*Jcom(j>>1,j&1));
					for(const auto &[ind, a] : Jreduc[f][i]) if(ind != ReductionMatrix::ONE) {
						const double mul = d2J*a;
						for(const auto &[ind2, a2] : Jreduc[f][j]) if(ind2 != ReductionMatrix::ONE)
							SS.addCoeff(ind, ind2, mul*a2);
					}
				}
			}
		}
		g.assign(SS.getB().begin(), SS.getB().end());
		SS.precondion();
		solve(SS, dx);
		constexpr double alpha = .25, beta = .5;
		double t = 1., gdx = 0.;
		x0.assign(x.begin(), x.end());
		for(const int i : Range(x.size())) {
			x[i] -= dx[i];
			gdx += g[i]*dx[i];
		}
		while(getF() > F - alpha * t * gdx) {
			t *= beta;
			for(const int i : Range(x.size())) x[i] = x0[i] - t*dx[i];
		}
		cerr << step << ": " << minDet << ' ' << sqrt(eps2) << ' ' << F << ' ' << gdx/F << " " << t << endl;
		if(minDet > 0. && gdx/F < delta) break;
	}
}

void ES2(const Mesh &m, const vector<vec2> &U, const ReductionMatrix &M, vector<double> &x, const double lambda=.05, const double delta=1.e-9) {
	vector<array<ReductionMatrix::Row, 4>> Jreduc(m.nfacets());
	vector<double> areas(m.nfacets());
	double eps2 = 1e9;

	// precompute
	for(const int f : m.facets()) {
		const mat2 dU(U[m.h(f,1)] - U[m.h(f,0)], U[m.h(f,2)] - U[m.h(f,0)]);
		const mat2 inv_dU = dU.inv();
		areas[f] = uvArea(m, U, f);
		for(const int j : {0,1}) for(const int i : {0,1})
			Jreduc[f][(j<<1)|i] =
				  inv_dU(i,0) * (M.getRow(2*m.h(f,1)+j) - M.getRow(2*m.h(f,0)+j))
				+ inv_dU(i,1) * (M.getRow(2*m.h(f,2)+j) - M.getRow(2*m.h(f,0)+j));
	}
	const auto getJ = [&](int f) {
		mat2 J;
		for(const int i : {0,1}) for(const int j : {0,1}) {
			double val = 0.;
			for(const auto &[ind, a] : Jreduc[f][(i<<1)|j])
				if(ind == ReductionMatrix::ONE) val += a;
				else val += a * x[ind];
			J(i,j) = val;
		}
		return J;
	};
	const auto getF = [&]() {
		double F = 0.;
		const double C = lambda*lambda / ((1. + lambda) * sqrt(eps2 + 1.) - 1.);
		const double B = (2. * (1. + C) * sqrt(eps2 + 1.) - 2.) / C - 1.;
		for(const int f : m.facets()) {
			const mat2 J = getJ(f);
			const double detJ = J.det();
			const double detJ2 = detJ*detJ;
			const double eabs_detJ = sqrt(eps2 + detJ2);
			const double khi = .5 * (detJ > 0. ? detJ + eabs_detJ : eps2 / (eabs_detJ - detJ));
			F += areas[f] * (J.norm2() + C * (detJ2 + B)) / khi;
		}
		return F;
	};
	
	const int MH = 8;
	const int N = x.size();
	vector<double> rho(MH), alpha(MH), s(MH*N), y(MH*N);
	vector<double> dx(N), g(N), x0(N), g0(N), diag(N);
	for(int step = 0; step < 2000; ++step) {
		double minDet = numeric_limits<double>::max();
		for(const int f : m.facets()) minDet = min(minDet, getJ(f).det());
		if(step%16==0 && eps2*1e16 > 1.0001) {
			const double detEPS = 1.e-14 + 5.e-2 * pow(min(0., minDet), 2);
			eps2 = min(max(0.01*detEPS, eps2*.85), 1.e-14 + 5.e-2 * pow(min(0., minDet), 2));
			if(!step && eps2 > 1.) eps2 = max(1., .1*detEPS);
			diag.assign(N, 1.);
		}
		double F = 0.;
		const double C = lambda*lambda / ((1. + lambda) * sqrt(eps2 + 1.) - 1.);
		const double B = (2. * (1. + C) * sqrt(eps2 + 1.) - 2.) / C - 1.;
		g.assign(N, 0.);
		for(const int f : m.facets()) {
			const mat2 J = getJ(f);
			const double detJ = J.det();
			const double detJ2 = detJ*detJ;
			const double eabs_detJ = sqrt(eps2 + detJ2);
			const double khi = .5 * (detJ > 0. ? detJ + eabs_detJ : eps2 / (eabs_detJ - detJ));
			const double F_num = J.norm2() + C * (detJ2 + B);
			const double area = areas[f];
			F += area * F_num / khi;
			const mat2 dJ = (area / khi) * (2.*J + (2.*C*detJ - F_num / eabs_detJ) * J.com());
			for(const int i : Range(4)) for(const auto &[ind, a] : Jreduc[f][i]) if(ind != ReductionMatrix::ONE) g[ind] += a * dJ(i>>1,i&1);
		}
		// LBFGS
		dx.assign(g.begin(), g.end());
		if(step) {
			const int k0 = step%MH;
			const int i0 = k0*N;
			const int K = min(step, MH);
			double* const s0 = &s[i0];
			double* const y0 = &y[i0];
			double yy = 0., ss = 0., ys = 0.;
			for(const int i : Range(N)) {
				s0[i] = x[i] - x0[i];
				y0[i] = g[i] - g0[i];
				ss += s0[i]*s0[i]/diag[i];
				yy += y0[i]*y0[i]*diag[i];
				ys += s0[i]*y0[i];
			}
			rho[k0] = 1. / ys;
			for(int k = 0; k < K; ++k) {
				const int kk = (k0+MH-k)%MH;
				const int ik = kk*N;
				double sdx = 0.;
				const double* const sk = &s[ik];
				for(const int i : Range(N)) sdx += sk[i]*dx[i];
				alpha[kk] = rho[kk] * sdx;
				const double* const yk = &y[ik];
				for(const int i : Range(N)) dx[i] -= alpha[kk]*yk[i];
			}
			for(const int i : Range(N)) {
				diag[i] *= ys / (yy + y0[i]*y0[i]*diag[i] - yy*s0[i]*s0[i] / (ss*diag[i]));
				dx[i] *= diag[i];
			}
			for(int k = K-1; k >= 0; --k) {
				const int kk = (k0+MH-k)%MH;
				const int ik = kk*N;
				double beta = 0.;
				const double* const yk = &y[ik];
				for(const int i : Range(N)) beta += yk[i]*dx[i];
				beta = alpha[kk] - rho[kk]*beta;
				const double* const sk = &s[ik];
				for(const int i : Range(N)) dx[i] += beta*sk[i];
			}
		}
		x0.assign(x.begin(), x.end());
		g0.assign(g.begin(), g.end());
		//=======
		constexpr double alpha = .25, beta = .5;
		double t = 1., gdx = 0.;
		for(const int i : Range(x.size())) {
			x[i] -= dx[i];
			gdx += g[i]*dx[i];
		}
		while(t > 1e-10 && getF() > F - alpha * t * gdx) {
			t *= beta;
			for(const int i : Range(x.size())) x[i] = x0[i] - t*dx[i];
		}
		if((step%10)==0 || step==1999 || (minDet > 0. && eps2 < 1e-15 && gdx/F < delta)) cerr << step << ": " << minDet << ' ' << sqrt(eps2) << ' ' << F << ' ' << gdx/F << ' ' << t << endl;
		// cerr << step << ": " << minDet << ' ' << sqrt(eps2) << ' ' << F << ' ' << gdx/F << ' ' << t << endl;
		if(minDet > 0. && eps2 < 1e-15 && gdx/F < delta) break;
	}
}

vector<vec2> reembed(const Mesh &m, const vector<vec2> &uv, const CutGraph &cg, const vector<bool> &feature,
							const Mesh &qm, const vector<vec2> &quv, const CutGraph &qcg, const vector<bool> &qfeature,
							const SparseMatrix &D) {
	// Create system
	ReductionMatrix M(2 * m.ncorners());

	// Add feature equalities and zero jumps across non cutgraph edges
	for(const int h : m.corners()) {
		int opp = m.opp(h);
		if(feature[h] && opp < h) {
			const vec2 a = uv[m.next(h)] - uv[h];
			const int d = abs(a.y) < abs(a.x);
			M.addEquality(M.getRow(2*h+d) - M.getRow(2*m.next(h)+d));
		}
		if(opp == -1 || cg.isCut[h]) continue;
		opp = m.next(opp);
		for(const int d : {0,1}) M.addEquality(M.getRow(2*h+d) - M.getRow(2*opp+d));
	}

	// Add cut graph jumps
	for(const int v : Range(M.size())) M.simplify(v);
	for(const int h : m.corners()) {
		M.simplify(2*h);
		M.simplify(2*h+1);
		if(!cg.isCut[h]) continue;
		const int opp = m.opp(h);
		if(opp < h) continue;
		const int nh = m.next(h);
		const int nopp = m.next(opp);
		vec2T<ReductionMatrix::Row> oppVec(
			M.getRow(2*nopp) - M.getRow(2*opp),
			M.getRow(2*nopp+1) - M.getRow(2*opp+1)
		);
		oppVec.rotate90(cg.jump[h]);
		M.addEquality(M.getRow(2*nh  ) - M.getRow(2*h  ) + oppVec.x);
		M.addEquality(M.getRow(2*nh+1) - M.getRow(2*h+1) + oppVec.y);
	}

	// Construct DM
	vector<ReductionMatrix::Row> DM(D.size());
	for(const int d : {0,1}) for(const int h : qm.corners()) {
		DM[2*h+d].emplace_back(ReductionMatrix::ONE, -quv[h][d]);
		for(const auto &[j, a] : D.getRow(2*h+d)) DM[2*h+d] += a * M.getRow(j);
	}
	
	// Relax non-integer variables
	int nNewVar = 0;
	for(const int v : qm.verts()) {
		const vec2b h2bi = qcg.getHas2BeInt(quv, qfeature, v);
		for(const int d0 : {0, 1}) if(!h2bi[d0]) {
			int d = d0;
			for(const int h : qm.one_ring(v)) {
				DM[2*h + (d&1)].emplace_back(M.size() + nNewVar, (d&2) ? -1. : 1.);
				d = (d + 4 - qcg.jump[qm.prev(h)])&3;
			}
			++ nNewVar;
		}
	}
	M.addNewVar(nNewVar);

	// Add Quantization constraints
	for(const int d : {0,1}) for(const int h : qm.corners())
		M.addEquality(DM[2*qm.next(h)+d] - DM[2*h+d]);
	
	// Fix singularities
	for(const int v : qm.verts()) {
		const vec2b h2bi = qcg.getHas2BeInt(quv, qfeature, v);
		for(const int d : {0,1}) if(h2bi[d] && !M.isFixed(DM[2*qm.v2h(v)+d]))
			M.addEquality(DM[2*qm.v2h(v)+d]);
	}
	M.reIndex();

	// Rescale objective param
	vector<vec2> U = uv;
	double area0 = 0., areaQ = 0.;
	for(const int f : m.facets()) area0 += uvArea(m, uv, f);
	for(const int f : qm.facets()) areaQ += uvArea(qm, quv, f);
	const double scale_ratio = sqrt(areaQ / area0);
	for(vec2 &v : U) v *= scale_ratio;
	
	// Solve Least Square System
	LeastSquareSystem LS(M.m(), 1e-8);
	for(const int f : m.facets()) {
		const mat2 dU(U[m.h(f,1)] - U[m.h(f,0)], U[m.h(f,2)] - U[m.h(f,0)]);
		const double sqrt_Area = sqrt(dU(0,0)*dU(1,1) - dU(0,1)*dU(1,0));
		const double inv_sqrt_Area = 1. / sqrt_Area;
		const mat2 scaled_inv_dU(inv_sqrt_Area*vec2(dU(1,1), -dU(0,1)), inv_sqrt_Area*vec2(-dU(1,0), dU(0,0)));
		for(const int i : {0,1}) for(const int j : {0,1}) {
			const ReductionMatrix::Row Jji =
				  scaled_inv_dU(i,0) * (M.getRow(2*m.h(f,1)+j) - M.getRow(2*m.h(f,0)+j))
				+ scaled_inv_dU(i,1) * (M.getRow(2*m.h(f,2)+j) - M.getRow(2*m.h(f,0)+j));
			double rhs = i==j ? sqrt_Area : 0.;
			for(const auto &[ind, a] : Jji)
				if(ind == ReductionMatrix::ONE) rhs -= a;
				else LS.addCoeff(ind, a);
			LS.addRow(rhs);
		}
	}
	vector<double> x;
	LS.compact();
	solve(LS, x, &cout);
	ES2(m, U, M, x);

	// Translate to uv
	for(const int h : m.corners()) for(const int d : {0, 1}) {
		double u = 0.;
		for(const auto &[i, a] : M.getRow(2*h+d))
			if(i == ReductionMatrix::ONE) u += a;
			else u += a*x[i];
		U[h][d] = u;
	}

	return U;
}

//=====================================
//== Reduction Matrix Implementation ==
//=====================================

void ReductionMatrix::sub(Row &e, uint32_t ref) {
	Row new_e;
	for(const auto &[i, a] : e) {
		if(isVar(i)) new_e.emplace_back(i, a);
		else {
			sub(_rows[i], refSize[i]);
			refSize[i] -= ref;
			transform(_rows[i].begin(), _rows[i].end(), back_inserter(new_e), [a=a](const Row::value_type &c) {
				return Row::value_type{c.first, a*c.second};
			});
		}
	}
	sort(new_e.begin(), new_e.end());
	int s = 0, i = 0;
	while(i < (int) new_e.size()) {
		new_e[s] = new_e[i++];
		while(i < (int) new_e.size() && new_e[i].first == new_e[s].first) {
			new_e[s].second += new_e[i++].second;
			if(new_e[s].first != ONE) refSize[new_e[s].first] -= ref;
		}
		if(!Row::isNull(new_e[s].second)) ++ s;
		else if(new_e[s].first != ONE) refSize[new_e[s].first] -= ref;
	}
	new_e.resize(s);
	e = std::move(new_e);
}

struct ImpossibleEquality : std::exception {
	ImpossibleEquality(double v) { sprintf(msg, "You tried to set an expression to zero while its value is already set to %f", v); }
	const char* what() const noexcept override { return msg; }
private:
	char msg[100];
};

void ReductionMatrix::addEquality(Row &e) {
	sub(e);
	if(e.empty()) return;
	if(e.size() == 1 && e[0].first == ONE) throw ImpossibleEquality(e[0].second);
	int v = e[0].first == ONE;
	for(int w = v+1; w < (int) e.size(); ++w)
		if(refSize[e[w].first] < refSize[e[v].first])
			v = w;
	const int ind = e[v].first;
	const double mul = - 1. / e[v].second;
	for(int w = v+1; w < (int) e.size(); ++w) e[w-1] = e[w];
	e.pop_back();
	for(auto &[i, a] : e) {
		a *= mul;
		if(i != ONE) refSize[i] += refSize[ind];
	}
	_rows[ind] = move(e);
}