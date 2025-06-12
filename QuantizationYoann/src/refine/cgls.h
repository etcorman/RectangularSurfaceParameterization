// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#pragma once

#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

class LinearSystem {
public:
	virtual int N() const = 0;
	virtual int NNZ() const = 0;
	virtual int MSize() const = 0;
	virtual std::tuple<const double*, const int*, const int*, const double*> Mb() const = 0;
	virtual double threshold() const = 0;

	virtual void preProcess(std::vector<double> &x) const = 0;
	virtual void postProcess(std::vector<double> &x) const = 0;
};

void solve(const LinearSystem &LS, std::vector<double> &x, std::ostream* out_stream = nullptr);

class BasicLinearSystem : public LinearSystem {
protected:
	int n, nnz, Msize;
	double eps, eps2;
	std::vector<double> P, b;
	double* Mcoeff = nullptr;
	int* Min;
	int* off;

	inline void clean() {
		if(!Mcoeff) return;
		delete[] Mcoeff;
		delete[] Min;
		delete[] off;
		Mcoeff = nullptr;
	}

public:
	BasicLinearSystem(int n, double eps): n(n), eps(eps) {}
	~BasicLinearSystem() { clean(); }

	inline int N() const { return n; }
	inline int NNZ() const { return nnz; }
	inline int MSize() const { return Msize; }
	inline std::tuple<const double*, const int*, const int*, const double*> Mb() const { return {Mcoeff,Min,off,b.data()}; };
	inline double threshold() const { return eps2; };

	inline void preProcess(std::vector<double> &x) const { for(int i = 0; i < n; ++i) x[i] /= P[i]; }
	inline void postProcess(std::vector<double> &x) const { for(int i = 0; i < n; ++i) x[i] *= P[i]; }
};

class SymetricSystem : public BasicLinearSystem {
private:
	std::vector<std::vector<int>> A;
public:
	SymetricSystem(int n, double eps=1e-8): BasicLinearSystem(n, eps), A(n) {}

	inline void addCoeff(int i, int j) { if(std::find(A[i].begin(), A[i].end(), j) == A[i].end()) A[i].push_back(j); }
	void setPattern();
	inline void addB(int i, double v) { b[i] += v; }
	inline void addCoeff(int i, int j, double val) { if(j < i) Mcoeff[std::lower_bound(&Min[off[i]], &Min[off[i+1]], j)-Min] += val; else if(j == i) P[i] += val; }
	inline const std::vector<double>& getB() { return b; }
	void precondion();
	inline void resetMb() { std::fill_n(Mcoeff, Msize, 0.); P.assign(n, 0.); b.assign(n, 0.); }
};

class LeastSquareSystem : public BasicLinearSystem {
private:
	std::vector<std::pair<int, double>> row;
	std::vector<std::vector<std::pair<int, double>>> A;
public:
	LeastSquareSystem(int n, double eps=1e-8): BasicLinearSystem(n, eps), A(n) { P.assign(n, 0.); b.assign(n, 0.); }
	~LeastSquareSystem() { clean(); }

	void compact();
	inline void addCoeff(int i, double val) { row.emplace_back(i, val); }
	void addRow(double rhs=0.);
};