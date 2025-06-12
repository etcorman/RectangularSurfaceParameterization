// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#include "cgls.h"

#include <chrono>
#include <cmath>
#include <condition_variable>
#include <iomanip>
#include <mutex>
#include <thread>

#include "vectorize.h"

using namespace std;

void SymetricSystem::setPattern() {
	nnz = Msize = 0;
	for(vector<int> &a : A) {
		sort(a.begin(), a.end());
		nnz += a.size();
		Msize += (4 - (a.size()&3))&3;
	}
	clean();
	Mcoeff = new(v4d::ALI) double[nnz+Msize];
	Min = new(v4i::ALI) int[nnz+Msize];
	off = new int[n+1];
	off[0] = 0;
	Msize = 0;
	for(int i = 0; i < n; ++i) {
		for(int j : A[i]) Min[Msize++] = j;
		while(Msize&3) {
			Min[Msize] = Min[Msize-1];
			++Msize;
		}
		off[i+1] = Msize;
	}
}

void SymetricSystem::precondion() {
	eps2 = 0.;
	for(int i = 0; i < n; ++i) {
		if(P[i] <= 0.) P[i] = 1.;
		else P[i] = 1. / sqrt(P[i]);
		b[i] *= P[i];
		eps2 += b[i]*b[i];
		for(int j = off[i]; j < off[i+1]; ++j) {
			const int k = Min[j];
			if(k >= i) {
				if(k == i) Mcoeff[j] = 1.;
				break;
			}
			Mcoeff[j] *= P[i]*P[k];
			Mcoeff[lower_bound(&Min[off[k]], &Min[off[k+1]], i)-Min] = Mcoeff[j];
		}
	}
	eps2 *= eps*eps;
}

void LeastSquareSystem::addRow(double rhs) {
	for(const auto &[i, ai] : row) {
		for(const auto &[j, aj] : row) if(i <= j) {
			const double a = ai*aj;
			if(i == j) {
				P[i] += a;
				continue;
			}
			const auto it = find_if(A[i].begin(), A[i].end(), [j=j](const auto &c){ return c.first == j; });
			if(it == A[i].end()) A[i].emplace_back(j, a);
			else it->second += a;
		}
		b[i] += rhs * ai;
	}
	row.clear();
}

void LeastSquareSystem::compact() {
	eps2 = 0.;
	for(int i = 0; i < N(); ++i) {
		if(P[i] <= 0.) { P[i] = 1.; continue; }
		A[i].emplace_back(i, 1.);
		P[i] = 1. / sqrt(P[i]);
		b[i] *= P[i];
		eps2 += b[i]*b[i];
	}
	eps2 *= eps*eps;
	nnz = 0;
	Msize = 0;
	for(int i = 0; i < N(); ++i) {
		const double Pi = P[i];
		sort(A[i].begin(), A[i].end());
		for(int k = A[i].size(); k--;) {
			auto &[j, a] = A[i][k];
			if(j <= i) break;
			a *= Pi*P[j];
			A[j].emplace_back(i, a);
		}
		nnz += A[i].size();
		Msize += (4 - (A[i].size()&3))&3;
	}
	clean();
	Mcoeff = new(v4d::ALI) double[nnz+Msize];
	Min = new(v4i::ALI) int[nnz+Msize];
	off = new int[n+1];
	off[0] = 0;
	Msize = 0;
	for(int i = 0; i < N(); ++i) {
		for(const auto &[j, a] : A[i]) {
			Mcoeff[Msize] = a;
			Min[Msize] = j;
			++Msize;
		}
		for(; Msize&3; ++Msize) {
			Mcoeff[Msize] = 0.;
			Min[Msize] = Min[Msize-1];
		}
		off[i+1] = Msize;
	}
}

void solve(const LinearSystem &LS, vector<double> &x0, std::ostream* out_stream) {
	const int n = LS.N();
	const int n2 = n&3 ? n+4-(n&3) : n;
	chrono::high_resolution_clock::time_point t0;
	ios::ios_base::fmtflags out_flags;
	bool stdio_sync = true;
	if(out_stream) {
		stdio_sync = ios::sync_with_stdio(false);
		out_flags = out_stream->flags();
		*out_stream << scientific;
		*out_stream << "Solving system with " << n << " variables and " << LS.NNZ() << " non zero coeffs\n";
		t0 = chrono::high_resolution_clock::now();
	}

	if((int) x0.size() < n) { x0.reserve(n2); x0.resize(n, 0.); }
	LS.preProcess(x0);

	double *x = ((uintptr_t)x0.data()%(size_t)v4d::ALI) || ((int) x0.capacity() < n2) ? new(v4d::ALI) double[n2] : x0.data();
	double *r = new(v4d::ALI) double[n2];
	double *d = new(v4d::ALI) double[n2];
	double *h = new(v4d::ALI) double[n2];
	if(x != x0.data()) copy_n(x0.begin(), n, x);
	for(int i = n; i < n2; ++i) r[i] = h[i] = 0.;
    double curr_err = 0., alpha, beta;
	int k = 0;

	const int NT = thread::hardware_concurrency();
	int num_working = NT;
	int step = 0;
	mutex redution_mutex, working_mutex;
	condition_variable start_cv;

	const auto synchrof = [&](int &local_step, const auto f) {
		bool last = false;
		working_mutex.lock();
		if(--num_working == 0) last = true;
		working_mutex.unlock();
		if(last) {
			f();
			working_mutex.lock();
			num_working = NT;
			++ step;
			working_mutex.unlock();
			start_cv.notify_all();
		} 
		++ local_step;
		unique_lock lock(working_mutex);
		start_cv.wait(lock, [&](){ return step == local_step; });
	};
	const auto synchro = [&](int &local_step) { synchrof(local_step, [](){}); };

	const auto work = [&](const int t) {
		const auto [Mcoeff, Min, off, b] = LS.Mb();
		const int R0 = t == 0 ? 0 : lower_bound(off, off+n+1, (t*(LS.MSize()+3*n))/NT, [off=off](const int &a, const int b) { return a + 3*(&a-off) < b; }) - off;
		const int R1 = t == NT-1 ? n : lower_bound(off, off+n+1, ((t+1)*(LS.MSize()+3*n))/NT, [off=off](const int &a, const int b) { return a + 3*(&a-off) < b; }) - off;
		const int B0 = ((t*n2)/NT) & (~(0b11));
		const int B1 = (((t+1)*n2)/NT) & (~(0b11));
		double local_sum = 0.;
		int local_step = 0;
		const auto initReduct = [&]() { local_sum = 0.; };
		const auto reduct = [&](double &x) { redution_mutex.lock(); x += local_sum; redution_mutex.unlock(); };

		initReduct();
		for(int i = R0, j = off[i]; i < R1; ++i) {
			v4d tmp(-b[i], 0., 0., 0.);
			for(; j < off[i+1]; j += 4)
				tmp = fma(v4d(&Mcoeff[j], true), v4d(x, v4i(&Min[j], true)), tmp);
			const double stmp = tmp.sum();
			r[i] = d[i] = stmp;
			local_sum += stmp*stmp;
		}
		reduct(curr_err);
		synchro(local_step);
		if(curr_err < LS.threshold()) goto xTimesP;

		while(k < n) {
			synchrof(local_step, [&]() {
				if(out_stream && (k%100) == 0) *out_stream << setw(4) << k << " : " << curr_err << " -- " << LS.threshold() << '\n';
				alpha = 0.;
				++ k;
			});

			initReduct();
			for(int i = R0, j = off[i]; i < R1; ++i) {
				v4d tmp(0.);
				for(; j < off[i+1]; j += 4)
					tmp = fma(v4d(&Mcoeff[j], true), v4d(d, v4i(&Min[j], true)), tmp);
				const double stmp = tmp.sum();
				h[i] = stmp;
				local_sum += d[i]*stmp;
			}
			reduct(alpha);
			synchrof(local_step, [&]() {
				alpha = curr_err / alpha;
				beta = 1. / curr_err;
				curr_err = 0.;
			});

			const v4d A4(alpha);
			v4d E4(0.);
			for(int i = B0; i < B1; i += 4) {
				const v4d tmp = fnma(A4, v4d(&h[i], true), v4d(&r[i], true));
				fnma(A4, v4d(&d[i], true), v4d(&x[i], true)).store(&x[i], true);
				tmp.store(&r[i], true);
				E4 = fma(tmp, tmp, E4);
			}
			redution_mutex.lock(); curr_err += E4.sum(); redution_mutex.unlock();			
			synchro(local_step);
			if(curr_err < LS.threshold()) break;
			const v4d B4(beta * curr_err);
			for(int i = B0; i < B1; i += 4)
				fma(B4, v4d(&d[i], true), v4d(&r[i], true)).store(&d[i], true);
		}

		xTimesP:
		for(int i = B0; i < min(B1, n); ++i) x0[i] = x[i];
	};

	vector<thread> threads; threads.reserve(NT-1);
	for(int i = 1; i < NT; ++i) threads.emplace_back(work, i);
	work(0);
	for(int i = 1; i < NT; ++i) threads[i-1].join();
	LS.postProcess(x0);

	if(x != x0.data()) delete[] x;
	delete[] r;
	delete[] d;
	delete[] h;

	if(out_stream) {
		*out_stream << "Converged in " << k << " iterations and "
			<< chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now()-t0).count()
			<< "ms\n||Ax - b|| = " << sqrt(curr_err) << '\n';
		out_stream->setf(out_flags);
		ios::sync_with_stdio(stdio_sync);
	}
}