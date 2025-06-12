// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#include "mesh.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <unordered_set>

using namespace std;

void Mesh::compute_opp() {
	_opp.assign(ncorners(), -1);
	_v2h.resize(nverts());
	vector<tuple<int, int, int>> edges;
	edges.reserve(ncorners());
	for(const int h : corners()) {
		int a = from(h), b = to(h);
		_v2h[a] = h;
		if(a > b) swap(a, b);
		edges.emplace_back(a, b, h);
	}
	sort(edges.begin(), edges.end());
	for(int i = 0; i+1 < (int) edges.size(); ++i) {
		if(get<0>(edges[i]) != get<0>(edges[i+1]) || get<1>(edges[i]) != get<1>(edges[i+1])) continue;
		_opp[get<2>(edges[i])] = get<2>(edges[i+1]);
		_opp[get<2>(edges[i+1])] = get<2>(edges[i]);
		++ i;
	}
	for(const int v : verts()) {
		const int h0 = _v2h[v];
		while(true) {
			const int h = prev_around_vertex(_v2h[v]);
			if(h == -1 || h == h0) break;
			_v2h[v] = h;
		}
	}
}

template<> struct std::hash<pair<int, int>> {
	size_t operator()(const pair<int, int> &p) const {
		return (size_t(p.first) << (8*(sizeof(size_t) - sizeof(int)))) ^ p.second;
	}
};

tuple<Mesh, vector<vec2>, vector<bool>> readObj(const string filename) {
	vector<vec2> VT, uv;
	unordered_set<pair<int,int>> edges;
	Mesh m;
	ifstream in(filename);
	if(in.fail()) {
		cerr << "Failed to open " << filename << endl;
		exit(1);
	}
	string line;
	while(!in.eof()) {
		getline(in, line);
		istringstream iss(line);
		string type_token;
		iss >> type_token;
		if(type_token == "v") {
			m.points.emplace_back();
			for(int i=0; i < 3; ++i) iss >> m.points.back()[i];
		} else if(type_token == "vt") {
			VT.emplace_back();
			for(int i=0; i < 2; ++i) iss >> VT.back()[i];
		} else if(type_token == "f") {
			while(true) { // in wavefront obj all indices start at 1, not zero
				while(!iss.eof() && !isdigit(iss.peek())) iss.get(); // skip (esp. trailing) white space
				if(iss.eof()) break;
				int tmp;
				iss >> tmp;
				m.h2v.push_back(tmp - 1);
				uv.emplace_back();
				if(iss.peek() == '/') {
					iss.get();
					if(iss.peek() == '/') {
						iss.get();
						iss >> tmp;
					} else {
						iss >> tmp;
						uv.back() = VT[tmp-1];
						if(iss.peek() == '/') {
							iss.get();
							iss >> tmp;
						}
					}
				}
			}
		} else if(type_token == "l") {
			int i, j;
			iss >> i;
			for(; iss >> j; i = j) edges.emplace(min(i, j)-1, max(i, j)-1);
		}
	}
	in.close();
	vector<bool> feature(m.ncorners(), false);
	for(const int h : m.corners()) feature[h] = edges.count(pair<int, int>(min(m.from(h), m.to(h)), max(m.from(h), m.to(h))));
	m.compute_opp();
	for(const int h : m.corners()) if(m.opp(h) == -1) feature[h] = true;
	return make_tuple(move(m), move(uv), move(feature));
}

void writeObj(const std::string filename, const Mesh &m, const vector<vec2> &uv, const vector<bool> &feature) {
	ofstream out(filename);
	if(out.fail()) throw runtime_error("Failed to open " + filename);
	out << setprecision(numeric_limits<double>::max_digits10);
	out << "mtllib texture.mtl\n";
	for(const auto& p : m.points) out << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
	for(const auto& p : uv) out << "vt " << p[0] << ' ' << p[1] << '\n';
	for(const int f : m.facets()) {
		out << "f ";
		for(int i : Range(3)) out << m.v(f, i) + 1 << "/" << m.h(f, i) + 1 << " ";
		out << '\n';
	}
	if(!feature.empty()) for(const int h : m.corners()) if(feature[h] && m.opp(h) < h) out << "l " << m.from(h)+1 << ' ' << m.to(h)+1 << '\n';
	out.close();
}
