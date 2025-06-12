// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#include "decimation.h"
#include "quantization.h"
#include "refine/imprint.h"
#include "refine/reembed.h"
#include <fstream>

using namespace std;

void help(const char* prog) {
	cerr << "Usage: " << prog << " [OPTION] seamless_map.obj" << endl;
	cerr << "OPTION:" << endl;
	cerr << "\t-c: for coarse quantization" << endl;
	cerr << "\t-d `matrix file path`: output the matrix in the file given in parameter if `scale` is a numeric value. Otherwhise if `scale` is \"a\", an automatic scale value is computed." << endl;
	cerr << "\t-h: show this help." << endl;
	cerr << "\t-i: imprint the quantized map in the original mesh." << endl;
	cerr << "\t-r: re-embed the quantized map in the original mesh." << endl;
	cerr << "\t-o: `ouput path`: change the file in which the output will be written. Default is \"out.obj\"." << endl;
	cerr << "\t-s: `scale`: scale the input seamless map by `scale` if `scale` is a numeric value. Otherwhise if `scale` is \"a\", an automatic scale value is computed." << endl;
	cerr << "\t-sa: `scale` the auto scale." << endl;
}

enum OUTPUT {
	DECIMATED,
	IMPRINT,
	REEMBED
};

int main(int argc, const char* argv[]) {
	// Read arguments
	if(argc < 2) {
		cerr << "This executable requires at least one argument: the input seamless map file in obj format." << endl;
		help(argv[0]);
		return 1;
	}
	double scale = 1.;
	double scale_auto = 1.;
	OUTPUT output = DECIMATED;
	bool coarse_energy = false;
	string out_name = "out.obj";
	string matrix_name = "";
	for(int i = 1; i+1 < argc; ++i) {
		if(string(argv[i]) == "-s") {
			if(string(argv[i+1]) == "a") scale = -1, ++i;
			else scale = atof(argv[++i]);
		} else if(string(argv[i]) == "-h") {
			help(argv[0]);
			return 0;
		} else if(string(argv[i]) == "-i") output = IMPRINT;
		else if(string(argv[i]) == "-sa") scale_auto = atof(argv[++i]);
		else if(string(argv[i]) == "-r") output = REEMBED;
		else if(string(argv[i]) == "-o") out_name = argv[++i];
		else if(string(argv[i]) == "-c") coarse_energy = true;
		else if(string(argv[i]) == "-d") matrix_name = argv[++i];
	}

	// Read input
	auto [m, seamless, feature] = readObj(argv[argc-1]);
	cerr << "INPUT: \t" << m.nverts() << " verts \t" << m.nfacets() << " facets" << endl; 
	if(scale == -1) {
		double surface_area = 0.;
		for(const int f : m.facets()) surface_area += uvArea(m, seamless, f);
		scale = .3 * sqrt(m.nfacets() / surface_area);
		scale *= scale_auto;
		cerr << "AUTOMATIC RESCALE: " << scale << endl;
	}
	for(vec2 &v : seamless) v *= scale;

	// Check validity of input
	for(const int f : m.facets())
		if((seamless[m.h(f, 1)].x - seamless[m.h(f, 0)].x)*(seamless[m.h(f, 2)].y - seamless[m.h(f, 0)].y)
				- (seamless[m.h(f, 1)].y - seamless[m.h(f, 0)].y)*(seamless[m.h(f, 2)].x - seamless[m.h(f, 0)].x) <= 0.) {
			cerr << "[ERROR] The input is not valid! Facet " << f << " has negative det." << endl;
			return 1;
		}

	// Get cut graph
	CutGraph cg = compute_cut_graph(m, seamless);

	// Save initial mesh
	const Mesh fine = m;
	const vector<vec2> f_uv = seamless;
	const vector<bool> f_feature = feature;
	const CutGraph f_cg = cg;

	// Decimation
	SparseMatrix D;
	if(!matrix_name.empty() || output == REEMBED) D.setIdentity(2*m.ncorners());
	do { makeDelaunay(m, cg, seamless, feature, D); } while(collapse(m, cg, seamless, feature, D));
	cerr << "DECIMATION: \t" << m.nverts() << " verts \t" << m.nfacets() << " facets" << endl; 
	if(!matrix_name.empty()) {
		ofstream dfile(matrix_name);
		dfile << D.size() << '\n';
		for(const SparseMatrix::Row &r : D.rows()) dfile << r << '\n';
		dfile.close();
	}

	// Quantization
	cerr << "QUANTIZATION" << endl;
	const vector<vec2i> quv0 = quantize(m, cg, seamless, feature, coarse_energy);
	const vector<vec2> quv(quv0.begin(), quv0.end());

	// Output
	switch(output) {
	case DECIMATED:
		writeObj(out_name, m, quv, feature);
		break;
	case IMPRINT: {
		cerr << "IMPRINTING..." << endl;
		const auto [fineQ, uvQ] = transfert_coarse_to_fine(fine, f_uv, m, seamless, quv);
		writeObj(out_name, fineQ, uvQ);
		break;
	} case REEMBED: {
		cerr << "RE-EMBEDING..." << endl;
		const vector<vec2> U = reembed(fine, f_uv, f_cg, f_feature, m, quv, cg, feature, D);
		writeObj(out_name, fine, U, f_feature);
		break;
	}
	}

	return 0;
}
