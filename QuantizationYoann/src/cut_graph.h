// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#pragma once

#include "mesh.h"

struct CutGraph {
	std::vector<bool> isCut;
	std::vector<int> jump;
	std::vector<int> valence;
	CutGraph(const Mesh &m): isCut(m.ncorners(), true), jump(m.ncorners()), valence(m.nverts()), m(m) {}
	inline bool isSingularity(int v) const { return valence[v] != (m.is_boundary_vert(v) ? 2 : 4); }
	vec2b getHas2BeInt(const std::vector<vec2> &uv, const std::vector<bool> &feature, int v) const;
	inline bool has2BeInt(const std::vector<vec2> &uv, const std::vector<bool> &feature, int v) const {
		const vec2b h2bi = getHas2BeInt(uv, feature, v); return h2bi.x && h2bi.y; 
	}
private:
	const Mesh &m;
};

CutGraph compute_cut_graph(const Mesh &m, std::vector<vec2> &uv);