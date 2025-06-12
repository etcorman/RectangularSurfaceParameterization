// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#include "mesh.h"

std::pair<Mesh, std::vector<vec2>> transfert_coarse_to_fine(
	const Mesh &fine,
	const std::vector<vec2> &f_uv,
	const Mesh &coarse, 
	const std::vector<vec2> &c_uv,
	const std::vector<vec2> &c_uv2
);
