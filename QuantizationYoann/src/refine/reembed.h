// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#pragma once

#include "decimation.h"

std::vector<vec2> reembed(
	const Mesh &m,
	const std::vector<vec2> &uv,
	const CutGraph &cg,
	const std::vector<bool> &feature,
	const Mesh &qm,
	const std::vector<vec2> &quv,
	const CutGraph &qcg,
	const std::vector<bool> &qfeature,
	const SparseMatrix &D
);