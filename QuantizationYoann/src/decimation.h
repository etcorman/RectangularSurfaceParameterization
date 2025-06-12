// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#pragma once

#include "mesh.h"
#include "cut_graph.h"
#include "matrix.h"

void makeDelaunay(Mesh &m, CutGraph &cg, std::vector<vec2> &uv, std::vector<bool> &feature, SparseMatrix &D);

bool collapse(Mesh &m, CutGraph &cg, std::vector<vec2> &uv, std::vector<bool> &feature, SparseMatrix &D);