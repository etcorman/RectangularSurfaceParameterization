// Copyright (C) 2022, Coudert--Osmont Yoann
// SPDX-License-Identifier: AGPL-3.0-or-later
// See <https://www.gnu.org/licenses/>

#pragma once

#include "mesh.h"
#include "cut_graph.h"

std::vector<vec2i> quantize(const Mesh &m, const CutGraph &cg, const std::vector<vec2> &uv, const std::vector<bool> &feature, bool coarsest=false);