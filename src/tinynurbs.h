/**
 * Import the entire library into the tinynurbs namespace.
 * 
 * Use of this source code is governed by a BSD-style license that can be found in
 * the LICENSE file.
 */

#include "core/basis.h"
#include "core/check.h"
#include "core/curve.h"
#include "core/evaluate.h"
#include "core/modify.h"
#include "core/model.h"
#include "core/refine.h"
#include "core/surface.h"
#include "io/obj.h"
#include "io/ionurbs.h"
#include "util/array2.h"
#include "util/coord.h"
#include "iga/stiffness.h"
#include "iga/assembly.h"
