/**
 * Import the entire library into the tinynurbs namespace.
 * 
 * Use of this source code is governed by a BSD-style license that can be found in
 * the LICENSE file.
 */

#include "nurbs/basis.h"
#include "nurbs/check.h"
#include "nurbs/curve.h"
#include "nurbs/evaluate.h"
#include "nurbs/modify.h"
#include "nurbs/model.h"
#include "nurbs/refine.h"
#include "nurbs/surface.h"
#include "io/obj.h"
#include "io/ionurbs.h"
#include "util/array2.h"
#include "util/coord.h"
#include "iga/stiffness.h"
#include "iga/assembly.h"
