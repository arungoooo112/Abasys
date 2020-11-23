/**
 * Import the entire library into the tinynurbs namespace.
 * 
 * Use of this source code is governed by a BSD-style license that can be found in
 * the LICENSE file.
 */

#include "src/nurbs/basis.h"
#include "src/nurbs/check.h"
#include "src/nurbs/curve.h"
#include "src/nurbs/evaluate.h"
#include "src/nurbs/modify.h"
#include "src/nurbs/model.h"
#include "src/nurbs/refine.h"
#include "src/nurbs/surface.h"
#include "src/io/obj.h"
#include "src/io/ionurbs.h"
#include "src/util/array2.h"
#include "src/util/coord.h"
#include "src/iga/stiffness.h"
#include "src/iga/assembly.h"
