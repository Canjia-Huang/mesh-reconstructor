#ifndef MESH_RECONSTRUCTOR_EXPORT_H_
#define MESH_RECONSTRUCTOR_EXPORT_H_

#include <iostream>

#ifdef DEBUG_VERBOSE
#define DEBUG_ONLY_COUT(x) std::cout << __FUNCTION__ << " " << x << std::endl
#else
#define DEBUG_ONLY_COUT(x)
#endif

#ifdef DEBUG_OUTPUT
#define DEBUG_OUTPUT_POISSON_MODEL  "Poisson_Model.obj"
#endif

#define MR_EPS 1e-12

#endif
