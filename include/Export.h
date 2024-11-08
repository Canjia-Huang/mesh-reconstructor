#ifndef MESH_RECONSTRUCTOR_EXPORT_H_
#define MESH_RECONSTRUCTOR_EXPORT_H_

#include <iostream>

#ifdef DEBUG_VERBOSE
#define DEBUG_ONLY_COUT(x) std::cout << __FUNCTION__ << " " << x << std::endl
#else
#define DEBUG_ONLY_COUT(x)
#endif

#define MR_EPS 1e-12

#endif
