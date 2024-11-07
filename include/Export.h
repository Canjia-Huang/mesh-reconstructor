#ifndef MESH_RECONSTRUCTOR_EXPORT_H_
#define MESH_RECONSTRUCTOR_EXPORT_H_

#include <iostream>

#ifdef RPD_VERBOSE
#define RPD_DEBUG_ONLY(x) std::cout << x << std::endl;
#else
#define RPD_DEBUG_ONLY(x)
#endif

#define MR_EPS 1e-12

#endif
