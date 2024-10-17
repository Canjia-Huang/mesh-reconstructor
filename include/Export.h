#ifndef MESH_RECONSTRUCTOR_EXPORT_H_
#define MESH_RECONSTRUCTOR_EXPORT_H_

# ifdef MESH_RECONSTRUCTOR_EXPORTS_
/* We are building this library */
#      define MESH_RECONSTRUCTOR_EXPORT __declspec(dllexport)
# else
/* We are using this library */
#      define MESH_RECONSTRUCTOR_EXPORT __declspec(dllimport)
# endif

#ifdef RPD_VERBOSE
#define RPD_DEBUG_ONLY(x) x
#else
#define RPD_DEBUG_ONLY(x)
#endif

#define MR_EPS 1e-12

#endif
