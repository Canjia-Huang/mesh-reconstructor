#ifndef MESH_RECONSTRUCTOR_EXPORT_H
#define MESH_RECONSTRUCTOR_EXPORT_H

#ifndef MESH_RECONSTRUCTOR_EXPORT
#  ifdef MESH_RECONSTRUCTOR_EXPORTS
/* We are building this library */
#      define MESH_RECONSTRUCTOR_EXPORT __declspec(dllexport)
#  else
/* We are using this library */
#      define MESH_RECONSTRUCTOR_EXPORT __declspec(dllimport)
#  endif
#endif

#endif
