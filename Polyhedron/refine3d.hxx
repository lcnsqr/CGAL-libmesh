#include <string.h>
#include <math.h>

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <array>
#include <vector>

// libmesh include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/boundary_info.h"
// Formato ExodusII
#include "libmesh/exodusII_io.h"

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Mesh_polyhedron_3<Kernel>::type                Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel> Mesh_domain;
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
