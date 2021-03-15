// libMesh
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/elem.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/node.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/boundary_info.h"
// Formato ExodusII
#include "libmesh/exodusII_io.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <vector>


#define RAND ((double)(rand() >> 1)/((RAND_MAX >> 1) + 1))

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
// O info do ponto armazena o id de um nó (se houver) na casca (libmesh) original
typedef CGAL::Triangulation_vertex_base_with_info_3<libMesh::dof_id_type, K>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
//typedef Delaunay::Point                                             Point;

// Classe para construir uma malha libMesh usando CGAL 
// a partir de uma malha original com contornos definidos
namespace mesh3D {

  class Triangulation {

    public:
    Triangulation(libMesh::Mesh * in_mesh, libMesh::Mesh * out_mesh);

    void remesh();

    libMesh::Mesh * _in_mesh;
    libMesh::Mesh * _out_mesh;

  };

}
