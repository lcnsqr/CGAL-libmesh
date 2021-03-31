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

#include <sys/random.h>
#define RAND ((double)(rand() >> 1)/((RAND_MAX >> 1) + 1))

// Classe para construir uma malha libMesh usando CGAL 
// a partir de uma malha original com contornos definidos
namespace mesh3D {

  // Origem do ponto inserido na nova malha
  enum NODE_SOURCE {
    ORIGINAL,
    ADDED
  };

  // Classe do campo de informação usado em cada nó da malha
  class VertexInfo {
    public:
    enum NODE_SOURCE source;
    unsigned long int id;
  };

  class Triangulation {

    public:
    Triangulation(libMesh::Mesh * in_mesh, libMesh::Mesh * out_mesh);

    void remesh();

    libMesh::Mesh * _in_mesh;
    libMesh::Mesh * _out_mesh;

  };

}

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
// O info do ponto armazena o id de um nó (se houver) na casca (libmesh) original
typedef CGAL::Triangulation_vertex_base_with_info_3<mesh3D::VertexInfo, K>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
//typedef Delaunay::Point                                             Point;
