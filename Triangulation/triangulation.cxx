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

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
//typedef Delaunay::Point                                             Point;

// libMesh namespace
using namespace libMesh;

int main(int argc, char **argv) {

  // Iniciar libMesh e bibliotecas subjacentes (PETSc)
  LibMeshInit init (argc, argv);

  // Instanciar uma nova malha
  Mesh mesh(init.comm());

  // Dimensão da malha
  int dim = 3;

  // Gerar uma malha uniforme num cubo de lado 2
  Real halfside = 1.;
  MeshTools::Generation::build_cube (mesh, 1, 1, 1, -halfside, halfside, -halfside, halfside, -halfside, halfside, TET4);

  // Exibir informações da malha
  mesh.print_info();

  // Todos os nós que estão no contorno
  std::unordered_set<dof_id_type> boundary_nodes;
  boundary_nodes = libMesh::MeshTools::find_boundary_nodes(mesh);

  /*
  std::vector< std::pair<Delaunay::Point,unsigned> > points;
  points.push_back( std::make_pair(Delaunay::Point(0,0,0),0) );
  points.push_back( std::make_pair(Delaunay::Point(1,0,0),1) );
  points.push_back( std::make_pair(Delaunay::Point(0,1,0),2) );
  points.push_back( std::make_pair(Delaunay::Point(0,0,1),3) );
  points.push_back( std::make_pair(Delaunay::Point(2,2,2),4) );
  points.push_back( std::make_pair(Delaunay::Point(-1,0,1),5) );
  points.push_back( std::make_pair(Delaunay::Point(-1,2,1),6) );

  // [CGAL] 3D Delaunay triangulation of points set
  Delaunay T( points.begin(), points.end() );
  */

  // [CGAL] Instanciar uma triangulação Delaunay 3D vazia
  Delaunay T;


  // Passar por todos os nós ativos da malha original
  MeshBase::const_node_iterator node = mesh.active_nodes_begin(), end_node = mesh.active_nodes_end();
  for ( ; node != end_node ; node++ )
	{
    // Ponteiro para o nó
    const Node * nptr = *node;

    // Coordenadas do nó
    libMesh::Point p = *nptr;

    // Inserir ponto na triangulação
    Delaunay::Vertex_handle vh = T.insert(Delaunay::Point(p(0),p(1),p(2)));

    /*
    // Verificar se nó está na face externa do elemento
    if (boundary_nodes.find(nptr->id()) != boundary_nodes.end()){
      // Nó de contorno
      libMesh::out << "Nó de contorno (" << nptr->id() << "): " << p(0) << ", " << p(1) << ", " << p(2) << std::endl;
      vh->info() = 1;
		}
    else {
      // Nó interno
      libMesh::out << "Nó interno (" << nptr->id() << "): " << p(0) << ", " << p(1) << ", " << p(2) << std::endl;
      vh->info() = 0;
    }
    */

    // Marcar o nó inserido com seu ID na malha original
    vh->info() = nptr->id();

	}

	// Limpar malha original
	mesh.clear();
	mesh.set_mesh_dimension(dim);
	mesh.set_spatial_dimension(dim);

  // Número de células na malha nova
  std::cout << "T.number_of_finite_cells() = " << T.number_of_finite_cells() << std::endl;

	unsigned int elem_id = 0;

  // Percorrer todas as células da malha nova
  for (Delaunay::Finite_cells_iterator cit = T.finite_cells_begin(); cit != T.finite_cells_end(); cit++)
  {
    std::cout << "Cell" << std::endl;
    std::cout << cit->vertex(0)->info() << " -> " << cit->vertex(0)->point() << std::endl;
    std::cout << cit->vertex(1)->info() << " -> " << cit->vertex(1)->point() << std::endl;
    std::cout << cit->vertex(2)->info() << " -> " << cit->vertex(2)->point() << std::endl;
    std::cout << cit->vertex(3)->info() << " -> " << cit->vertex(3)->point() << std::endl;

		/*
		Elementos que compartilham nós fazem com que 
		os nós compartilhados sejam duplicados na malha
		*/

    /*
    Usar um bitmap para marcar os nós que já foram inseridos.
    Usar um vetor com os Node* que foram inseridos.
    Colocar no elemento os nós que estão no vetor.
    */


    /*
		// Novo elemento na malha nova
		libMesh::Elem * elem = new libMesh::Tet4;
		elem->set_id(elem_id++);

		// Percorre os 4 pontos do tetraedro e cria os nós correspondentes
		std::vector< std::pair<libMesh::Node *, unsigned> > nodes;
		for (int n = 0; n < 4; n++)
		{

			nodes.push_back( std::make_pair(

				// Coordenadas
				mesh.add_point(libMesh::Point(
					cit->vertex(n)->point().x(), 
					cit->vertex(n)->point().y(), 
					cit->vertex(n)->point().z()
				)),

				// Informação de contorno
				cit->vertex(n)->info())

			);

			// Definir n-ésimo nó do elemento
			elem->set_node(n) = nodes[n].first;
		}

		// Incluir elemento na malha nova
		elem = mesh.add_elem(elem);
    */




		/*

		libMesh::Node * n0 =vertices_map[it->vertex(0)];
		libMesh::Node * n1 =vertices_map[it->vertex(1)];
		libMesh::Node * n2 =vertices_map[it->vertex(2)];
		libMesh::Node * n4 =vertices_map[it->vertex(4)];


		delete elem;
		*/

  }

	// Exibir informações da malha reconstruída
  mesh.print_info();

  return 0;
}
