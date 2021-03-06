#include "Boundaries.hxx"
#include "Triangulation.hxx"

// libMesh namespace
using namespace libMesh;

mesh3D::Triangulation::Triangulation(libMesh::Mesh * in_mesh, libMesh::Mesh * out_mesh){
  this->_in_mesh = in_mesh;
  this->_out_mesh = out_mesh;
}

void mesh3D::Triangulation::remesh(){
  // Contorno do casco, assumindo que existem nós associados a contornos
  BoundaryInfo & in_boundary_info = _in_mesh->get_boundary_info();

  // Exibir informações da malha
#ifdef DEBUG
  _in_mesh->print_info();
  in_boundary_info.print_summary();
#endif

  // Vetor de tuplas dof_id/boundary_id referente a nós e contornos
  typedef std::vector<std::tuple<dof_id_type,boundary_id_type>> nodes_boundaries_type;
  nodes_boundaries_type nodes_boundaries;
  nodes_boundaries = in_boundary_info.build_node_list();

  // Mapa com cada nó original e seus contornos correspondentes
  std::map<dof_id_type, std::vector<boundary_id_type>> input_node_boundaries;

  // Preencher o mapa 
  for (nodes_boundaries_type::iterator nb_it = nodes_boundaries.begin(); 
    nb_it != nodes_boundaries.end();
    ++nb_it)
    {
       input_node_boundaries[std::get<0>(*nb_it)].push_back(std::get<1>(*nb_it));
    }

  // [CGAL] Instanciar uma triangulação Delaunay 3D vazia
  Delaunay T;

  // [CGAL] Nós inseridos na malha CGAL
  std::map<dof_id_type, bool> cgal_added_nodes;

  // Inserir nós originais na malha CGAL
  for (auto & elem : _in_mesh->element_ptr_range())
    for (auto s : elem->side_index_range())
      if (elem->neighbor_ptr(s) == nullptr)
      {
        // Em cada face de elemento de contorno, gerar um mapa 
        // de contornos e seus nós correspondentes
        std::map<boundary_id_type, std::vector<dof_id_type>> input_boundary_nodes;

        // Face do elemento
        std::unique_ptr<Elem> side = elem->side_ptr(s);

        // Percorrer cada nó na face
        for (auto n : side->node_index_range())
        {
          // Identificar boundary IDs associados ao nó
          for (std::vector<boundary_id_type>::iterator bid_it = input_node_boundaries[side->node_id(n)].begin();
               bid_it != input_node_boundaries[side->node_id(n)].end();
               ++bid_it)
          {
            // Incluir nó no vetor correspondente a cada um de seus boundary ID
            input_boundary_nodes[*bid_it].push_back(side->node_id(n));
            // Se algum boundary ID totalizar o número de
            // nós total da face, a face está no contorno
            if ( input_boundary_nodes[*bid_it].size() == side->n_nodes() )
            {
              // Todos os nós deste lado estão no contorno *bid_it
              for (std::vector<dof_id_type>::iterator node_it = input_boundary_nodes[*bid_it].begin();
                   node_it != input_boundary_nodes[*bid_it].end();
                   ++node_it)
              {

                // Ponteiro para o nó
                const Node * nptr = _in_mesh->node_ptr(*node_it);

                // Coordenadas do nó
                libMesh::Point p = *nptr;

                // Inserir ponto na triangulação
                if ( ! cgal_added_nodes[*node_it] )
                {
                  Delaunay::Vertex_handle vh = T.insert(Delaunay::Point(p(0),p(1),p(2)));
                  // Id do nó na malha libmesh original + 1 (evita usar 0 para identificação)
                  vh->info().id = *node_it + 1;
                  // Nó da malha original
                  vh->info().source = mesh3D::ORIGINAL;
                  // Marcar nó como incluído
                  cgal_added_nodes[*node_it] = true;
                }

              }
            }
          }
        }
      }
      else
      {
        // Face não está no contorno
        std::unique_ptr<Elem> side = elem->side_ptr(s);

        // Percorrer cada nó na face
        for (auto n : side->node_index_range())
        {
          // Ponteiro para o nó
          const Node * nptr = _in_mesh->node_ptr(side->node_id(n));

          // Coordenadas do nó
          libMesh::Point p = *nptr;

          // Inserir ponto na triangulação
          if ( ! cgal_added_nodes[side->node_id(n)] )
          {
            Delaunay::Vertex_handle vh = T.insert(Delaunay::Point(p(0),p(1),p(2)));
            // Id do nó na malha libmesh original + 1 (evita usar 0 para identificação)
            vh->info().id = side->node_id(n) + 1;
            // Nó da malha original
            vh->info().source = mesh3D::ORIGINAL;
            // Marcar nó como incluído
            cgal_added_nodes[side->node_id(n)] = true;
          }
        }
      }


  // Mapa com cada nó adicional e seus contornos correspondentes
  std::map<unsigned long int, std::vector<boundary_id_type>> added_node_boundaries;

  // Adicionar nós na zona de interesse
  unsigned int seed;
  getrandom((void *)&seed, sizeof(unsigned int), 0);
  srand(seed);

  /*
  // Pontos em MIN_Z
  for (unsigned long int added_node_id = 0 ; added_node_id < 50; added_node_id++)
  {
    Delaunay::Vertex_handle vh = T.insert(Delaunay::Point(-.1 + .2*RAND, -.1 + .2*RAND, -.05));
    // Id do nó na malha libmesh original + 1 (evita usar 0 para identificação)
    vh->info().id = added_node_id + 1;
    // Nó adicional
    vh->info().source = mesh3D::ADDED;
    // Associar ao contorno correspondente
    added_node_boundaries[added_node_id].push_back(BOUNDARY_ID_MIN_Z);
  }

  // Pontos em MAX_Z
  for (unsigned long int added_node_id = 0 ; added_node_id < 50; added_node_id++)
  {
    Delaunay::Vertex_handle vh = T.insert(Delaunay::Point(-.1 + .2*RAND, -.1 + .2*RAND, .05));
    // Id do nó na malha libmesh original + 1 (evita usar 0 para identificação)
    vh->info().id = added_node_id + 1;
    // Nó adicional
    vh->info().source = mesh3D::ADDED;
    // Associar ao contorno correspondente
    added_node_boundaries[added_node_id].push_back(BOUNDARY_ID_MAX_Z);
    // Se próximo à área de impacto, associar ao contorno apropriado
    if ( sqrt(pow(vh->point().x(), 2.)+pow(vh->point().y(), 2.)) < .1 )
      added_node_boundaries[added_node_id].push_back(PUSH_BOUNDARY_ID);
  }
  */

  assert(T.is_valid());

	// Nova malha libmesh
  _out_mesh->clear();
  // Preservar dimensão da malha
	_out_mesh->set_mesh_dimension(_in_mesh->mesh_dimension());

  _out_mesh->reserve_nodes(T.number_of_vertices());
  _out_mesh->reserve_elem(T.number_of_finite_cells());

  // Nós inseridos na nova malha libmesh
  std::map<Delaunay::Vertex_handle, bool> mesh_added_nodes;

  // Mapa vertex_handle -> novo node id
  std::map<Delaunay::Vertex_handle, dof_id_type> cgal_vertex_to_new_node;

  // Mapa com cada nó da nova malha libmesh e seus contornos correspondentes
  std::map<dof_id_type, std::vector<boundary_id_type>> node_boundaries;

  // Número sequencial de identificação do novo elemento
	unsigned int elem_id = 0;

  // [CGAL] Percorrer todas as células da malha gerada
  for (Delaunay::Finite_cells_iterator cit = T.finite_cells_begin(); cit != T.finite_cells_end(); cit++)
  {

		// Novo elemento na malha libmesh nova
		libMesh::Elem * elem = new libMesh::Tet4;
		elem->set_id(elem_id++);

		// Percorre os 4 pontos do tetraedro e insere os nós correspondentes
		for (unsigned int n = 0; n < elem->n_nodes(); n++)
		{

      // Nó a inserir na nova malha
      dof_id_type new_node_id;
      Node *new_node;

      if ( ! mesh_added_nodes[cit->vertex(n)] )
      {
        // Ausente na malha nova, inserir nó
				new_node = _out_mesh->add_point(libMesh::Point(
					cit->vertex(n)->point().x(), 
					cit->vertex(n)->point().y(), 
					cit->vertex(n)->point().z()
				));
        new_node_id = new_node->id();
        // Marcar como inserido no mapa
        mesh_added_nodes[cit->vertex(n)] = true;
        // Associar ao nó inserido
        cgal_vertex_to_new_node[cit->vertex(n)] = new_node_id;

        // Se vertex refere-se a um nó de contorno, recuperar seus boundary IDs
        if ( cit->vertex(n)->info().id > 0 )
        {
          // Associar novo node -> boundary ids
          if ( cit->vertex(n)->info().source == mesh3D::ORIGINAL )
          {
            // Usar mapeamento dos nós originais
            node_boundaries[new_node_id] = input_node_boundaries[cit->vertex(n)->info().id - 1];
          }
          else // mesh3D::ADDED
          {
            // Usar mapeamento dos nós adicionais
            node_boundaries[new_node_id] = added_node_boundaries[cit->vertex(n)->info().id - 1];
          }
        }

      }
      else
      {
        // Nó já foi inserido
        new_node_id = cgal_vertex_to_new_node[cit->vertex(n)];
      }

			// Definir n-ésimo nó do elemento
			elem->set_node(n) = _out_mesh->node_ptr(new_node_id);
		}

		// Incluir elemento na malha nova
		elem = _out_mesh->add_elem(elem);

  }

  // Contorno a reconstruir
  BoundaryInfo & out_boundary_info = _out_mesh->get_boundary_info();
  out_boundary_info.clear();

  // Construir Boundary info da malha nova
  for (auto & elem : _out_mesh->element_ptr_range())
    for (auto s : elem->side_index_range())
      if (elem->neighbor_ptr(s) == nullptr)
      {
        // Em cada face de elemento de contorno, gerar um mapa 
        // de contornos e seus nós correspondentes
        std::map<boundary_id_type, std::vector<dof_id_type>> boundary_nodes;

        // Face do elemento
        std::unique_ptr<Elem> side = elem->side_ptr(s);

        // Percorrer cada nó na face
        for (auto n : side->node_index_range())
        {
          // Identificar boundary IDs associados ao nó
          for (std::vector<boundary_id_type>::iterator bid_it = node_boundaries[side->node_id(n)].begin();
               bid_it != node_boundaries[side->node_id(n)].end();
               ++bid_it)
          {
            // Incluir nó no contador da boundary ID
            boundary_nodes[*bid_it].push_back(side->node_id(n));
            // Se algum boundary ID totalizar o número de
            // nós total da face, a face está no contorno
            if ( boundary_nodes[*bid_it].size() == side->n_nodes() )
            {
              // Todos os nós deste lado estão no contorno *bid_it
              out_boundary_info.add_side(elem, s, *bid_it);
              // Inserir nós no contorno
              for (auto m : side->node_index_range())
              {
                out_boundary_info.add_node(side->node_ptr(m), *bid_it);
              }
            }
          }
        }
      }

  // Antes, transformar malha para 2a. ordem.
  
  _out_mesh->prepare_for_use();

#ifdef DEBUG
  _out_mesh->print_info();
  out_boundary_info.print_summary();
#endif

}
