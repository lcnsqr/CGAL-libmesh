#include "triangulation.hxx"

// libMesh namespace
using namespace libMesh;

int main(int argc, char **argv) {

  // Iniciar libMesh e bibliotecas subjacentes (PETSc)
  LibMeshInit init (argc, argv);

  // Malha usada como casco e que também 
  // fornece os IDs de contorno e nós associados
  Mesh mesh_hull(init.comm());

  // Gerar uma malha uniforme num cubo de lado 2 e aresta com 15 elementos.
  // Cada face do cubo resultante é associada a um boundary ID, do 0 ao 5.
  Real halfside = 1.;
  int edgediv = 15;
  MeshTools::Generation::build_cube (mesh_hull, edgediv, edgediv, edgediv, -halfside, halfside, -halfside, halfside, -halfside, halfside, TET4);

  //mesh_hull.read ("3D.off");

  mesh_hull.write("antes.e");

  // Malha de destino
  Mesh mesh(init.comm());

  // Processar a nova malha
  mesh3D::Triangulation trng(&mesh_hull, &mesh);
  trng.remesh();

  mesh.write("depois.e");

  return 0;
}

