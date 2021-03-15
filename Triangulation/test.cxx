#include "triangulation.hxx"

// libMesh namespace
using namespace libMesh;

int main(int argc, char **argv) {

  // Iniciar libMesh e bibliotecas subjacentes (PETSc)
  LibMeshInit init (argc, argv);

  // Instanciar uma nova malha
  Mesh mesh_hull(init.comm());

  // Gerar uma malha uniforme num cubo de lado 2 e aresta com 15 elementos
  Real halfside = 1.;
  MeshTools::Generation::build_cube (mesh_hull, 15, 15, 15, -halfside, halfside, -halfside, halfside, -halfside, halfside, TET4);
  //mesh.read ("3D.off");

  triangulate(mesh_hull);

  return 0;
}

