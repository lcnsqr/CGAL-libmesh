#include "convmesh.h"
#include <vector>

using namespace convmesh;

int main(int argc, char *argv[])
{
	std::vector<Triangulo<Ponto<double>>> casco;

  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,-1,-1),Ponto<double>(0,0,-1),Ponto<double>(-1,1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,1,-1),Ponto<double>(0,0,-1),Ponto<double>(1,1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,1,-1),Ponto<double>(0,0,-1),Ponto<double>(1,-1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,-1,-1),Ponto<double>(0,0,-1),Ponto<double>(-1,-1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,-1,-1),Ponto<double>(0,-1,0),Ponto<double>(1,-1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,-1,-1),Ponto<double>(0,-1,0),Ponto<double>(1,-1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,-1,1),Ponto<double>(0,-1,0),Ponto<double>(-1,-1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,-1,1),Ponto<double>(0,-1,0),Ponto<double>(-1,-1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,-1,-1),Ponto<double>(1,0,0),Ponto<double>(1,1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,1,-1),Ponto<double>(1,0,0),Ponto<double>(1,1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,1,1),Ponto<double>(1,0,0),Ponto<double>(1,-1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,-1,1),Ponto<double>(1,0,0),Ponto<double>(1,-1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,1,-1),Ponto<double>(0,1,0),Ponto<double>(-1,1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,1,-1),Ponto<double>(0,1,0),Ponto<double>(-1,1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,1,1),Ponto<double>(0,1,0),Ponto<double>(1,1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,1,1),Ponto<double>(0,1,0),Ponto<double>(1,1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,1,-1),Ponto<double>(-1,0,0),Ponto<double>(-1,-1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,-1,-1),Ponto<double>(-1,0,0),Ponto<double>(-1,-1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,-1,1),Ponto<double>(-1,0,0),Ponto<double>(-1,1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,1,1),Ponto<double>(-1,0,0),Ponto<double>(-1,1,-1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,-1,1),Ponto<double>(0,0,1),Ponto<double>(1,-1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,-1,1),Ponto<double>(0,0,1),Ponto<double>(1,1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(1,1,1),Ponto<double>(0,0,1),Ponto<double>(-1,1,1)));
  casco.push_back(Triangulo<Ponto<double>>(Ponto<double>(-1,1,1),Ponto<double>(0,0,1),Ponto<double>(-1,-1,1)));

	mesh_hull(casco);
}
