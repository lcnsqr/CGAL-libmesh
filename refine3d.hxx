#include <string.h>
#include <array>
#include <vector>

namespace refine3d {

  template <typename T> class Ponto {

    public:
		std::array<T, 3> pos;

		Ponto();
		Ponto(const T, const T, const T);

  };

  template <typename T> Ponto<T>::Ponto(){
  }

  template <typename T> Ponto<T>::Ponto(const T x, const T y, const T z){
		this->pos[0] = x;
		this->pos[1] = y;
		this->pos[2] = z;
  }

	// Triângulo que compõe o casco da malha
  template <typename T> class Triangulo {
    public:
		std::array<T, 3> vertices;

		Triangulo();
		Triangulo(const T, const T, const T);
	};

  template <typename T> Triangulo<T>::Triangulo(){
  }

  template <typename T> Triangulo<T>::Triangulo(const T a, const T b, const T c){
		this->vertices[0] = a;
		this->vertices[1] = b;
		this->vertices[2] = c;
  }


	void mesh_hull(std::vector<Triangulo<Ponto<double>>>);
}
