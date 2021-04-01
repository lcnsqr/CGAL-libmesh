#include <stdlib.h>
#include <math.h>

namespace boxMuller {
  // Um par de variáveis aleatórias com distribuição normal,
  // esperança zero e variância unitária (método Box-Muller)
  void parNormal(double *par);

  // n variáveis aleatórias com distribuição 
  // normal, esperança mu e variância sigma
  void normal(double *p, double mu, double sigma, int n);
}
