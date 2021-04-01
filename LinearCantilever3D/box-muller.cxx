#include "box-muller.hxx"

// Um par de variáveis aleatórias com distribuição normal,
// esperança zero e variância unitária (método Box-Muller)
void boxMuller::parNormal(double *par){
	// Gerar um ponto uniforme dentro do disco 
	// de raio unitário centrado na origem
	double r;
	par[0] = -1.0 + 2.0*((double)rand()/RAND_MAX);
	par[1] = -1.0 + 2.0*((double)rand()/RAND_MAX);
	r = pow(par[0],2) + pow(par[1],2);
	// Gerar novamente caso ponto fora 
	// do disco ou exatamente na origem
	while ( r > 1 || r == 0 ){
		par[0] = -1.0 + 2.0*((double)rand()/RAND_MAX);
		par[1] = -1.0 + 2.0*((double)rand()/RAND_MAX);
		r = pow(par[0],2) + pow(par[1],2);
	}
	// Fator comum
	double f = sqrt(-2*log(r)/r);
	// Par normal
	par[0] *= f;
	par[1] *= f;
}

// n variáveis aleatórias com distribuição 
// normal, esperança mu e variância sigma
void boxMuller::normal(double *p, double mu, double sigma, int n){
	double par[2];
	if ( n == 1 ){
		parNormal(par);
		// Um valor do par é descartado
		p[0] = mu + sigma * par[0];
		return;
	}
	for (int i = 0; i < n - 1; i += 2){
		parNormal(par);
		p[i] = mu + sigma * par[0];
		p[i+1] = mu + sigma * par[1];
	}
	// Se n ímpar, falta gerar o último
	if ( n % 2 == 1 ){
		parNormal(par);
		// Um valor do par é descartado
		p[n-1] = mu + sigma * par[0];
	}
}
