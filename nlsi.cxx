#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);
void step(cmplx* const psi1, cmplx* const psi0,
          const double dt, const double dx,
          const int N);
void step2(cmplx* const psiteeeeh,
	  cmplx* const psi0,
          const double dt,
          const int N);

//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];

	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin);


	for (int i = 1; i <= Na; i++) {

		for (int j = 1; j <= Nk-1; j++) {
		step(psi1,psi0,dt,dx,Nx);
		step2(psi1,psi0,dt,Nx);  
		}
		
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}

	return 0;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}
void step(cmplx* const psi1, cmplx* const psi0,
          const double dt, const double dx,
          const int N)
{

  cmplx* d=new cmplx[N];
  cmplx* u=new cmplx[N];
  cmplx* l=new cmplx[N];

  for(int i=0;i<N;i++) d[i] = 1.0+2.0*cmplx(0.0,-dt/(dx*dx));
  for(int i=0;i<N;i++) u[i] = cmplx(0.0,dt/(dx*dx));
  for(int i=0;i<N;i++) l[i] = cmplx(0.0,dt/(dx*dx));

  for(int j = 0; j < N-1; j++){
    d[j+1] -= u[j] * l[j+1] / d[j];
    psi0[j+1]-= psi0[j] * l[j+1] / d[j];  
  }
  
  psi1[N-1] = psi0[N-1] / d[N-1];
  
  for(int j = N-2; j >= 0; j--){
    psi1[j] = (psi0[j] - u[j] * psi1[j+1]) / d[j];
  }

  delete[] d;
  delete[] u;
  delete[] l;
}

void step2(cmplx* const psiteeeeh,
	  cmplx* const psi0,
          const double dt,
          const int N)
{
  double* normquad = new double[N];
  for(int i = 0; i<N; i++){
    normquad[i]=abs(psiteeeeh[i])*abs(psiteeeeh[i]);  
    psi0[i]=psiteeeeh[i]*exp(cmplx(0.0,-normquad[i]*dt));
  }
  delete[] normquad;
}