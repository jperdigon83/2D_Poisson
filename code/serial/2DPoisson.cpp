#include <iostream>
#include <cmath>
#include <fstream>
#include<chrono>


/* Fonction récursive pour calculer la FFT d'un vecteur x (algorithme de Cooley-Tukey).
 */
inline double ** FFT(std::size_t N, double * x){

  double ** z = new double*[2]; // Vecteur complexe
  z[0] = new double[N]; // partie réelle
  z[1] = new double[N]; // partie immaginaire

  if(N==1){
    z[0][0] = x[0];
    z[1][0] = 0;
    
  }else{
    
    double x_even[N/2]; // Eléments pairs 
    double x_odd[N/2]; // Eléments impairs
    
    for(std::size_t i =0; i<(N/2); ++i){
      x_even[i] = x[2*i];
      x_odd[i] = x[2*i+1];
    }

    double ** z_even = FFT( N/2, x_even); // FFT des éléments pairs
    double ** z_odd = FFT( N/2, x_odd); // FFT des éléments imparis

    for(std::size_t k =0; k<(N/2); ++k){
      
      z[0][k] = z_even[0][k] + cos(2*M_PI*k/N)*z_odd[0][k] + sin(2*M_PI*k/N)*z_odd[1][k];
      z[1][k] = z_even[1][k] + cos(2*M_PI*k/N)*z_odd[1][k] - sin(2*M_PI*k/N)*z_odd[0][k];

      z[0][k+N/2] = z_even[0][k] - cos(2*M_PI*k/N)*z_odd[0][k] - sin(2*M_PI*k/N)*z_odd[1][k];
      z[1][k+N/2] = z_even[1][k] - cos(2*M_PI*k/N)*z_odd[1][k] + sin(2*M_PI*k/N)*z_odd[0][k];
      
    }

    delete [] z_even[0];
    delete [] z_even[1];   
    delete [] z_even;

    delete [] z_odd[0];
    delete [] z_odd[1];   
    delete [] z_odd;  
  }

  return z;
  
}

/* Fonction pour calculer la DST d'un vecteur x à partir de la FFT.
 */

inline void DST(std::size_t N, double * x){

  int Nl = 2*N+2;
  
  double xl[Nl];
  
  xl[0] = 0;
  xl[N+1] = 0;
  for(std::size_t i=1; i<N+1; ++i){
    xl[i] = x[i-1];
    xl[N+1+i] = -x[N-i];
  }
  
  double ** z = FFT(Nl,xl);

  for(std::size_t i=0; i<N; ++i){
    x[i] = -0.5*z[1][i+1];
  }
  
  delete [] z[0];
  delete [] z[1];   
  delete [] z;

}

/* Classe pour calculer la solution du problème de Poisson à 2D sur un domaine rectangulaire [0,Lx] x [0,Ly]
 */

class TwoDPoisson{

public:
  
  // Constructeur

  TwoDPoisson(const std::size_t N, const double Lx, const double Ly){

    m_ = N - 1;
    h_ = 1./N;
    
    x_ = new double[m_];
    y_ = new double[m_];
    lambda_ = new double[m_];
    u_ = new double*[m_];
    
    for(std::size_t i=0; i<m_; ++i){

      x_[i] = (i+1)*h_*Lx;
      y_[i] = (i+1)*h_*Ly;
      lambda_[i] = 4*std::pow(sin(0.5*(i+1)*M_PI*h_),2);
      u_[i] = new double[m_];
      
      for(std::size_t j=0; j<m_; ++j){
	u_[i][j] = 1.0; 
      }
    }

    // Produit matriciel: SF = S * F
    
    dst_col(u_);
    
    // Produit matriciel: Ftilde = SFS = SF * S

    dst_row(u_);
    
    // Produit Matrice*Matrice = UTilde_ij =  4*h^4 * Lx**2 * Ly**2 *FTilde_ij / (lambda_i + lambda_j)

    double Lx2 = Lx*Lx;
    double Ly2 = Ly*Ly;
    double rhs_factor_ = 4 * std::pow(h_,4)* Lx2 * Ly2;
    
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	u_[i][j] = rhs_factor_ * u_[i][j] / (Ly2* lambda_[i] + Lx2 * lambda_[j]);
      }
    }
    
    // Produit matriciel UtildeS  =  Utilde * S
    
    dst_col(u_);

    // Produit matriciel U = S * UtildeS
    
    dst_row(u_);
    
  }
  
  // Destructeur
  
  ~TwoDPoisson(){

    for(int i=0; i<m_; ++i){  
      delete [] u_[i];
    }

    delete [] x_;
    delete [] y_;
    delete [] lambda_;
    delete [] u_;   
  }

  // Fonctions

  double getU(const std::size_t i, const std::size_t j){
    return u_[i][j];
  }

  double getx(const std::size_t i){
    return x_[i];
  }

  double gety(const std::size_t j){
    return y_[j];
  }
  
  void saveU(const char* filename) const{

    std::ofstream file(filename ,std::ios::trunc);

    file << "#x y u" << std::endl;
  
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){

	file << x_[i] << " " << y_[j] << " " << u_[i][j] << std::endl;
      }
    }
    
    file.close();    
  }
  
private:

  // Variables
  
  std::size_t m_;
  double h_;

  double * x_;
  double * y_;
  double * lambda_;

  double ** u_;

  // Fonctions
  
  void dst_col(double ** A) const{
    
    for(std::size_t j=0; j<m_; ++j){
       
      double x[m_];
      for(std::size_t i=0; i<m_; ++i){
	x[i] = A[i][j];
      }
       
      DST(m_,x);

      for(std::size_t i=0; i<m_; ++i){
	A[i][j] = x[i];	 
      }
    }
  }
  
  void dst_row(double ** A) const{
    
    for(std::size_t i=0; i<m_; ++i){
       
      double x[m_];
      
      for(std::size_t j=0; j<m_; ++j){
	x[j] = A[i][j];
      }
       
      DST(m_,x);

      for(std::size_t j=0; j<m_; ++j){
	A[i][j] = x[j];	 
      }   
    }
  }
  
};




/* Solution analytique de l'équation de Poisson  avec  f=-1,  sur le domaine [0,PI] x[0,PI], avec conditions aux bord de Dirichelt u(0,y) = u(x,0) = u(PI,y) = u(x,PI) = 0. Comme la solution est sous la forme d'une série, la solution est tronquée aux termes (M,N);
 */

inline double Utheorique(const std::size_t M, const std::size_t N, const double x, const double y){

  double u = 0;
  double coeff;
  
  for(std::size_t m=1; m<=M; ++m){
    for(std::size_t n=1; n<=N; ++n){
      coeff = 4*(1-std::pow(-1,m))*(1-std::pow(-1,n))/((std::pow(m,2)+std::pow(n,2))*m*n*std::pow(M_PI,2));
      u += coeff * sin(m*x) * sin(n*y);
    }
  }

  return u;  
}



int main(int argc, char** argv){

  const double  Lx = M_PI;
  const double  Ly = M_PI;

  /* // Measure and save in a file the serial computation time as a function of N, the poisson equation with constant f=-1.
 
  
     std::ofstream file("time_serial.dat" ,std::ios::trunc);
     file << "#N t(s)" << std::endl;

     for(std::size_t k=3; k<=9; ++k){

     const std::size_t N = std::pow(2,k);

     std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

     TwoDPoisson U(N, Lx, Ly);
    
     std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

     auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

     file << N  << " " << duration/1e6 << std::endl;
    
     }
     file.close();

  */

  /* // This part measures the serial truncature error as a function of N, for the poisson equation with constant f=-1. 
   
  
     std::ofstream file("serial_error.dat" ,std::ios::trunc);
     file << "#N errL2" << std::endl;
   
     for(std::size_t k=3; k<=9; ++k){
    
     const std::size_t N = std::pow(2,k);
    
     TwoDPoisson U(N, Lx, Ly);
   
     const std::size_t M = N-1;
     double L2err = 0;
     for(std::size_t i=0; i<M; ++i){
     for(std::size_t j=0; j<M; ++j){
     L2err += std::pow(U.getU(i,j) - Utheorique(200,200,U.getx(i),U.gety(j)),2);
     }
     }
     L2err = sqrt( L2err / (M*M) );

     std::cout << "N = " << N << " L2err = " <<  L2err << std::endl;
     file << N  << " " << L2err  << std::endl;

     }

     file.close();

  */
  
  return 0; 
}
