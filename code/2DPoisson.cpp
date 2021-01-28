#include <iostream>
#include <cmath>
#include <fstream>



inline double ** FFT(std::size_t N, double * x){

  double ** z = new double*[2]; // Complex vector
  z[0] = new double[N]; // real part
  z[1] = new double[N]; // Imaginary part

  if(N==1){
    z[0][0] = x[0];
    z[1][0] = 0;
    
  }else{
    
    double x_even[N/2];
    double x_odd[N/2];
    
    for(std::size_t i =0; i<(N/2); ++i){
      x_even[i] = x[2*i];
      x_odd[i] = x[2*i+1];
    }

    double ** z_even = FFT( N/2, x_even);
    double ** z_odd = FFT( N/2, x_odd);

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

inline double utheorique(int K, double x, double y, double Lx, double Ly){

  double u = 0;
  
  for(int k=0; k<K; ++k){
    u += sin((2*k+1)*M_PI*x/Lx)*( sinh((2*k+1)*M_PI*y/Ly) + sinh((2*k+1)*M_PI*(1-y/Ly)) )/(std::pow(2*k+1,3)*sinh(M_PI*(2*k+1)));
  }
  u = 0.5*(1-std::pow(2*x/Lx-1,2)) - 16*u/std::pow(M_PI,3);
  
  return u; 
}

class TwoDPoisson{

public:
  
  // Constructeur

  TwoDPoisson(const std::size_t N, const double Lx, const double Ly){

    m_ = N - 1;

    u_ = new double*[m_];  
    for(int i=0; i<m_; ++i){   
      u_[i] = new double[m_];     
    }

    double h = 1./N;
    double lambda[m_];
    
    for(std::size_t i=0; i<m_; ++i){
      lambda[i] = 4*std::pow(sin(0.5*(i+1)*M_PI*h),2);
      for(std::size_t j=0; j<m_; ++j){
	u_[i][j] = 1.0; 
      }
    }

    // Produit matriciel FS = F * S
    
    dst_col(u_);
    
    // Produit matriciel  Ftilde = SFS = S * FS

    dst_row(u_);
    
    // Produit Matrice*Matrice = UTilde_ij = 4*h^4 * FTilde_ij / (lambda_i + lambda_j)

    double Lx2 = Lx*Lx;
    double Ly2 = Ly*Ly;
    double rhs_factor_ = 4 * std::pow(h,4)* Lx2 * Ly2;
    
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	u_[i][j] = rhs_factor_ * u_[i][j] / (Ly2* lambda[i] + Lx2 * lambda[j]);
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

    delete [] u_;   
  }

  // Fonctions
  
  double ** get_U(){
    return u_;
  }
  
private:

  // Variables
  
  std::size_t m_;
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

int main(int argc, char** argv){
  
  const std::size_t N = atoi(argv[1]);
  const std::size_t Lx = atof(argv[2]);
  const std::size_t Ly = atof(argv[3]);
  
  TwoDPoisson test(N, Lx, Ly);

  double ** sol = test.get_U();

  
  // for(int i=0; i<N-1; ++i){  
  //   delete [] sol[i];
  // }

  // delete [] sol;   

  
  
  
  return 0; 
}
