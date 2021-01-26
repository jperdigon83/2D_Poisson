#include <iostream>
#include <cmath>
#include <fstream>

class TwoDPoisson{

public:
  
  // Constructeur

  TwoDPoisson(const std::size_t m, const double Lx, const double Ly){
    
    m_ = m;
    h_ = 1./(m_+1);
    Lx2_ = Lx*Lx;
    Ly2_ = Ly*Ly;
    rhs_factor_ = 4 * std::pow(h_,4)* Lx2_ * Ly2_;

    x_ = new double[m_];
    y_ = new double[m_]; 
    lambda_ = new double[m_];
   
    u_ = new double*[m_];
    f_ = new double*[m_];
    buff_ = new double*[m_];
    s_ = new double*[m_];
    
    for(int i=0; i<m_; ++i){
      
      u_[i] = new double[m_]; 
      f_[i] = new double[m_];
      buff_[i] = new double[m_];
      s_[i] = new double[m_];
     
    }

    // Remplissage des matrices F, S et vecteur lambda
    
    for(std::size_t i=0; i<m_; ++i){
      
      x_[i] = (i+1)*h_*Lx;
      y_[i] = (i+1)*h_*Ly;
      lambda_[i] = 4*std::pow(sin(0.5*(i+1)*M_PI*h_),2);
      
      for(std::size_t j=0; j<m_; ++j){
	s_[i][j] = sin((i+1)*(j+1)*M_PI*h_);
	f_[i][j] = 1.; // Homogenous Boundary conditions with constant f
      }
    }

    // Produit matriciel FS = F * S
    
    matricialDot(u_, f_, s_);
    
    // Produit matriciel  Ftilde = SFS = S * FS
    
     matricialDot(buff_, s_, u_);

    // Produit Matrice*Matrice = UTilde_ij = 4*h^4 * FTilde_ij / (lambda_i + lambda_j)
   
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	u_[i][j] = rhs_factor_ * buff_[i][j] / (Ly2_* lambda_[i] + Lx2_ * lambda_[j]);	
      }
    }
    
    // Produit matriciel UtildeS  =  Utilde * S

    matricialDot(buff_, u_, s_);

    // Produit matriciel U = S * UtildeS

    matricialDot(u_, s_, buff_);
    
  }
  
  // Destructeur
  
  ~TwoDPoisson(){

    delete [] x_;
    delete [] y_;
    delete [] lambda_;

    for(int i=0; i<m_; ++i){
      
      delete [] u_[i];
      delete [] f_[i];
      delete [] buff_[i];
      delete [] s_[i];
   
    }

    delete [] u_;
    delete [] f_;
    delete [] buff_;
    delete [] s_;
    
  }

  // Fonctions de classe

  void display_U() const{ // affiche dans la console la solution sous sa forme matricielle
    
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	
	std::cout << u_[i][j] << " "; 
      }
      std::cout << std::endl; 
    }
  }

  double display_errL2(){ // retourne l'erreur L2 

    double errL2 = 0;
    
    for(std::size_t i=1; i<m_-1; ++i){
      for(std::size_t j=1; j<m_-1; ++j){
	errL2 +=  std::pow( 1./Lx2_*(-u_[i+1][j] + 2*u_[i][j] - u_[i-1][j]) + 1./Ly2_*(-u_[i][j+1] + 2*u_[i][j] - u_[i][j-1])- h_*h_*f_[i][j] , 2);
      }
    }
    return sqrt(h_*errL2);
  }

  void saveU(const char* filename) const{ // Enregiste la solution dans un fichier

    std::ofstream file(filename ,std::ios::trunc);

    
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	
	file << x_[i] << " " << y_[j] << " " << u_[i][j] << std::endl;
      }
    }
    
    file.close();
  }
  
private:

  std::size_t m_;
  double h_;
  double Lx2_;
  double Ly2_;
  double rhs_factor_;
 
  double * x_;
  double * y_;
  double * lambda_;
  
  double ** u_;
  double ** f_;
  double ** buff_;
  double ** s_;

  // Fonctions de classes (priv�es);

  void matricialDot(double ** &C, double ** &A, double ** &B) const{

    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){

	C[i][j] = 0;
	
  	for(std::size_t k=0; k<m_; ++k){
  	  C[i][j] += A[i][k]*B[k][j];
  	}
      }
    }  
  }
  
};

int main(int argc, char** argv){
  
  const std::size_t M = atoi(argv[1]);
  const std::size_t Lx = atof(argv[2]);
  const std::size_t Ly = atof(argv[3]);
  
  TwoDPoisson test(M, Lx, Ly);

  const char* filename = "../result.dat";
  test.saveU(filename);

  return 0; 
}
