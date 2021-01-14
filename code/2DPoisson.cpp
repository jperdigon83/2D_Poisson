#include <iostream>
#include <cmath>
#include <fstream>

class TwoDPoisson{

public:
  // Constructeur

  TwoDPoisson(const std::size_t m){
    
    m_ = m;
    h_ = 1./(m_+1);
    

    lambda_ = new double[m_];
    
    u_ = new double*[m_];
    uTilde_ = new double*[m_];
    
    f_ = new double*[m_];
    fTilde_ = new double*[m_];
    
    s_ = new double*[m_];

    buffer_ = new double*[m_];
  
    for(int i=0; i<m_; ++i){
      
      u_[i] = new double[m_];
      uTilde_[i] = new double[m_];
      
      f_[i] = new double[m_];
      fTilde_[i] = new double[m_];
      
      s_[i] = new double[m_];

      buffer_[i] = new double[m_];
    }

    // Remplissage des matrices F, S et vecteur lambda
    
    for(std::size_t i=0; i<m_; ++i){
      
      lambda_[i] = 4*std::pow(sin(0.5*(i+1)*M_PI*h_),2);
      
      for(std::size_t j=0; j<m_; ++j){
	s_[i][j] = sin((i+1)*(j+1)*M_PI*h_);
	f_[i][j] = 1.;
      }
    }

    // Produit matriciel buffer = F * S

    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){

	buffer_[i][j] = 0;
  	for(std::size_t k=0; k<m_; ++k){
  	  buffer_[i][j] += f_[i][k]*s_[k][j];
  	}
      }
    }

    // Produit matriciel fTilde = S * buffer

    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){

	fTilde_[i][j] = 0;
  	for(std::size_t k=0; k<m_; ++k){
  	  fTilde_[i][j] += s_[i][k]*buffer_[k][j];
  	}
      }
    }

    // Produit matrice*Matrice = UTilde_ij = 4*h^4 * fTilde_ij / (lambda_i + lambda_j)
   
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	uTilde_[i][j] = 4*std::pow(h_,4)*fTilde_[i][j] / (lambda_[i] + lambda_[j]);	
      }
    }


    // Produit matriciel buffer =  utilde * S

    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){

	buffer_[i][j] = 0;
  	for(std::size_t k=0; k<m_; ++k){
  	  buffer_[i][j] += uTilde_[i][k]*s_[k][j];
  	}
      }
    }

    // Produit matriciel U = S * buffer

    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){

	u_[i][j] = 0;
  	for(std::size_t k=0; k<m_; ++k){
  	  u_[i][j] += s_[i][k]*buffer_[k][j];
  	}
      }
    }

    
  }


  // Destructeur
  
  ~TwoDPoisson(){
    
    delete [] lambda_;

    for(int i=0; i<m_; ++i){
      
      delete [] u_[i];
      delete [] uTilde_[i];

      delete [] f_[i];
      delete [] fTilde_[i];
    
      delete [] s_[i];

      delete [] buffer_[i];
    }

    delete [] u_;
    delete [] uTilde_;

    delete [] f_;
    delete [] fTilde_;
  
    delete [] s_;

    delete [] buffer_;
    
  }

  // Fonctions de classe

  void displayF() const{
    
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	
	std::cout << f_[i][j] << " "; 
      }

      std::cout << std::endl; 
    }
  }

  void displayS() const{
    
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	
	std::cout << s_[i][j] << " "; 
      }

      std::cout << std::endl; 
    }
  }

  void displayUTilde() const{
    
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	
	std::cout << uTilde_[i][j] << " "; 
      }

      std::cout << std::endl; 
    }
  }

  void displayU() const{
    
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	
	std::cout << u_[i][j] << " "; 
      }

      std::cout << std::endl; 
    }
  }

  void saveU(const char* filename) const{

    std::ofstream file(filename ,std::ios::trunc);
    
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	
	file << u_[i][j] << std::endl;; 
      }
    }
    
    file.close();
  }
  
private:

  std::size_t m_;
  double h_;

  double * lambda_;
  
  double ** u_;
  double ** uTilde_;
  
  double ** f_;
  double ** fTilde_;
 
  double ** s_;

  double ** buffer_;

  // Fonctions de classes (privées);

};



int main(){

  const std::size_t n = 128;

  TwoDPoisson test(n);

  const char* filename = "sol_f=1.dat";
  test.saveU(filename);
  
  
}
