#include <iostream>
#include <cmath>
#include <fstream>
#include<mpi.h>

/* Fonction récursive pour calculer la FFT d'un vecteur x (algorithme de Cooley-Tukey).
 */

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

/* Classe parallélisée pour calculer la solution du problème de Poisson à 2D sur un domaine rectangulaire [0,Lx] x [0,Ly]. Chaque processeur possède une matrice u_ de taille (m_*m_) mais calcule uniquement la DST sur une sous partie des lignes / colonnes. Entre chaque produit matriciel, il y a communication entre les processeurs afin de mettre à jour l'ensemble de la matrice, sur chaque processeur).
 */


class TwoDPoisson{

public:
  
  // Constructeur

  TwoDPoisson(const std::size_t N, const double Lx, const double Ly, const int procId, const int nProc){
    
    m_ = N - 1;
    h_ = 1./N;
    
    nProc_ = nProc;
    procId_ = procId;
    
    M_ = m_/nProc_;
    Mr_ = m_%nProc_;
       
    x_ = new double[m_];
    y_ = new double[m_];
    lambda_ = new double[m_];
    u_ = new double*[m_];
    
    for(std::size_t i=0; i<m_; ++i){

      x_[i] = (i+1)*h_*Lx;
      y_[i] = (i+1)*h_*Ly;
      lambda_[i] = 4*std::pow(sin(0.5*(i+1)*M_PI*h_),2);
      u_[i] = new double[m_];
    }

    
    for(std::size_t i=0; i<m_; ++i){ 
      for(std::size_t j=0; j<m_; ++j){	
	//u_[i][j] = 1.0; // Constant source
	u_[i][j] = exp(x_[i])*sin(2*M_PI*x_[i])*sin(2*M_PI*y_[j]);
	//u_[i][j] = sin(2*M_PI*x_[i])*sin(2*M_PI*y_[j]);

	// Point sources
	/*
	if(i==m_/4 && j==m_/4){
	  u_[i][j] = 1;
	}else if(i==3*m_/4 && j==3*m_/4){
	  u_[i][j] = -1;
	}else{u_[i][j] = 0;}
	*/
      }    
    }
    
    // Produit matriciel: SF = S * F
    
    dst_col(u_); 
    
    communication_col(u_); 
    
    MPI_Barrier ( MPI_COMM_WORLD );

    // Produit matriciel: Ftilde = SFS = SF * S

    dst_row(u_);
    
    communication_row(u_);
    
    MPI_Barrier ( MPI_COMM_WORLD );

    
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
   
    communication_col(u_);
    
    MPI_Barrier (MPI_COMM_WORLD );


    // Produit matriciel U = S * UtildeS
    
    dst_row(u_);
   
    communication_row(u_);
    
    MPI_Barrier ( MPI_COMM_WORLD );

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
  std::size_t procId_;
  std::size_t nProc_;
  std::size_t M_;
  std::size_t Mr_;
  
  double h_;

  double * x_;
  double * y_;
  double * lambda_;

  double ** u_;

  // Fonctions

  void communication_col(double ** A) const{

    int displs[nProc_];
    int recv_counts[nProc_];

    for (std::size_t i=0; i<nProc_-1; ++i) {
      displs[i] = i * M_ * m_;
      recv_counts[i] = M_ * m_;
    }
    
    displs[nProc_-1] = (nProc_-1) * M_ * m_;
    recv_counts[nProc_-1] = (M_+Mr_) * m_;

    double data[m_*m_];
 
    for(std::size_t j=0; j<m_; ++j){
      for(std::size_t i=0; i<m_; ++i){
	data[j*m_+i] = A[i][j];
      }
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &(data[0]), recv_counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    
    for(std::size_t j=0; j<m_; ++j){
      for(std::size_t i=0; i<m_; ++i){
	 A[i][j] = data[j*m_+i];
      }
    }
    
  }

  void communication_row(double ** A) const{
    
    int displs[nProc_];
    int recv_counts[nProc_];

    for (std::size_t i=0; i<nProc_-1; ++i) {
      displs[i] = i * M_ * m_;
      recv_counts[i] = M_ * m_;
    }
    
    displs[nProc_-1] = (nProc_-1) * M_ * m_;
    recv_counts[nProc_-1] = (M_+Mr_) * m_;

    double data[m_*m_];
 
    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	data[i*m_+j] = A[i][j];
      }
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &(data[0]), recv_counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    for(std::size_t i=0; i<m_; ++i){
      for(std::size_t j=0; j<m_; ++j){
	A[i][j] = data[i*m_+j];
      }
    }
    
  }
  // Fonction parallélisée pour le calcul de la DST des colonnes de la matrice. Les colonnes sont réparties parmi les processeurs et le reste est pris en charge par le dernier processeur.
  
  void dst_col(double ** A) const{

    if (procId_ != (nProc_ -1)){

      for(std::size_t j=M_*procId_; j<M_*(procId_+1); ++j){
       
      double x[m_];
      for(std::size_t i=0; i<m_; ++i){
	x[i] = A[i][j];
      }
       
      DST(m_,x);

      for(std::size_t i=0; i<m_; ++i){
	A[i][j] = x[i];	 
      }
    }
     
    }else{

      for(std::size_t j=(M_*procId_); j<M_*(procId_+1)+Mr_; ++j){
       
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
  
  }
  
  void dst_row(double ** A) const{

     if (procId_ != (nProc_ -1)){

       for(std::size_t i=(M_*procId_); i<M_*(procId_+1); ++i){
       
	 double x[m_];
      
	 for(std::size_t j=0; j<m_; ++j){
	   x[j] = A[i][j];
	 }
       
	 DST(m_,x);

	 for(std::size_t j=0; j<m_; ++j){
	   A[i][j] = x[j];	 
	 }   
       }
   
     }else{

       for(std::size_t i=(M_*procId_); i<M_*(procId_+1)+Mr_; ++i){
       
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
     
  }
  
};

/* Solution analytique de l'équation de Poisson  avec  f=-1,  sur le domaine [0,PI] x[0,PI], avec conditions aux bord de Dirichelt u(0,y) = u(x,0) = u(PI,y) = u(x,PI) = 0. Comme la solution est sous la forme d'une série, la solution est tronquée aux termes (M,N);

Comme amélioration notable, on peut ultérieurement paralléliser cette fonction ( si le temps le permet :) )
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

  MPI_Status status;
  int procId;
  int nProc;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);

  const double  N = atoi(argv[1]);
  const double  Lx = atof(argv[2]);
  const double  Ly = atof(argv[3]);

  TwoDPoisson U(N, Lx, Ly, procId, nProc);
  if(procId == 0){
    U.saveU("res.dat");
  }


  /*

  // This part measures the serial truncature error as a function of N, for the poisson equation with constant f=-1. 

  const double  Lx = M_PI;
  const double  Ly = M_PI;
   
  std::ofstream file("error_parallel_P=1.dat" ,std::ios::trunc);

  if(procId == 0){
  file << "#N errL2" << std::endl;
  }
   
  for(std::size_t k=3; k<=9; ++k){
    
  const std::size_t N = std::pow(2,k);
    
  TwoDPoisson U(N, Lx, Ly, procId, nProc);
   
  const std::size_t M = N-1;
  double L2err = 0;
  for(std::size_t i=0; i<M; ++i){
  for(std::size_t j=0; j<M; ++j){
  L2err += std::pow(U.getU(i,j) - Utheorique(200,200,U.getx(i),U.gety(j)),2);
  }
  }
  L2err = sqrt( L2err / (M*M) );

  if(procId == 0){
  std::cout << "N = " << N << " L2err = " <<  L2err << std::endl;
  file << N  << " " << L2err  << std::endl;
  }

  }

  file.close();

  */

  


  
  

  /* 

  // For a given number of procs, This part measures and save in a file the computation time as a function of N, for the poisson equation with constant f=-1. 
  

  std::ofstream file("time_parallel.dat" ,std::ios::trunc);

  if(procId == 0){
  file << "#N  t" << std::endl ; 
  }
   
  for(std::size_t k=3; k<=9; ++k){
    
  const std::size_t N = std::pow(2,k);

  double times,timee;
  times = MPI_Wtime();
     
  TwoDPoisson U(N, Lx, Ly, procId, nProc);

  timee = MPI_Wtime();
  double time = timee - times;

  double res;
  MPI_Reduce(&time, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(procId == 0){
  std::cout << "N = " << N << " Elapsed time = "  << res   << " s" << std::endl;
  file << N  << " " << res << std::endl;
  }   
  }

  file.close();

  */




  
  MPI_Finalize();
  
  return 0; 
}
