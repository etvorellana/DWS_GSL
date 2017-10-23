/*
 ============================================================================
 Name        : libdws.c
 Author      : Esbel T. Valero Orellana (evalero@uesc.br)
 Version     : 0.1
 Copyright   : Copyright 2017 Esbel Tomas Valero Orellana <evalero@uesc.br>
 	 	 	 	 This program is free software; you can redistribute it and/or
 	 	 	 	 modify  it under the terms of the GNU General Public License
 	 	 	 	 as published by  the Free Software Foundation; either version
 	 	 	 	 2 of the License, or  (at your option) any later version.

 	 	 	 	 This program is distributed in the hope that it will be useful,
 	 	 	 	 but WITHOUT ANY WARRANTY; without even the implied warranty of
 	 	 	 	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 	 	 	 	 GNU General Public License for more details.

 	 	 	 	 You should have received a copy of the GNU General Public
 	 	 	 	 License along with this program; if not, write to the Free
 	 	 	 	 Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 	 	 	 	 Boston, MA 02110-1301, USA.

 Description : TODO
 ============================================================================
 */

#include "dws.h"

/*
============================================================================
*/

gsl_matrix* setLinSpace(double T_0, double dt, double T_F, int *ls_size)
{
	*ls_size = (int) round((T_F - T_0) / dt) + 1;
	//double* ls = (double*) malloc(*ls_size * sizeof(double));
  gsl_matrix* ls = gsl_matrix_alloc (*ls_size, 1);
	for(int i = 0; i < *ls_size; i++)
		gsl_matrix_set(ls, i, 0, i*dt);
	return ls;
}

/*
============================================================================
*/

gsl_matrix* setGrid(int dim, double Ro, gsl_matrix* ThObs, int Nthe,
                                        gsl_matrix* Phi, int Nphi)
{
	gsl_matrix* grid = gsl_matrix_alloc(Nthe,Nphi);
  double ThObs_, Phi_;
	for(int i = 0; i < Nthe; i++){
		for(int j = 0; j < Nphi; j++){
      ThObs_  = gsl_matrix_get(ThObs, i, 0);
      Phi_ = gsl_matrix_get(Phi, j, 0);
			switch(dim){
			case 0: //Xobs=Ro*sin(ThObs)'*cos(Phi);
				gsl_matrix_set(grid, i, j, Ro*sin(ThObs_)*cos(Phi_));
				break;
			case 1: //Yobs=Ro*sin(ThObs)'*sin(Phi);
				gsl_matrix_set(grid, i, j, Ro*sin(ThObs_)*sin(Phi_));
				break;
			default: //Zobs=Ro*cos(ThObs)'*ones(1,length(Phi));
				gsl_matrix_set(grid, i, j, Ro*cos(ThObs_));
			}
		}
	}
	return grid;
}

/*
============================================================================
*/

gsl_matrix* normrnd(double mu, double sigma,int rows, int cols)
{
	gsl_matrix *dist = gsl_matrix_alloc(rows,cols);
	gsl_rng * r;
	// create random number generator
  r = gsl_rng_alloc (gsl_rng_mt19937);
	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			gsl_matrix_set(dist, i, j, mu + gsl_ran_gaussian(r, sigma));
	gsl_rng_free (r);
	return dist;
}

/*
============================================================================
*/

gsl_matrix* atanVet(gsl_matrix* Vd, double omTd, gsl_matrix* d, int N)
{
	gsl_matrix *phid = gsl_matrix_alloc(N, 1);
  double Vd_, d_;
	for(int i = 0; i < N; i++){
    Vd_ = gsl_matrix_get(Vd, i,0);
    d_ = gsl_matrix_get(d, i, 0);
		gsl_matrix_set(phid, i, 0, atan2(-1*Vd_/omTd, d_));
  }

	return phid;
}

/*
============================================================================
*/

gsl_matrix* setD0(gsl_matrix* D, gsl_matrix* PhiD, int N)
{
	gsl_matrix *D0 = gsl_matrix_alloc(N,1);
  double D_, PhiD_;
	for(int i = 0; i < N; i++){
    D_ = gsl_matrix_get(D, i, 0);
    PhiD_ = gsl_matrix_get(PhiD, i, 0);
		gsl_matrix_set(D0, i, 0, D_/cos(PhiD_));
  }
	return D0;
}

/*
============================================================================
*/

gsl_matrix* ones(int rows, int cols)
{
  gsl_matrix *matrix = gsl_matrix_alloc(rows, cols);
  gsl_matrix_set_all(matrix, 1.0);
  return matrix;
}



/*
============================================================================
*/

double complex Elaser(double x, double y, double z, double E0, double k)
{
	return E0*cexp(k*I*z);
}







// % %***************************************
// ? kernel3(gsl_matrix* Xt,gsl_matrix* Yt,gsl_matrix* Zt, int N)
// {
//   gsl_matrix* ones1N = ones(1,N);
//   gsl_matrix* Xjm = gsl_matrix_calloc(N,N);
//   gsl_matrix* Yjm = gsl_matrix_calloc(N,N);
//   gsl_matrix* Xjm = gsl_matrix_calloc(N,N);
//   gsl_matrix* Rjm = gsl_matrix_calloc(N,N);
//   //Xjm=Xt*ones(1,N)-ones(N,1)*Xt';
//   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Xt, ones1N, 0.0, Xjm);
//   gsl_blas_dgemm(CblasTrans, CblasTrans, -1.0, ones1N, Xt, 1.0, Xjm);
//   //Yjm=Yt*ones(1,N)-ones(N,1)*Yt';
//   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Yt, ones1N, 0.0, Yjm);
//   gsl_blas_dgemm(CblasTrans, CblasTrans, -1.0, ones1N, Yt, 1.0, Yjm);
//   //Zjm=Zt*ones(1,N)-ones(N,1)*Zt';
//   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Zt, ones1N, 0.0, Zjm);
//   gsl_blas_dgemm(CblasTrans, CblasTrans, -1.0, ones1N, Zt, 1.0, Zjm);
//   //Rjm=sqrt(Xjm.^2+Yjm.^2+Zjm.^2)+diag(ones(N,1));
//   double sum;
//   for(int i = 0; i < N; i++){
//     for(int j = 0; j < i; j++){
//       sum = gsl_matrix_get(Xjm, i, j)*gsl_matrix_get(Xjm, i, j);
//       sum += gsl_matrix_get(Yjm, i, j)*gsl_matrix_get(Yjm, i, j);
//       sum += gsl_matrix_get(Zjm, i, j)*gsl_matrix_get(Zjm, i, j);
//       gsl_matrix_set(Rjm, i, j, sqrt(sum);
//     }
//     sum = 1.0;
//     sum += gsl_matrix_get(Xjm, i, i)*gsl_matrix_get(Xjm, i, i);
//     sum += gsl_matrix_get(Yjm, i, i)*gsl_matrix_get(Yjm, i, i);
//     sum += gsl_matrix_get(Zjm, i, i)*gsl_matrix_get(Zjm, i, i);
//     gsl_matrix_set(Rjm, i, i, sqrt(sum));
//     for(int j = i+1; j < N; j++){
//       sum = gsl_matrix_get(Xjm, i, j)*gsl_matrix_get(Xjm, i, j);
//       sum += gsl_matrix_get(Yjm, i, j)*gsl_matrix_get(Yjm, i, j);
//       sum += gsl_matrix_get(Zjm, i, j)*gsl_matrix_get(Zjm, i, j);
//       gsl_matrix_set(Rjm, i, j, sqrt(sum));
//     }
//   }
//
//   gsl_matrix_complex* Ker = gsl_matrix_complex_calloc(3*N, 3*N);
//
//
//
//   gsl_matrix_complex_free(Ker);
//   gsl_matrix_free(ones1N);
//   gsl_matrix_free(Xjm);
//   gsl_matrix_free(Yjm);
//   gsl_matrix_free(Zjm);

//
//     c1=exp(1i*k*Rjm)./(1i*k*Rjm.^3);
//     c2=1i./(k*Rjm)-1./(k*Rjm).^2;
//
//     Ker(1:N,1:N)=c1.*(Rjm.^2-Xjm.^2+c2.*(Rjm.^2-3*Xjm.^2));
//     Ker(N+1:2*N,N+1:2*N)=c1.*(Rjm.^2-Yjm.^2+c2.*(Rjm.^2-3*Yjm.^2));
//     Ker(2*N+1:3*N,2*N+1:3*N)=c1.*(Rjm.^2-Zjm.^2+c2.*(Rjm.^2-3*Zjm.^2));
//     Ker(1:N,N+1:2*N)=c1.*(-Xjm.*Yjm-3*c2.*Xjm.*Yjm);
//     Ker(N+1:2*N,1:N)=Ker(1:N,N+1:2*N);
//     Ker(1:N,2*N+1:3*N)=c1.*(-Xjm.*Zjm-3*c2.*Xjm.*Zjm);
//     Ker(2*N+1:3*N,1:N)=Ker(1:N,2*N+1:3*N);
//     Ker(N+1:2*N,2*N+1:3*N)=c1.*(-Yjm.*Zjm-3*c2.*Yjm.*Zjm);
//     Ker(2*N+1:3*N,N+1:2*N)=Ker(N+1:2*N,2*N+1:3*N);
//
//     Ker=Ker-diag(diag(Ker)); % Removing artificially non-zero diagonal terms
//     kn=(-Gamma/2+1i*Delta)*eye(3*N)-Gamma/2*Ker;
// end
//
