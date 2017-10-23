/*
 ============================================================================
 Name        : DWS_C.c
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

int main(void) {

	// *** Physical constants ***
	double lambda = 780e-9;
	double k = 2.0 * M_PIl / lambda;
	double m = 86.9*1.66e-27;
	double kB = 1.38e-23; // Boltzmann constant
	double Gamma = 2*M_PIl * 6.07e6;
	double hbar = 1.05457148e-34;
	double epsilon0 = 8.85418782e-12;
	double ddip = sqrt((hbar * 2*M_PIl) * epsilon0/(k*k*k*Gamma));
	// ***************************

	// *** Simulation parameters ***
	double Tf = 10e-6; // Final time in second
	double dt = 5e-9; // Time step in second
	double T = 200e-6; // Temperature in K
	double Nreal = 1; //

	// #TODO: here an input file for simulation parameters
	int Ntimeg;
	gsl_matrix *timeg = setLinSpace(0.0, dt, Tf, &Ntimeg); // Time grid for g2
	printf("length(timeg) = %d\n", Ntimeg);

	double b0 = 0.5;
	double delta = 0; //(-16:2:16);
	int	N=300;

	double sigma = sqrt(3*N/b0); //./(1+4*delta^2)
	double R = sigma / k;
	double Ro = 50.0 * R; // distance of observation

	int Nphi = 32;
	double dphi = 2*M_PIl / Nphi;

	gsl_matrix *Phi = setLinSpace(0.0, dphi, 2*M_PIl - dphi, &Nphi);
	printf("length(Phi) = %d\n", Nphi);

	//ThObs=[(0:5/sigma/100:10/sigma) (pi-10/sigma:5/sigma/100:pi)];
	double Dthe = 5.0/(sigma*100);
	int Nthe = 2*((int)round((10.0/sigma) / Dthe) + 1);
	gsl_matrix* ThObs = gsl_matrix_alloc(Nthe, 1);
	printf("length(ThObs) = %d\n", Nthe);

	gsl_matrix_set(ThObs, 0, 0, 0.0);
	gsl_matrix_set(ThObs, Nthe/2, 0, M_PIl - 10.0/sigma);
	for(int i = 1; i < Nthe/2; i++){
		gsl_matrix_set(ThObs, i, 0, i*Dthe);
		gsl_matrix_set(ThObs, i + (Nthe/2), 0, gsl_matrix_get(ThObs, Nthe/2, 0) + i*Dthe);
	}

	gsl_matrix *Xobs = setGrid(0, Ro, ThObs, Nthe, Phi, Nphi);
	gsl_matrix *Yobs = setGrid(1, Ro, ThObs, Nthe, Phi, Nphi);
	gsl_matrix *Zobs = setGrid(2, Ro, ThObs, Nthe, Phi, Nphi);

	FILE *grid3d;
	grid3d = fopen("./grid1.dat", "w");
	for(int i = 0; i < Nthe/2; i++)
		for(int j = 0; j < Nphi; j++)
			fprintf(grid3d, "%.12lf\t%.12lf\t%.12lf\n",
														gsl_matrix_get(Xobs, i, j),
														gsl_matrix_get(Yobs, i, j),
														gsl_matrix_get(Zobs, i, j));
	fclose(grid3d);

	grid3d = fopen("./grid2.dat", "w");

	for(int i = Nthe/2; i < Nthe; i++)
		for(int j = 0; j < Nphi; j++)
		fprintf(grid3d, "%.12lf\t%.12lf\t%.12lf\n",
													gsl_matrix_get(Xobs, i, j),
													gsl_matrix_get(Yobs, i, j),
													gsl_matrix_get(Zobs, i, j));
	fclose(grid3d);

	// Dynamics
	double h=1e-5; // Numerical tolerance
	//options = odeset('AbsTol',h,'RelTol',h);

	double E0 = 1.0;
	// double complex out = Elaser(1.0, 2.0, 3.0, E0, k);
	// printf("Complex: %.5lf + (%.5lfi)\n", creal(out), cimag(out));

	//Trap
	double omTz = sqrt(kB*T/m)/R; // \omega_trap
	double omTx = omTz*M_PIl/3.0; // Factor pi/3 to break the system periodicity
	double omTy = omTz*3.0/M_PIl;

	/*double g2x=gpuArray(zeros(length(timeg),length(ThObs)));
	g12x=g2x;
	g2y=g2x;
	g12y=g2x;
	g2z=g2x;
	g12z=g2x;*/

	gsl_matrix* g2x = gsl_matrix_calloc(Ntimeg, Nthe);
	gsl_matrix* g12x = gsl_matrix_calloc(Ntimeg, Nthe);
	gsl_matrix* g2y = gsl_matrix_calloc(Ntimeg, Nthe);
	gsl_matrix* g12y = gsl_matrix_calloc(Ntimeg, Nthe);
	gsl_matrix* g2z = gsl_matrix_calloc(Ntimeg, Nthe);
	gsl_matrix* g12z = gsl_matrix_calloc(Ntimeg, Nthe);

	for(int nr = 1; nr <= Nreal; nr++){

		//double* normrnd(double mu, double sigma,int rows, int cols)
		//Cloud. The sqrt(2) is justified by the fact that the position of the
		// particles will be X*cos(omega*t+PhiX)
		gsl_matrix *X = normrnd(0.0,R,N,1);
		gsl_matrix *Y = normrnd(0.0,R,N,1);
		gsl_matrix *Z = normrnd(0.0,R,N,1);
		// initial phase calculated with a gaussian distribution for the velocity
		gsl_matrix *VX = normrnd(0.0,sqrt(kB*T/m),N,1);
		gsl_matrix *VY = normrnd(0.0,sqrt(kB*T/m),N,1);
		gsl_matrix *VZ = normrnd(0.0,sqrt(kB*T/m),N,1);

		// FILE *normDist;
		// double *pdf = (double*) malloc(N*sizeof(double));
		// normDist = fopen("./norm.dat", "w");
		//
		// for(int i = 0; i < N; i++){
		// 	pdf[i] = gsl_ran_gaussian_pdf(X[i], R);
		// 	fprintf(normDist, "%.12lf \t %.12lf \n", X[i], pdf[i]);
		// }
		// fclose(normDist);
		// free(pdf);

		gsl_matrix *PhiX = atanVet(VX, omTx, X, N);
		gsl_matrix *PhiY = atanVet(VY, omTy, Y, N);
		gsl_matrix *PhiZ = atanVet(VZ, omTz, Z, N);

		gsl_matrix *X0 = setD0(X, PhiX, N);
		gsl_matrix *Y0 = setD0(Y, PhiY, N);
		gsl_matrix *Z0 = setD0(Z, PhiZ, N);

		double Delta=delta*Gamma;

		//Initial state, taken to be the steady-state of the system for the static system with X(t=0), etc
		// M0=kernel(X,Y,Z,-omTz*Z.*tan(PhiZ));
		//M03=kernel3(X0,Y0,Z0);

		gsl_matrix_free(X0);
		gsl_matrix_free(Y0);
		gsl_matrix_free(Z0);

		gsl_matrix_free(PhiX);
		gsl_matrix_free(PhiY);
		gsl_matrix_free(PhiZ);

		gsl_matrix_free(VX);
		gsl_matrix_free(VY);
		gsl_matrix_free(VZ);

		gsl_matrix_free(X);
		gsl_matrix_free(Y);
		gsl_matrix_free(Z);
	}

	gsl_matrix_free(g2x);
	gsl_matrix_free(g12x);
	gsl_matrix_free(g2y);
	gsl_matrix_free(g12y);
	gsl_matrix_free(g2z);
	gsl_matrix_free(g12z);

	gsl_matrix_free(Xobs);
	gsl_matrix_free(Yobs);
	gsl_matrix_free(Zobs);

	gsl_matrix_free(ThObs);
	gsl_matrix_free(Phi);
	gsl_matrix_free(timeg);
	return EXIT_SUCCESS;
}
