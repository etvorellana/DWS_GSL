/*
 * dws.h
 *
 *  Created on: 16 de out de 2017
 *      Author: evalero
 */

#include <stdio.h>
#include <stdlib.h>

#define __USE_GNU
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

#include <complex.h>

#ifndef DWS_H_
#define DWS_H_

gsl_matrix* setLinSpace(double T_0, double dt, double T_F, int *ls_size);
gsl_matrix* setGrid(int dim, double Ro, gsl_matrix* ThObs, int Nthe,
                                        gsl_matrix* Phi, int Nphi);
gsl_matrix* normrnd(double mu, double sigma,int rows, int cols);
gsl_matrix* atanVet(gsl_matrix* Vd, double omTd, gsl_matrix* d, int N);
gsl_matrix* setD0(gsl_matrix* D, gsl_matrix* PhiD, int N);
gsl_matrix* ones(int rows, int cols);
double complex Elaser(double x, double y, double z, double E0, double k);

#endif /* DWS_H_ */
