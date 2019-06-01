#ifndef FDM_BAR_H
#define FDM_BAR_H
#include <stdio.h>
#include <gsl/gsl_linalg.h>

const int model_time = 2;
const int x_nodes = 8;
const int y_nodes = 6;
const int vector_N = 8 * 6;
const double dt = 1;
const double left_bound = 200;
const double corner = 50;
const double init_condition = 0;
const double a = 3;
FILE *gp;


void initialization(gsl_vector *);

void FDM(gsl_vector *, gsl_matrix *);

void print_LAE(const gsl_vector *, const gsl_matrix *);

void print_temperature(gsl_vector *);

void get_LAE(gsl_matrix *, gsl_vector *);

void boundary_condition(gsl_vector *);

void create_plot(gsl_vector *);


#endif //FDM_BAR_H
