#ifndef FDM_BAR_H
#define FDM_BAR_H
#include <stdio.h>
#include <gsl/gsl_linalg.h>

const int model_time = 5;
const int x_nodes = 8;
const int y_nodes = 3;
const int vector_N = x_nodes * y_nodes;
const double lenght = 8;
const double hight = 3;
const double dx = lenght / x_nodes + 1;
//const double dx = 1;
const double dy = hight / y_nodes + 1;
//const double dy = 1;
const double dt = 1;
const double left_bound = 400;
const double top_bound = 0;
const double bottom_bound = 400;
const double top_right_bound = 0;
const double bottom_right_bound = 400;
const double init_condition = 10;
const double a = 2;
FILE *gp;


void initialization(gsl_vector *);

void FDM(gsl_vector *, gsl_matrix *);

void print_LAE(const gsl_vector *, const gsl_matrix *);

void print_temperature(gsl_vector *);

void get_LAE(gsl_matrix *);

void boundary_condition(gsl_vector *);

void create_plot(gsl_vector *);


#endif //FDM_BAR_H
