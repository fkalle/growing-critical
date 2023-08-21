////////////////////////////////////////////////////////////////////////////////////////
// Kalle Kossio et al. 2018, Phys. Rev. Lett. 121, 058301,
// https://doi.org/10.1103/PhysRevLett.121.058301
////////////////////////////////////////////////////////////////////////////////////////

#ifndef NETWORK_H 
#define NETWORK_H 

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <math.h>

#define none -1 /* spiked neuron = none, if no neurons spiked yet */

typedef struct neuron{
  double r;             /* radius */ 
  double x;             /* x-coordinate */ 
  double y;             /* y-coordinate */ 
  double f;             /* spike rate */
  double t_ref_left;    /* remaining time in refractory state */
  double t_ref;         /* refractory period */
  double f_0;           /* spontaneous spike rate */
  double tau;           /* decay time constant of spike rate */
  double k;             /* growth rate */
  double g;             /* coupling constant */
  double f_sat;         /* desired spike rate */
} neuron;

double overlap(neuron* n1, neuron* n2);
double spike_time(neuron* n, const gsl_rng* r);
void save_state(const char file_name[], const int N, const neuron net[N], 
               const double time, gsl_rng* r);
void load_state(const char file_name[], int* N, neuron** net, 
                double* time, gsl_rng* r);
void evolve_one_step(double* time_sle, int* spiked_nrn, int N, neuron net[N], const gsl_rng* r, int growing);
void sample_radii(const int sample_freq, FILE* fp, unsigned int* counter,
                          const int spiked_neuron, double* time, 
                          const int N, neuron net[N]);
void sample_total_overlap(const int sample_freq, FILE* fp, unsigned int* counter,
                          const int spiked_neuron, double* time, 
                          const int N, neuron net[N]);
                          
void show_progress(double now, double end);
#endif
