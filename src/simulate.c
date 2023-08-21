////////////////////////////////////////////////////////////////////////////////////////
// Kalle Kossio et al. 2018, Phys. Rev. Lett. 121, 058301,
// https://doi.org/10.1103/PhysRevLett.121.058301
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include "network.h"

// Update progress bar every [milliseconds, simulation time]
#define PROGRESS_BAR_UPDATE 1.0e3 

int main(int argc, char* argv[])
{
  double K = atof(argv[1]);     
  printf("Neuron growth rate [1/milliseconds]: %g\n", K);
  double F_SAT = atof(argv[2]); 
  printf("Equilibrium firing rate [1/milliseconds]: %g\n", F_SAT);
  double F_0 = atof(argv[3]);   
  printf("Spontaneous firing rate [1/milliseconds]: %g\n", F_0);
  double G = atof(argv[4]);  
  printf("Coupling proportionality constant [1/milliseconds]: %g\n", G);
  double TAU = atof(argv[5]); 
  printf("Time constant [milliseconds]: %g\n", TAU);
  double T_REF = atof(argv[6]); 
  printf("Absolute refractory period [milliseconds]: %g\n", T_REF);
  int N = atoi(argv[7]);    
  printf("Number of neurons: %i\n", N);
  double BOX_SIDE_LENGTH = atof(argv[8]); 
  printf("Maximum value for x and y coordinates of soma center [unitless]: %g\n", BOX_SIDE_LENGTH);
  double MAX_INIT_R = atof(argv[9]);
  printf("Maximum value for inital neuron radius [unitless]: %g\n", MAX_INIT_R);
  double SIM_TIME = atof(argv[10]);
  printf("Total simulation time [milliseconds]: %g\n", SIM_TIME);
  unsigned long int RNG_SEED = atoi(argv[11]); 
  printf("Random number generator seed: %lu\n", RNG_SEED);
  int SAMP = atoi(argv[12]);
  printf("How often to sample neuron radii/overlap [spikes]: %i\n", SAMP);
  char* FNAME_SPIKETIME = argv[13];
  printf("File for writing spike times: %s\n", FNAME_SPIKETIME);
  char* FNAME_NEURON = argv[14];
  printf("File for writing spiking neuron index: %s\n", FNAME_NEURON);
  char* FNAME_RADII = argv[15];
  printf("File for writing radii: %s\n", FNAME_RADII);
  char* FNAME_OVERLAP = argv[16];
  printf("File for writing total overlap: %s\n", FNAME_OVERLAP);
  char* FNAME_INITSTATE = argv[17];
  printf("File for writing initial state of the simulation: %s\n", FNAME_INITSTATE);
  char* FNAME_FINSTATE = argv[18];
  printf("File for writing final state of the simulation: %s\n", FNAME_FINSTATE);

  // Creating output files
  FILE* fp_spike_time;
  fp_spike_time = fopen(FNAME_SPIKETIME, "w");

  FILE* fp_neuron;
  fp_neuron = fopen(FNAME_NEURON, "w");
  
  FILE* fp_radii;
  FILE* fp_overlap;
  if(SAMP > 0){
    fp_radii = fopen(FNAME_RADII, "w");
    fp_overlap = fopen(FNAME_OVERLAP, "w");
  }
   
  // Initialising random number generator with these
  gsl_rng* r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set(r, RNG_SEED);
  
  // Number of neurons
  int n = N;

  // Initializing network 
  neuron net[n];
  int i;
  for(i=0; i<n; i++)
  {
    net[i].r = MAX_INIT_R * gsl_rng_uniform(r); 
    net[i].x = BOX_SIDE_LENGTH * gsl_rng_uniform(r);
    net[i].y = BOX_SIDE_LENGTH * gsl_rng_uniform(r);
    net[i].t_ref_left = 0.0; // Initially no neurons are refractory
    net[i].t_ref = T_REF;
    net[i].f_0 = F_0;
    net[i].f = net[i].f_0; // Inital firing rate is the spontaneous rate
    net[i].tau = TAU;
    net[i].k = K;
    net[i].g = G;
    net[i].f_sat = F_SAT;
  }

  double time = 0.0; // Current simulation time  
  save_state(FNAME_INITSTATE, n, net, time, r); // Saving initial state
  
  double next_progress_update = 0.0; // Next time for updating progress bar
  double time_sls = 0.0;  // Time elapesed Since Last Spike
  int spiked_neuron = none; // Last neuron that spiked
  unsigned int samp_counter_radii = 0; // Counter for sampling radii
  unsigned int samp_counter_overlap = 0; // Counter for sampling total overlap

  // Main simulation loop
  while(time <= SIM_TIME)
  {
    evolve_one_step(&time_sls, &spiked_neuron, n, net, r, 1);
    time += time_sls;
    
    // Writing data to files
    if(spiked_neuron != none){
      fwrite(&time, sizeof(double), 1, fp_spike_time);
      fwrite(&spiked_neuron, sizeof(int), 1, fp_neuron);
    }
    
    // Writing radii and total overlap area to files
    if(SAMP > 0){
      sample_radii(SAMP, fp_radii, &samp_counter_radii, spiked_neuron, &time, n, net);
      sample_total_overlap(SAMP, fp_overlap, &samp_counter_overlap, spiked_neuron, &time, n, net);
    }

    // Update progress bar
    if(time >= next_progress_update){
      next_progress_update += PROGRESS_BAR_UPDATE;
      show_progress(time, SIM_TIME);
    }
  }
  
  save_state(FNAME_FINSTATE, n, net, time, r); // Saving final state of the simulation
  
  gsl_rng_free(r); // Closing random number generator
  
  // Closing files
  fclose(fp_spike_time);
  fclose(fp_neuron);
  if(SAMP > 0){
    fclose(fp_radii);
    fclose(fp_overlap);
  }
  return 0;
}

