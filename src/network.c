////////////////////////////////////////////////////////////////////////////////////////
// Kalle Kossio et al. 2018, Phys. Rev. Lett. 121, 058301,
// https://doi.org/10.1103/PhysRevLett.121.058301
////////////////////////////////////////////////////////////////////////////////////////

#include "network.h"

void evolve_one_step(double* time_sle, int* spiked_nrn, int N, 
                     neuron net[N], const gsl_rng* r, int growing){
  /* The core simulation function. Evolves the given network one step until the next spike
   * (or "nonspike" if T_REF!=0). Returns time to the next spike (or "nonspike" if T_REF!=0)
   * and a neuron that spiked. Bullean variable growing tells if there should be change
   * in radia of neurons. Growing = 0 is equivalent to K = 0. */

  /* Loop counter */
  unsigned int i;

  /* Picking next neuron to spike */
  double time_sls = spike_time(&net[0], r);
  int spiked_neuron = 0;

  for(i=1; i<N; i++){ 
    double t = spike_time(&net[i], r);  
    if( t < time_sls ) {spiked_neuron = i; time_sls = t;}
  } 

  if( net[spiked_neuron].t_ref_left >= time_sls )
    spiked_neuron = none;

  /* Refractory */
  for(i=0; i<N; i++)
    if(net[i].t_ref_left > 0.0)
      net[i].t_ref_left -= time_sls;
  if(spiked_neuron != none)
    net[spiked_neuron].t_ref_left = net[spiked_neuron].t_ref;


  /* Adjusting radia and firing rates after spike */
  for(i=0; i<N; i++)
  {
    net[i].f = net[i].f_0 + (net[i].f - net[i].f_0) * exp(-time_sls/net[i].tau);
    net[i].r += net[i].k * time_sls;
    if(spiked_neuron != none && i != spiked_neuron )
      net[i].f += net[i].g * overlap(&net[spiked_neuron], &net[i]);
  }

  if(spiked_neuron != none && growing){
    net[spiked_neuron].r -= net[spiked_neuron].k / net[spiked_neuron].f_sat;
    if (net[spiked_neuron].r < 0.0)
      net[spiked_neuron].r = 0.0;
  }

  *time_sle = time_sls;
  *spiked_nrn = spiked_neuron;
}

void sample_radii(const int sample_freq, FILE* fp, unsigned int* counter,
                          const int spiked_neuron, double* time, 
                          const int N, neuron net[N]){
  if(spiked_neuron != none)
    (*counter)++;
  if(*counter == sample_freq){
    *counter = 0;
    fwrite(time, sizeof(double), 1, fp);
    int i;
    for (i = 0; i<N; i++){
      fwrite(&net[i].r, sizeof(double), 1, fp);
    } 
  }
}

void sample_total_overlap(const int sample_freq, FILE* fp, unsigned int* counter,
                          const int spiked_neuron, double* time, 
                          const int N, neuron net[N]){
  if(spiked_neuron != none)
    (*counter)++;
  if(*counter == sample_freq){
    int i,j;
    double tot_overlap[N];

    for(i=0; i<N; i++){
      tot_overlap[i] = 0.0;

      for(j=0; j<N; j++)
        if(j != i)
          tot_overlap[i] += overlap(&net[i], &net[j]);
    }

    fwrite(time, sizeof(double), 1, fp);
    fwrite(tot_overlap, sizeof(double), N, fp);
    *counter = 0;
  }
}
 
double overlap(neuron* n1, neuron* n2){
  /* Calculating area of overlap between two neurons */
  double d = gsl_hypot(n1->x - n2->x, n1->y - n2->y);
  double r1 = n1->r;
  double r2 = n2->r;

  if(d >= r1+r2)
    return 0.0;
  else if(r1+d <= r2)
    return r1*r1*M_PI;
  else if(r2+d <= r1)
    return r2*r2*M_PI;
  else {
    double theta1 = 2.0*acos((r1*r1+d*d-r2*r2)/(2.0*r1*d));
    double theta2 = 2.0*acos((r2*r2+d*d-r1*r1)/(2.0*r2*d));

    double tri1_area = 0.5*r1*r1*sin(0.5*theta1);
    double tri2_area = 0.5*r2*r2*sin(0.5*theta2);
    
    double A1 = 0.5*theta1*r1*r1-tri1_area;
    double A2 = 0.5*theta2*r2*r2-tri2_area;

    return A1+A2;}
}

double spike_time(neuron* n, const gsl_rng* r){
  /* Calculate next interspike interval due to f_0 rate alone, and if firing 
   * rate f is equal to f_0 just return the result. Otherwise 
   * calculate next interspike interval due to f-f_0 rate alone
   * and return smaller of the two intervals. */
  
  double t1;
  t1 = - log( gsl_rng_uniform( r ) ) / n->f_0;

  if(0 == gsl_fcmp(n->f - n->f_0, 0.0, 1e-10))
    return t1;
  
  double k = 1.0 + log(gsl_rng_uniform(r))/(n->tau * (n->f - n->f_0));
  if(k <= 0.0)
    return t1;
  
  double t2 = - n->tau * log(k); 
  if(t1 <= t2)
    return t1;
  
  return t2;
}

void save_state(const char file_name[], const int N, const neuron net[N], 
                const double time, gsl_rng* r){

  /* Function for saving state of the simulation in binary format */
  FILE* fp;
  fp = fopen(file_name, "w");
  
  fwrite(&N, sizeof(int), 1, fp);
  fwrite(&time, sizeof(double), 1, fp);

  int i;
  for (i = 0; i<N; i++){
    fwrite( &net[i].r, sizeof(double), 1, fp );
    fwrite( &net[i].x, sizeof(double), 1, fp );
    fwrite( &net[i].y, sizeof(double), 1, fp );
    fwrite( &net[i].f, sizeof(double), 1, fp );
    fwrite( &net[i].t_ref_left, sizeof(double), 1, fp );
    fwrite( &net[i].t_ref, sizeof(double), 1, fp );
    fwrite( &net[i].f_0, sizeof(double), 1, fp );
    fwrite( &net[i].tau, sizeof(double), 1, fp );
    fwrite( &net[i].k, sizeof(double), 1, fp );
    fwrite( &net[i].g, sizeof(double), 1, fp );
    fwrite( &net[i].f_sat, sizeof(double), 1, fp );
  }
  
  gsl_rng_fwrite(fp, r);
 
  fclose(fp);
}

void load_state(const char file_name[], int* N, neuron** net, 
                double* time, gsl_rng* r){

  FILE* fp;
  fp = fopen(file_name, "r");

  size_t totl = 0; /*Total size of data loaded*/
  
  totl += fread(N, sizeof(int), 1, fp);
  totl += fread(time, sizeof(double), 1, fp);

  *net = (neuron *) malloc( (*N) * sizeof(neuron) );

  int i;
  for (i = 0; i<*N; i++){
    totl += fread( &(*net)[i].r, sizeof(double), 1, fp );
    totl += fread( &(*net)[i].x, sizeof(double), 1, fp );
    totl += fread( &(*net)[i].y, sizeof(double), 1, fp );
    totl += fread( &(*net)[i].f, sizeof(double), 1, fp );
    totl += fread( &(*net)[i].t_ref_left, sizeof(double), 1, fp );
    totl += fread( &(*net)[i].t_ref, sizeof(double), 1, fp );
    totl += fread( &(*net)[i].f_0, sizeof(double), 1, fp );
    totl += fread( &(*net)[i].tau, sizeof(double), 1, fp );
    totl += fread( &(*net)[i].k, sizeof(double), 1, fp );
    totl += fread( &(*net)[i].g, sizeof(double), 1, fp );
    totl += fread( &(*net)[i].f_sat, sizeof(double), 1, fp );
  }

  if(totl != 2+11*(*N))
    fprintf(stderr, "ERROR: There was some error loading network state.\n");
  
  gsl_rng_fread(fp, r);

  fclose(fp);
}

void show_progress(double now, double end){
  #define WIDTH 20
  double ratio = now/end;
  printf("%3d%% [", (int)(ratio*100.0));
  
  int i;
  for(i=0; i<(int)(ratio*WIDTH);i++)
    printf("=");
  
  for(i=(int)(ratio*WIDTH); i<WIDTH; i++)
    printf(" ");

  printf("]\n\033[F\033[J"); 
}

