#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <complex.h>
#include "potential.h"

#include "tools.h"

void
elementwise_addition(
                     double *res,
                     double *v1,
                     double *v2,
                     unsigned int len
                    )
{
	for(int ix = 0; ix<len;++ix)
		res[ix] = v1[ix]+v2[ix];
}

void
elementwise_multiplication(
                           double *res,
                           double *v1,
                           double *v2,
                           unsigned int len
                          )
{
	for(int ix = 0; ix<len;++ix)
		res[ix] = v1[ix]*v2[ix];

}

void
addition_with_constant(
                       double *res,
                       double *v,
                       double constant,
                       unsigned int len)
{
	for(int ix = 0; ix<len;++ix)
		res[ix] = v[ix]+constant;

}

void
multiplication_with_constant(
                             double *res,
                             double *v,
                             double constant,
                             unsigned int len)
{
	for(int ix = 0; ix<len;++ix)
		res[ix] = v[ix]*constant;

}

double
dot_product(
            double *v1,
            double *v2,
            unsigned int len
           )
{
	double sum = 0.;
	for(int ix = 0; ix<len; ++ix)
		sum += v1[ix]*v2[ix];
    return sum;
}

double **
create_2D_array(
                unsigned int row_size,
                unsigned int column_size
               )
{
	double *array_entries = malloc(sizeof(double)*row_size*column_size);
    double **array = malloc(sizeof(double*)*row_size);
	for(int ix = 0; ix<row_size; ++ix)
		array[ix]=&array_entries[ix*column_size];

    return array;
}

void
destroy_2D_array(
                 double **array 
                )

{
/*
	unsigned int n = sizeof(array)/sizeof(double);
	for(int ix = n-1; ix>=0; --ix)
		free(array[ix]);
*/
	free(array[0]);
	free(array);
}

void
matrix_vector_multiplication(
                             double *result,
                             double **A,
                             double *b,
                             unsigned int n,
                             unsigned int m
                            )
{
	for(int ix = 0; ix<n; ++ix)
		result[ix] = dot_product(A[ix],b,m);
}

void
matrix_matrix_multiplication(
                             double **result,
                             double **A,
                             double **B,
                             unsigned int n,
                             unsigned int m,
                             unsigned int k
                            )
{
	double sum;
	for(int ix = 0; ix < n; ++ix)
		for(int jx = 0; jx < k; ++jx){
			sum = 0.;
			for(int kx = 0; kx < m; ++kx)
				sum+=A[ix][kx]*B[kx][jx];
			result[ix][jx] = sum;
		}
}

double
vector_norm(
            double *v1,
            unsigned int len
           )
{
	double result = dot_product(v1,v1,len);
	result = sqrt(result);
    return result;
}


void
normalize_vector(
                 double *v1,
                 unsigned int len
                )
{
	double norm = 1/vector_norm(v1,len);
	multiplication_with_constant(v1,v1,norm,len);
}

double
average(
        double *v1,
        unsigned int len
       )
{
	double sum = 0.;
	for(int ix = 0; ix<len; ++ix)
		sum += v1[ix];
	sum/=len;
    return sum;
}


double
standard_deviation(
                       double *v1,
                       unsigned int len
                  )
{
	double avg = average(v1,len);
	addition_with_constant(v1,v1,-avg,len);
	double dev = dot_product(v1,v1,len);
	dev/=len;
    return sqrt(dev);
}

double
distance_between_vectors(
                         double *v1,
                         double *v2,
                         unsigned int len
                        )
{
	multiplication_with_constant(v2,v2,-1,len);
    elementwise_addition(v1,v1,v2,len);
	double distance = vector_norm(v1,len);
	
	return distance;
}

void
cumulative_integration(
                       double *res,
                       double *v,
                       double dx,
                       unsigned int v_len
                      )
{
	res[0] = 0;
	for(int ix = 1; ix<v_len;++ix)
		res[ix] = res[ix-1]+(v[ix]+v[ix-1])*dx/2;
}


void
write_xyz(
          FILE *fp,
          char *symbol,
          double **positions,
          double **velocities,
          double alat,
          int natoms)
{
    fprintf(fp, "%i\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" ", natoms, alat, alat, alat);
    fprintf(fp, "Properties=species:S:1:pos:R:3:vel:R:3 pbc=\"T T T\"\n");
    for(int i = 0; i < natoms; ++i){
        fprintf(fp,"%s %f %f %f %f %f %f\n",symbol,positions[i][0],positions[i][1],positions[i][2],
				velocities[i][0],velocities[i][1],velocities[i][2]);
    }
}

void fft_freq(
          double *res,
              int n,
              double timestep)
{
	int num_pos = (n-1)/2 + 1;
	double stepsize = 2*M_PI/(n*timestep);
	for(int ix = 0; ix<num_pos; ++ix)
		res[ix] = ix*stepsize;
	for(int ix = num_pos; ix<n; ++ix)
		res[ix] = (ix-n)*stepsize;
}

/* Freely given functions */
void
skip_line(FILE *fp)
{
    int c;
    while (c = fgetc(fp), c != '\n' && c != EOF);
}

void
read_xyz(
         FILE *fp,
         char *symbol,
         double **positions,
         double **velocities,
         double *alat)
{
    int natoms;
    if(fscanf(fp, "%i\nLattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 %lf\" ", &natoms, alat, alat, alat) == 0){
        perror("Error");
    }
    skip_line(fp);
    for(int i = 0; i < natoms; ++i){
        fscanf(fp, "%s %lf %lf %lf ",
                symbol, &positions[i][0], &positions[i][1], &positions[i][2]);
        fscanf(fp, "%lf %lf %lf\n",
                &velocities[i][0], &velocities[i][1], &velocities[i][2]);
    }
}

void powerspectrum(
           double *res,
           double *signal,
           int n,
                   double timestep)
{
    /* Declaration of variables */
    double *complex_coefficient = malloc(sizeof(double) * 2*n); // array for the complex fft data
    double *data_cp = malloc(sizeof(double) * n);

    /*make copy of data to avoid messing with data in the transform*/
    for (int i = 0; i < n; i++) {
    data_cp[i] = signal[i];
    }

    /* Declare wavetable and workspace for fft */
    gsl_fft_real_wavetable *real;
    gsl_fft_real_workspace *work;

    /* Allocate space for wavetable and workspace for fft */
    work = gsl_fft_real_workspace_alloc(n);
    real = gsl_fft_real_wavetable_alloc(n);

    /* Do the fft*/
    gsl_fft_real_transform(data_cp, 1, n, real, work);

    /* Unpack the output into array with alternating real and imaginary part */
    gsl_fft_halfcomplex_unpack(data_cp, complex_coefficient,1,n);

    /*fill the output powspec_data with the powerspectrum */
    for (int i = 0; i < n; i++) {
    res[i] = (complex_coefficient[2*i]*complex_coefficient[2*i]+complex_coefficient[2*i+1]*complex_coefficient[2*i+1]);
    res[i] *= timestep / n;
    }

    /* Free memory of wavetable and workspace */
    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
    free(complex_coefficient);
    free(data_cp);
}

void calculate_acceleration(
		double **accelerations,
		double **forces,
		double mass,
		int N){	
	for(int ix = 0; ix<N; ++ix)
		for(int jx = 0; jx<3;++jx)
			accelerations[ix][jx] = forces[ix][jx]/mass;
}

void velocity_verlet_one_step(
		double **accelerations, 
		double **positions, 
		double **velocities,
		double **forces,
        double mass, 
		double timestep,
		int n_atoms,
		int N,
		double lattice_const)
{

    // Write to accelerations, positions and velocities vectors
	for(int ix = 0; ix<n_atoms;++ix)
		for(int jx = 0; jx<3; ++jx){
		velocities[ix][jx] += timestep*accelerations[ix][jx]/2;
		positions[ix][jx] += velocities[ix][jx]*timestep;
	}
	get_forces_AL(forces,positions,N*lattice_const,n_atoms);
	calculate_acceleration(accelerations, forces, mass,n_atoms);

	for(int ix = 0; ix<n_atoms;++ix)
		for(int jx = 0; jx < 3; ++jx)
			velocities[ix][jx] += timestep*accelerations[ix][jx]/2;
}

double calculate_kinetic_energy(
		double **velocities,
		double mass,
		int n_atoms
		){
	double kinetic_energy = 0.;
	for(int ix = 0; ix<n_atoms; ++ix)
		for(int jx = 0; jx < 3; ++jx)
			kinetic_energy += velocities[ix][jx]*velocities[ix][jx]*mass/2.0;
	return kinetic_energy;
}
/*
void velocity_verlet_one_step_temp(
		double **accelerations, 
		double **positions, 
		double **velocities,
		double **forces,
		double temperature,
        double mass, 
		double timestep,
		int n_atoms,
		int N,
		double lattice_const,
		double tau_t,
		double T_eq)
{

    // Write to accelerations, positions and velocities vectors
	for(int ix = 0; ix<n_atoms;++ix)
		for(int jx = 0; jx<3; ++jx){
			temperature = update_velocity_half(velocities, temperature,timestep,tau_t,T_eq,n_atoms);
			positions[ix][jx] += velocities[ix][jx]*timestep;
	}
	get_forces_AL(forces,positions,N*lattice_const,n_atoms);
	calculate_acceleration(accelerations, forces, mass,n_atoms);

	for(int ix = 0; ix<n_atoms;++ix)
		for(int jx = 0; jx < 3; ++jx)
			velocities[ix][jx] += timestep*accelerations[ix][jx]/2;
}
*/
/*
double update_velocity_half(
		double **velocities,
		double temperature,
		double dt,
		double tau_t,
		double T_eq,
		int n_atoms
		){
	double alpha_T = 1. + dt*(T_eq-temperature)/(tau_t*temperature);
	for(int ix = 0; ix < n_atoms; ++ix)
		for(int jx = 0; jx < 3; ++jx)
			velocities[ix][jx] = sqrt(alpha_T)*velocities[ix][jx];
	temperature *= (1+alpha_T)/2.;
	return temperature;
}*/

void scale_velocity(
		double **velocities,
		double temperature,
		double T_eq,
		double tau_t,
		double dt,
		int n_atoms
		){
	double alpha_t = 1 + 2*dt*(T_eq-temperature)/(tau_t*temperature);
	for(int ix = 0; ix < n_atoms; ++ix)
		for(int jx = 0; jx < 3; ++jx)
			velocities[ix][jx] *= sqrt(alpha_t);
}

void scale_pressure(
		double **positions,
		double dt,
		double tau_P,
		double P_eq,
		double n_atoms,
		double pressure,
		double kappa_T,
		double *lattice_const
		){
	double alpha_P = 1-kappa_T*dt*(P_eq-pressure)/tau_P;
	for(int ix = 0; ix<n_atoms;++ix)
		for(int jx = 0; jx<3; ++jx)
			positions[ix][jx] *= pow(alpha_P,0.3333333);
	*lattice_const *= pow(alpha_P,0.3333333);
}

void
write_csv(
		FILE *file,
		double **array,
		int N,
		int M
		)
{
	for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            fprintf(file, "%f", array[i][j]);
            if (j < M - 1) {
                fprintf(file, ",");
            }
        }
        fprintf(file, "\n");
    }
	fclose(file);
}

double calculate_temperature(
		double **velocities,
		double mass,
		int n_atoms
		){
	double e_k = calculate_kinetic_energy(velocities,mass,n_atoms);
	double k_b = 8.617333*1e-5;
	return (2.*e_k)/(3.*n_atoms*k_b);
}

double calculate_pressure(
		double temperature,
		double virial,
		double total_length,
		int n_atoms
		){
	double volume = total_length*total_length*total_length;
	double k_b = 8.617333*1e-5;
	return (n_atoms*k_b*temperature + virial)/volume;
}
