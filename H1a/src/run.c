#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include <gsl/gsl_rng.h>
#include <time.h>

int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code
	int N = 4;
	int n_atoms = 4*N*N*N;
	double lattice_const[] = {4.0, 4.01, 4.02, 4.03, 4.04, 4.05, 4.06, 4.07, 4.08}; //Å
	int lattice_const_len = sizeof(lattice_const)/sizeof(lattice_const[0]);
	double pot_energy[lattice_const_len];
	double **positions = create_2D_array(n_atoms,3);

	
	for(int ix = 0; ix < lattice_const_len; ++ix){
		//initialize positions
		init_fcc(positions,N,lattice_const[ix]);

		//calculate potential energy
		pot_energy[ix] = get_energy_AL(positions, N*lattice_const[ix],n_atoms)/(N*N*N);
	}
	FILE *fp= fopen("ex1_pot_energy.csv", "w");
	for(int ix = 0; ix<lattice_const_len; ++ix){
		fprintf(fp,"%f",pot_energy[ix]);
		if(ix != lattice_const_len-1)
			fprintf(fp,",");
	}
	fclose(fp);
	destroy_2D_array(positions);
    return 0;
}

int
run2(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code
	int N = 4;
	int n_atoms = 4*N*N*N;
	double lattice_const = 4.03;
	double **positions = create_2D_array(n_atoms,3);
	double **accelerations = create_2D_array(n_atoms,3);
	double **velocities = create_2D_array(n_atoms,3);
	double **forces = create_2D_array(n_atoms,3);
	double displacement_const = 0.065;
	double mass_AL = 2.7963*1e-3;
	double total_time = 1;
	double dt = 1e-3;
	int timesteps = total_time/dt;
	double *total_energies = malloc(sizeof(double)*timesteps);
	double *kinetic_energies = malloc(sizeof(double)*timesteps);
	double *potential_energies = malloc(sizeof(double)*timesteps);

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, time(NULL)); 



	//initial values
	init_fcc(positions,N,lattice_const);
	for(int ix = 0; ix<n_atoms; ++ix)
		for(int jx = 0; jx < 3; ++jx){
			positions[ix][jx] += lattice_const*displacement_const*(-1.0+2*gsl_rng_uniform(r));
			velocities[ix][jx] = 0;
		}
	get_forces_AL(forces,positions,N*lattice_const,n_atoms);
	calculate_acceleration(accelerations,forces,mass_AL,n_atoms);
	kinetic_energies[0] = calculate_kinetic_energy(velocities,mass_AL,n_atoms);
	potential_energies[0] = get_energy_AL(positions, N*lattice_const,n_atoms);
	total_energies[0] = kinetic_energies[0]+potential_energies[0];

	for(int ix = 1; ix < timesteps; ++ix){
		velocity_verlet_one_step(accelerations,positions,velocities,forces,mass_AL,dt,n_atoms,N,lattice_const);
		kinetic_energies[ix] = calculate_kinetic_energy(velocities, mass_AL, n_atoms);
		potential_energies[ix] = get_energy_AL(positions,N*lattice_const,n_atoms);
		total_energies[ix] = potential_energies[ix]+kinetic_energies[ix];
	}
	FILE *fp= fopen("ex2_pot_energy.csv", "w");
	for(int ix = 0; ix<timesteps; ++ix){
		fprintf(fp,"%f",potential_energies[ix]);
		if(ix != timesteps-1)
			fprintf(fp,",");
	}
	fclose(fp);

	fp= fopen("ex2_kin_energy.csv", "w");
	for(int ix = 0; ix<timesteps; ++ix){
		fprintf(fp,"%f",kinetic_energies[ix]);
		if(ix != timesteps-1)
			fprintf(fp,",");
	}
	fclose(fp);

	fp = fopen("ex2_tot_energy.csv","w");
	for(int ix = 0; ix<timesteps; ++ix){
		fprintf(fp,"%f",total_energies[ix]);
		if(ix != timesteps-1)
			fprintf(fp,",");
	}
	fclose(fp);



	free(kinetic_energies);
	free(potential_energies);
	free(total_energies);
	destroy_2D_array(forces);
	destroy_2D_array(positions);
	destroy_2D_array(accelerations);
	destroy_2D_array(velocities);
	gsl_rng_free(r);
    return 0;
}

