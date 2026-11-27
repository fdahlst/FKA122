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
	double dt = 1e-1;
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
int
run3(
    int argc,
    char *argv[]
   )
{
	int N = 4;
	int n_atoms = 4*N*N*N;
	double lattice_const = 4.03;
	double **positions = create_2D_array(n_atoms,3);
	double **accelerations = create_2D_array(n_atoms,3);
	double **velocities = create_2D_array(n_atoms,3);
	double **forces = create_2D_array(n_atoms,3);
	double mass_AL = 2.7963*1e-3;
	double total_time = 20;
	double dt = 1e-3;
	double T_eq = 500.+ 273.15; //500*c
	double P_eq = 6.2415 * 1e-7; //1 bar
	double tau_T = 1.;
	double tau_P = 1.;
	double displacement_const = 0.065;
	int timesteps = total_time/dt;
	double *temperatures = malloc(sizeof(double)*timesteps);
	double *pressures = malloc(sizeof(double)*timesteps);
	double kappa_T = 2.219;
	double virial;
	
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
	temperatures[0] = 1.;
	pressures[0] = P_eq*temperatures[0]/T_eq;
	virial = get_virial_AL(positions,N*lattice_const,n_atoms);



	FILE *fp = fopen("data/task3/positions.csv","w");
	fprintf(fp, "%g, %g, %g, %g, %g, %g, %g, %g, %g \n",
			positions[0][0], positions[0][1],positions[0][2],
			positions[100][0], positions[100][1],positions[100][2],
			positions[n_atoms-1][0], positions[n_atoms-1][1],positions[n_atoms-1][2]);
	
	//run simulation
	for(int ix = 1; ix < timesteps; ++ix){
		velocity_verlet_one_step(accelerations,positions,velocities,forces,mass_AL,dt,n_atoms,N,lattice_const);
		temperatures[ix] = calculate_temperature(velocities,mass_AL,N*lattice_const);
		scale_velocity(velocities, temperatures[ix],T_eq,tau_T,dt,n_atoms);
		

		virial = get_virial_AL(positions, N*lattice_const,n_atoms);
		pressures[ix] = calculate_pressure(temperatures[ix],virial,N*lattice_const,n_atoms);
		scale_pressure(positions,dt,tau_P,P_eq,n_atoms,pressures[ix],kappa_T,&lattice_const);
		
		//write particle positions
		fprintf(fp, "%g, %g, %g, %g, %g, %g, %g, %g, %g \n",
				positions[0][0], positions[0][1],positions[0][2],
				positions[100][0], positions[100][1],positions[100][2],
				positions[n_atoms-1][0], positions[n_atoms-1][1],positions[n_atoms-1][2]);
		
	}
	fclose(fp);

	//write data

	fp = fopen("data/task3/ex3_temperatures.csv", "w");
	for(int ix = 0; ix<timesteps; ++ix){
		fprintf(fp,"%g",temperatures[ix]);
		if(ix != timesteps-1)
			fprintf(fp,",");
	}
	fclose(fp);

	fp = fopen("data/task3/ex3_pressures.csv", "w");
	for(int ix = 0; ix<timesteps; ++ix){
		fprintf(fp,"%g",pressures[ix]);
		if(ix != timesteps-1)
			fprintf(fp,",");
	}
	fclose(fp);
	
	fp = fopen("data/task3/end_positions.csv", "w");
	write_csv(fp,positions,n_atoms,3);

	printf("lattice constant: %g Å", lattice_const);

	destroy_2D_array(forces);
	destroy_2D_array(positions);
	destroy_2D_array(accelerations);
	destroy_2D_array(velocities);
	free(pressures);
	free(temperatures);
	gsl_rng_free(r);
    return 0;

}
