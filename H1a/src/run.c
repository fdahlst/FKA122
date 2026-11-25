#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"


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
	double **positions = create_2D_array(n_atoms,4);
	
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
