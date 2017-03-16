#include "header.h"

int main(int argc, char** argv)
{

/*************** Define Variables *******************/

	// Simulation Parameters
	int n_r = 129;
	int n_mu = 33; // mu = cos(theta) where theta is angle from axis of rotation
	int n_phi = 33;
	int L_max = 16; // number of legendre functions to keep
	int max_iter = 10; // max number of evolution steps for planet

	// Initial axes for Ellipsoid (a/b will be fixed in sim)
	double a = 9.6e5;
	double b = 7.7e5;
	double c = 8.95e5;
	double r_max = 16./15*a; // maximum r-value > R_e to allow for max r off equator
	
	// Minerology
	double rho_0 = 2600;  // kg/m^3
	double K0 = 65*1e9;  //K0 Olivine in mks
	double K0_prime = 6.54;  //K0_prime Olivine dimless	

	// Useful Information
	double G = 6.67e-11;
	double real_mass = 4*M_PI/3.*a*b*c*rho_0; // Benchmark for ellipsoid
	double real_W = 3*G*pow(real_mass,2)*2.714/1e7; // Benchmark for ellipsoid w/ a=9600km,b=7700km,c=4950km
//	double real_W = 3*G*pow(real_mass,2)/(5.*a);	// Benchmark for sphere


/**************** Initialize Planet *********************/

	Planet Haumea = Planet(n_r,n_mu,n_phi,L_max);  // create planet named haumea on the given grid
	Haumea.init_density(a,b,c,rho_0,r_max);  // initialize size and density
	Haumea.get_potential(L_max);  // gravitational potential at each cell
	Haumea.init_Q1();
	Haumea.init_Q2();
	Haumea.init_S1();
	Haumea.init_S2();
	Haumea.get_enthalpy(K0,K0_prime,rho_0);  // enthalpy

	// Calculate benchmark quantities using methods defined in "Planet" class
	Haumea.Mass = Haumea.get_mass();
	Haumea.W = Haumea.get_W();


/*************** Evolve Planet *********************/

	// Open File to Hold Simulation Parameters
	std::ofstream params_file("Params.txt");
	params_file << "n_r\t" << "n_mu\t" << "n_phi\t" << "a\t" << "b\t" << "c\t" 
	<< "max_iter" << std::endl;
	params_file << n_r << "\t" << n_mu << "\t" << n_phi << "\t" << a << "\t" << b << "\t" << c << "\t" 
	<< max_iter << std::endl;	
	params_file.close();

	// Open Files to Hold Data
	std::ofstream density_file("Density.txt");
	std::ofstream data_file("Data.txt");
	density_file << "iter\t" << "phi\t" << "mu\t" << "r\t" << "Density" << "\t" << "Potential" << std::endl;
	data_file << "iter\t" << "Period[hours]\t" << "Mass[g]\t" << std::endl;

	//Evolve Planet
	for(int i=0;i<max_iter;i++){
		Haumea.iter = i;
		density_file << Haumea; // write step to file
		std::cout << "Iter:\t" << i << "\tPeriod\t" << 2*M_PI/(sqrt(Haumea.omega_sq)*3600) << 
		"\tMass\t" << Haumea.Mass << std::endl;
		data_file << i << "\t" << 2*M_PI/(sqrt(Haumea.omega_sq)*3600) << "\t" << Haumea.Mass << std::endl;
		Haumea.evolve(K0,K0_prime,rho_0,L_max); // updata
	}

	density_file.close();
	data_file.close();
	std::cout << "Done" << std::endl;

}