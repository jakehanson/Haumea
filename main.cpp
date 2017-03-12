#include "header.h"
#include <fstream>

int main(int argc, char** argv)
{

/*************** Define Variables *******************/

	// Simulation Parameters
	int n_r = 129;
	int n_mu = 33; // mu = cos(theta) where theta is angle from axis of rotation
	int n_phi = 65;
	int L_max = 8; // number of legendre functions to keep
	int max_iter = 6;

	// Initial axes for Ellipsoid (a/b will be fixed in sim)
	double a = 9.6e7;
	double b = 7.7e7;
	double c = 4.95e7;
	double r_max = 16./15*a; // maximum r-value > R_e to allow for max r off equator

	// Minerology 
	double rho_0 = 3.81;
	double K0 = 126.54*1e10;  //K0 Olivine in cgs (Ba)
	double K0_prime = 4.3;  //K0_prime Olivine dimless	

	// Useful Information
	double G = 6.67e-8;
	double real_mass = 4*M_PI/3.*a*b*c*rho_0; // Benchmark for ellipsoid
	double real_W = 3*G*pow(real_mass,2)*2.714/1e9; // Benchmark for ellipsoid w/ a=9600km,b=7700km,c=4950km
//	double real_W = 3*G*pow(real_mass,2)/(5.*a);	// Benchmark for sphere


/**************** Build Planet *********************/

	Planet Haumea = Planet(n_r,n_mu,n_phi,L_max);  // create planet named haumea on the given grid
	Haumea.init_density(a,b,c,rho_0,r_max);  // initialize size and density
	
	// Initialize all relevant functions and values
	Haumea.init_P();
	Haumea.init_D1();
	Haumea.init_D2();
	Haumea.init_D3();
	Haumea.get_potential();  // gravitational potential at each cell
	Haumea.init_Q1();
	Haumea.init_Q2();
	Haumea.init_S1();
	Haumea.init_S2();
	Haumea.get_enthalpy(K0,K0_prime,rho_0);  // enthalpy

	// Calculate benchmark quantities using methods defined in "Planet" class
	Haumea.Mass = Haumea.get_mass();
	Haumea.W = Haumea.get_W();


/*************** Evolve Planet *********************/

	// Open File to Hold Density Data
	std::ofstream data_file("Density.txt");
	data_file << "iter\t" << "phi\t" << "mu\t" << "r\t" << "Density" << std::endl;

	// Evolve Planet
	for(int i=0;i<max_iter;i++){
		Haumea.iter = i;
		data_file << Haumea; // write step to file
		std::cout << "Iter:\t" << i << "\tPeriod\t" << 1./(sqrt(Haumea.omega_sq)*3600) << 
		"\tMass\t" << Haumea.Mass << std::endl;
		Haumea.evolve(K0,K0_prime,rho_0); // updata
	}

	data_file.close();


/****************** Output *************************/
	std::cout << "\nComputed Mass:\t" << Haumea.Mass << std::endl;
	std::cout << "Real Mass:\t" << 	real_mass << std::endl;
	std::cout << "Mass Ratio:\t" << real_mass/Haumea.Mass << std::endl;
	//std::cout << "Legendre Test:\t" << Haumea.P_array[0][0] << "\t" << Haumea.P_array[1][0] << "\t" << Haumea.P_array[4][30] << std::endl;
	//std::cout << "Factorial Test:\t" << std::tgamma(4+1) << std::endl;
	std::cout << "Computed Total Grav Potential:\t" << Haumea.W << std::endl;
	std::cout << "Real Total Grav Potential:\t" << real_W << std::endl;
	std::cout << "Ratio:\t" << real_W/Haumea.W << std::endl;
	std::cout << "Computed A Radius:\t" << Haumea.r_array[Haumea.A_r] << std::endl;
	std::cout << "Real A Radius:\t" << a << std::endl;
	std::cout << "Computed B Radius:\t" << Haumea.r_array[Haumea.B_r] << std::endl;
	std::cout << "Real B Radius:\t" << b << std::endl;
	std::cout << "Computed Omega Squared:\t" << Haumea.omega_sq << std::endl;
	std::cout << "Period [hours]:\t" << 1./(sqrt(Haumea.omega_sq)*3600)<< std::endl;
	std::cout << "Done" << std::endl;

}