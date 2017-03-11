#include "header.h"
#include <fstream>

int main(int argc, char** argv)
{

/*************** Define Variables *******************/

	// Spherical Grid Spacing
	int n_r = 129;
	int n_mu = 129; // mu = cos(theta)
	int n_phi = 129;
//	int n_mu = 129;
//	int n_phi = 129;
//	int n_r = 65;
//	int n_mu = 17;
//	int n_phi = 33;
//	int n_r = 129*5-1;
//	int n_mu = 129*4-1;
//	int n_phi = 129*4-1;
	double r_max = 16./15.*9.6e7; // maximum r-value > R_e to allow for max r off equator

	// Axes for Triaxial Ellipsoid
	double a = 9.6e7;
//	double b = 9.6e7;
//	double c = 9.6e7;
	double b = 7.7e7;
	double c = 4.95e7;

	// Planet Properties
	double rho_rock = 2.6;
	double rho_ice = 2.6;
	double Haumea_Mass = 0.;
	double Haumea_W = 0.;

	double rho_0 = 1.0; // typical density scale is rock
	double G = 6.67e-8;
	double real_mass = 4*M_PI/3.*a*b*c*rho_rock;
//	double real_W = 3*G*pow(real_mass,2)/(5.*a);
	double real_W = 3*G*pow(real_mass,2)*2.714/1e9;

	// Function Values
	int L_max = 32; // number of legendre functions to keep

/**************** Build Planet *********************/

	// Call Instance of Planet Class and Initialize it with the given composition
	Planet Haumea = Planet(n_r,n_mu,n_phi,L_max);
	Haumea.init_density(a,b,c,rho_rock,r_max);
	Haumea.init_P();
	Haumea.init_D1();
	Haumea.init_D2();
	Haumea.init_D3();
	Haumea.init_potential();
	Haumea.init_Q1();
	Haumea.init_Q2();
	Haumea.init_S1();
	Haumea.init_S2();

	// Calculate relevant quantities using methods defined in "Planet" class
	Haumea_Mass = Haumea.get_mass();
	//Haumea_Mass = Haumea_Mass*rho_0*pow(R_0,3);
	Haumea_W = Haumea.get_W();
//	Haumea_W = Haumea_W*rho_0*G*Haumea_Mass*pow(R_0,2);
//	Haumea_W = Haumea_W*G*pow(R_0,5)*pow(rho_max,2);

/****************** Output *************************/
	std::cout << "Computed Mass:\t" << Haumea_Mass << std::endl;
	std::cout << "Real Mass:\t" << 	real_mass << std::endl;
	std::cout << "Error:\t" << real_mass/Haumea_Mass << std::endl;
	//std::cout << "Legendre Test:\t" << Haumea.P_array[0][0] << "\t" << Haumea.P_array[1][0] << "\t" << Haumea.P_array[4][30] << std::endl;
	std::cout << "Computed Total Grav Potential:\t" << Haumea_W << std::endl;
	std::cout << "Real Total Grav Potential:\t" << real_W << std::endl;
	std::cout << "Error:\t" << real_W/Haumea_W << std::endl;
	std::cout << "Done" << std::endl;

}