#include "header.h"
#include <fstream>

int main(int argc, char** argv)
{

/*************** Define Variables *******************/

	// Spherical Grid Spacing
	int n_r = 129;
	//int n_r = 65;
	//int n_mu = 17;
	//int n_phi = 33;
	int n_mu = 33; // mu = cos(theta)
	int n_phi = 65;
	double r_max = 16./15; // maximum r-value > 1 to allow for max r off equator

	// Axes for Triaxial Ellipsoid
	double a = 1;
	double b = 1;
	double c = 1;

	// Planet Properties
	double rho_rock = 2.6;
	double rho_ice = 2.6;
	double Haumea_Mass = 0.;
	double Haumea_W = 0.;

	double R_0 = 2.9e8; // 2900km in cm
	double rho_0 = 1; // typical density scale is rock
	double G = 6.67e11;

	// Function Values
	int L_max = 16; // number of legendre functions to keep

/**************** Build Planet *********************/

	// Call Instance of Planet Class and Initialize it with the given composition
	Planet Haumea = Planet(n_r,n_mu,n_phi);
	Haumea.init_density(a,b,c,rho_0,rho_ice,r_max);
	Haumea.init_P(L_max);
	Haumea.init_potential(L_max);

	// Calculate relevant quantities using methods defined in "Planet" class
	Haumea_Mass = Haumea.get_mass();
	Haumea_W = Haumea.get_W();

/****************** Output *************************/
	std::cout << "Computed Mass:\t" << Haumea_Mass << std::endl;
	std::cout << "Real Mass:\t" << 	4*M_PI/3.*a*b*c*rho_rock << std::endl;
	std::cout << "Error:\t" << Haumea_Mass/(4*M_PI/3.*a*b*c*rho_rock) << std::endl;
	//std::cout << "Legendre Test:\t" << Haumea.P_array[0][0] << "\t" << Haumea.P_array[1][0] << "\t" << Haumea.P_array[4][30] << std::endl;
	std::cout << "Grav Potentials:\t" << Haumea.g_potential[0][0][0] << "\t" << Haumea.g_potential[128][0][0] << std::endl;
	std::cout << "Computed Total Grav Potential:\t" << Haumea_W*rho_0*G*pow(R_0,2) << std::endl;
	std::cout << "Real Total Grav Potential:\t" << 3*G*pow(4/3.*M_PI*pow(R_0,3)*rho_rock,2)/(5*R_0) << std::endl;
	std::cout << "Error:\t" << Haumea_W*rho_0*G*pow(R_0,2)/(3*G*pow(4/3.*M_PI*pow(R_0,3)*rho_rock,2)/(5*R_0)) << std::endl;
	std::cout << "Done" << std::endl;

}