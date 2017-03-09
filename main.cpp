#include "header.h"
#include <fstream>

int main(int argc, char** argv)
{

/*************** Define Variables *******************/

	// Spherical Grid Spacing
	int n_r = 129;
	int n_mu = 33; // mu = cos(theta)
	int n_phi = 65;
	double r_max = 16./15; // maximum r-value > 1 to allow for max r off equator

	// Axes for Triaxial Ellipsoid
	double a = 1;
	double b = 1;
	double c = 1;

	// Planet Properties
	double rho_rock = 2600;
	double rho_ice = 2600;
	double Haumea_Mass = 0.;
	double Haumea_W = 0.;

	// Function Values
	int L_max = 16; // number of legendre functions to keep

/**************** Build Planet *********************/

	// Call Instance of Planet Class and Initialize it with the given composition
	Planet Haumea = Planet(n_r,n_mu,n_phi);
	Haumea.init_density(a,b,c,rho_rock,rho_ice,r_max);
	Haumea.init_P(L_max);
	Haumea.init_potential(L_max);

	// Calculate relevant quantities using methods defined in "Planet" class
	Haumea_Mass = Haumea.get_mass();
	//Haumea_W = Haumea.get_W();

/****************** Output *************************/
	std::cout << "Computed Mass:\t" << Haumea_Mass << std::endl;
	std::cout << "Real Mass:\t" << 	4*M_PI/3.*a*b*c*rho_rock << std::endl;
	std::cout << "Error:\t" << Haumea_Mass/(4*M_PI/3.*a*b*c*rho_rock) << std::endl;
	//std::cout << "Legendre Test:\t" << Haumea.P_array[0][0] << "\t" << Haumea.P_array[1][0] << "\t" << Haumea.P_array[4][30] << std::endl;
	std::cout << "Grav Potentials:\t" << Haumea.g_potential[0][0][0] << "\t" << Haumea.g_potential[128][0][0] << std::endl;
	//std::cout << "Total Grav Potential:\t" << Haumea_W << std::endl;
	std::cout << "Done" << std::endl;

}