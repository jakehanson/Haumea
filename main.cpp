#include "header.h"
#include <fstream>

int main(int argc, char** argv)
{
	/* Define Variables */
	int n_r = 129;
	int n_mu = 33;
	int n_phi = 65;
	double a = 1;
	double b = 1;
	double c = 1;
	double r_max = 16./15;
	//double r_max = 1.0;
	double rho_rock = 2600;
	double rho_ice = 2600;

	Planet Haumea = Planet(n_r,n_mu,n_phi);
	Haumea.initialize(a,b,c,rho_rock,rho_ice,r_max);


	/* Run Code */
	std::cout << "R array" << '\t' << Haumea.r_array[120] << '\t' << Haumea.r_array[119] << std::endl;

	std::cout << 4*M_PI/3.*a*b*c*rho_rock << std::endl;
	std::cout << Haumea.get_mass() << std::endl;
	std::cout << Haumea.density[121][32][64] << std::endl;
	std::cout << Haumea.density[120][32][4] << std::endl;
	std::cout << "Done" << std::endl;

}