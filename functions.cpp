#include "header.h"
#include <stdio.h>
#include <cmath>
//#include <stdexcept>

/* Define our constructor. Pass Args to this to create an instance of Planet Class */
Planet::Planet(int n_r,int n_mu, int n_phi)
{
	/* Store information as member of structure */
	r_size = n_r;
	mu_size = n_mu;
	phi_size = n_phi;
	r_array = std::vector<double>(r_size);
	mu_array = std::vector<double>(mu_size);
	phi_array = std::vector<double>(phi_size);

	/* Generatre Grid */	
	density = std::vector<std::vector<std::vector<double>>>
	(r_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(phi_size,0.0)));

	enthalpy = std::vector<std::vector<std::vector<double>>>
	(r_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(phi_size,0.0)));

	g_potential = std::vector<std::vector<std::vector<double>>>
	(r_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(phi_size,0.0)));

}

/* Define function to initialize values inside grid */
void Planet::initialize(double a,double b,double c,double rho_rock, double rho_ice,double r_max)
{
	double r,mu,phi;
	double x,y,z;
	double f;

	/* Initialize Grid */
	for(int i=0;i<phi_size;i++){
		phi_array[i] = M_PI*i/(phi_size-1);
	}
	for(int j=0;j<mu_size;j++){
		mu_array[j] = j/(mu_size-1);
	}
	for(int k=0;k<r_size;k++){
		r_array[k] = r_max*k/(r_size-1);
	}

	/* Initialize Density */
	for(int k=0;k<phi_size;k++){
		phi = phi_array[k];
		for(int j=0;j<mu_size;j++){
			mu = mu_array[j];
			for(int i=0;i<r_size;i++){
				r = r_array[i];

				x = r*sqrt(1-pow(mu,2))*sin(phi);
				y = r*sqrt(1-pow(mu,2))*cos(phi);
				z = r*mu;

				f = pow(x/a,2)+pow(y/b,2)+pow(z/c,2);
				
				if(f<=1){
					density[i][j][k] = rho_rock;
				}
				else{
					density[i][j][k] = 0;
				}
			}
		}
	}

}

double Planet::Q_1(int j,int k){
	double Q_1 = 0.;
	for(int i=0;i<phi_size-2;i=i+2){
		Q_1 = Q_1 + 1/3.*(phi_array[i+2]-phi_array[i])*(density[k][j][i]+4*density[k][j][i]+density[k][j][i]);
	}
	return Q_1;
}

double Planet::Q_2(int k){
	double Q_2 = 0.;
	for(int j=0;j<mu_size-2;j=j+2){
		Q_2 = Q_2 + 1/3.*(mu_array[j+2]-mu_array[j])*(Q_1(j,k)+4*Q_1(j+1,k)+Q_1(j+2,k));
	}
	return Q_2;
}

/* Define Function to integrate density and find the total mass */
double Planet::get_mass(void){
	double M=0;
	for(int k=0;k<r_size-2;k=k+2){
		M = M + 1/6.*(r_array[k+2]-r_array[k])*(pow(r_array[k],2)*Q_2(k)+4*pow(r_array[k+1],2)*Q_2(k+1)+pow(r_array[k+2],2)*Q_2(k+2));
	}
	return M;	
}
