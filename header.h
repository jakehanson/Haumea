#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>


/* Structure to hold Planet Properties*/
struct Planet
{

/****************** Simulation Properties ******************/

	int r_size,mu_size,phi_size,L_max;
	int iter; // number of evolution iterations
	int A_phi,A_r,A_mu,B_phi,B_r,B_mu; // indices for fixed points
	std::vector<double> r_array; // r positions
	std::vector<double> mu_array; // mu=cos(theta) w/ theta measured from z axis)
	std::vector<double> phi_array; // angle from x axis


/******************** Planet Properties ********************/

	double Mass; // total mass
	double W; // total grav potential (binding energy)
	double C; // boundary condition (equation 17, Hachisu 1986)
	double omega_sq; // angular velocity squared
	std::vector<std::vector<std::vector<double>>> density; // 3d array for density
	std::vector<std::vector<std::vector<double>>> g_potential; // 3d array for G potential
	std::vector<std::vector<std::vector<double>>> enthalpy; // 3d array for enthalpy


/********** Algorithm and Finite Difference Arrays *********/

	std::vector<std::vector<double>> Q1_array;	// 2d array for calculating M
	std::vector<double> Q2_array;	// 1d array for calculating M
	std::vector<double> S2_array;	// 1d array for calculating W
	std::vector<std::vector<double>> S1_array; 	// 2d array for calculating W
	std::vector<std::vector<double>> D1_array; 	// 2d array for calculating potential
	std::vector<double> D2_array;	// 1d array for calculating potential
	std::vector<double> fl; // 1d array for calcuating potential

/************************* Methods *************************/

	Planet(int n_r,int n_mu, int n_phi,int L_max);  // Signature for constructor. Function to create instance of class
	void init_density(double a, double b, double c, double rho_0,double r_max); // Initialize grid and density
	void evolve(double K0, double K0_prime,double rho_0, int L_max); // update method
	void get_enthalpy(double K0, double K0_prime,double rho_0); // function to invert density to find enthalpy
	double newton_raphson(double H, double x_guess, double K0, double K0_prime,double rho_0); // root finder for enthalpy
	double root_func(double H, double x,double K0,double K0_prime,double rho_0); // function to get root of
	double root_func_deriv(double H,double x,double K0,double K0_prime,double rho_0); // derivative of root func

	//Mass and Supplemental Functions
	double get_mass(void);
	void init_Q1(void);
	void init_Q2(void);

	//Potential and Supplemental Functions
	void get_potential(int L_max);
	double factorial(int n);
	double plgndr(int l, int m, double x);

	// Binding Energy and Supplemental Functions
	double get_W(void);
	void init_S2(void);
	void init_S1(void);

};

/* Overload "<<" to easily write density data to file */
std::ostream &operator<<(std::ostream &out, Planet const &Haumea);

