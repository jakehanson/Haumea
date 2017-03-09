#pragma once
#include <vector>
#include <iostream>
#include <random>

/* Structure to hold Planet Properties*/
struct Planet
{

	//std::array<double> density; // Array to hold density

	/* Grid Properties */
	int r_size,mu_size,phi_size;
	std::vector<double> r_array;
	std::vector<double> mu_array;
	std::vector<double> phi_array;

	std::vector<std::vector<double>> P_array; // stores legendre functions

	/* Arrays to hold calculated values */
	std::vector<std::vector<std::vector<double>>> density; // 3d array for density
	std::vector<std::vector<std::vector<double>>> enthalpy; // 3d array for enthalpy
	std::vector<std::vector<std::vector<double>>> g_potential; // 3d array for G potential

	/* Define Methods */
	Planet(int n_r,int n_mu, int n_phi);  // Signature for constructor. Function to create instance of class
	void init_density(double a, double b, double c, double rho_rock, double rho_ice,double r_max); // Function to initialize values

	double get_mass(void);
	double Q_2(int k);
	double Q_1(int j, int k);
	
	void init_potential(int L_max);
	double D_3(int l, int m, int k);
	double D_2(int s, int l, int m);
	double D_1(int t, int s, int m);
	double f_l(double r1, double r2, int l);

	double get_W(void);
	double S_2(int k);
	double S_1(int j, int k);
	void init_P(int L_max); // Legendre Polynomial Functions
	

};

