#pragma once
#include <vector>
#include <iostream>
#include <random>

/* Structure to hold Planet Properties*/
struct Planet
{

	//std::array<double> density; // Array to hold density

	/* Grid Properties */
	int r_size,mu_size,phi_size,L_max;
	std::vector<double> r_array;
	std::vector<double> mu_array;
	std::vector<double> phi_array;

	std::vector<std::vector<double>> P_array; // stores legendre functions
	std::vector<std::vector<double>> Q1_array;
	std::vector<double> Q2_array;
	std::vector<double> S2_array;
	std::vector<std::vector<double>> S1_array; 

	/* Arrays to hold calculated values */
	std::vector<std::vector<std::vector<double>>> density; // 3d array for density
	std::vector<std::vector<std::vector<double>>> enthalpy; // 3d array for enthalpy
	std::vector<std::vector<std::vector<double>>> g_potential; // 3d array for G potential

	std::vector<std::vector<std::vector<double>>> D1_array; // 3d array for G potential
	std::vector<std::vector<std::vector<double>>> D2_array; // 3d array for G potential
	std::vector<std::vector<std::vector<double>>> D3_array; // 3d array for G potential



	/* Define Methods */
	Planet(int n_r,int n_mu, int n_phi,int L_max);  // Signature for constructor. Function to create instance of class
	void init_density(double a, double b, double c, double rho_rock,double r_max); // Function to initialize values

	double get_mass(void);
	void init_Q1(void);
	void init_Q2(void);

	void init_potential(void);
	void init_D3(void);
	void init_D2(void);
	void init_D1(void);
	double f_l(double r1, double r2, int l);

	double get_W(void);
	void init_S2(void);
	void init_S1(void);
	void init_P(void); // Legendre Polynomial Functions
	

};

