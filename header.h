#pragma once
#include <vector>
#include <iostream>
#include <random>

/* Structure to hold Planet Properties*/
struct Planet
{

	//std::array<double> density; // Array to hold density
	int r_size,mu_size,phi_size;

	std::vector<double> r_array;
	std::vector<double> mu_array;
	std::vector<double> phi_array;

	std::vector<std::vector<std::vector<double>>> density; // 3d array for desnity
	std::vector<std::vector<std::vector<double>>> enthalpy; // 3d array for desnity
	std::vector<std::vector<std::vector<double>>> g_potential; // 3d array for desnity

	Planet(int n_r,int n_mu, int n_phi);  // Signature for constructor. Function to create instance of class
	void initialize(double a, double b, double c, double rho_rock, double rho_ice,double r_max); // Function to initialize values
	double Q_1(int j, int k);
	double Q_2(int k);
	double get_mass(void);

};