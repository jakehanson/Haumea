#include "header.h"
#include <stdio.h>
//#include <stdexcept>

/* Define function to overload the operator "<<" such that if the input 
is std:ostream &out and &Planet_Name it prints the relevant info */
std::ostream &operator<<(std::ostream &out, Planet const &Haumea){
	double phi,r,mu;
	mu = 0.;
	for(int i=0;i<Haumea.phi_size;i++){
		phi = Haumea.phi_array[i];
		for(int k=0;k<Haumea.r_size;k++){
			r = Haumea.r_array[k];
			out << Haumea.iter << "\t" << phi << "\t" << mu << "\t" << r << "\t" << Haumea.density[i][0][k] << std::endl;
		}
	}
	return out;
}


/* Define our constructor. Pass Args to this to create an instance of Planet Class */
Planet::Planet(int n_r,int n_mu, int n_phi, int L)
{
	// Define structure members
	r_size = n_r;
	mu_size = n_mu;
	phi_size = n_phi;
	L_max = L;
	r_array = std::vector<double>(r_size);
	mu_array = std::vector<double>(mu_size);
	phi_array = std::vector<double>(phi_size);
	
	// Get angular indices for fixed points A and B (in xy plane)
	A_phi = 0;  // on x axis
	B_phi = phi_size-1;	// on y axis
	A_mu = 0;
	B_mu = 0;

	// Generate 3d arrays to hold values
	density = std::vector<std::vector<std::vector<double>>>
	(phi_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(r_size,0.0)));

	enthalpy = std::vector<std::vector<std::vector<double>>>
	(phi_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(r_size,0.0)));

	g_potential = std::vector<std::vector<std::vector<double>>>
	(phi_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(r_size,0.0)));
}


/* Define function to initialize grid and density */
void Planet::init_density(double a,double b,double c,double rho_0,double r_max)
{
	double r,mu,phi;
	double x,y,z;
	double f;

	// Initialize Grid
	for(int i=0;i<phi_size;i++){
		phi_array[i] = M_PI*i/(2*(phi_size-1.));
	}
	for(int j=0;j<mu_size;j++){
		mu_array[j] = j/(mu_size-1.);
	}
	r_array[0] = 0.;
	A_r = 0;
	B_r = 0;
	for(int k=1;k<r_size;k++){
		r_array[k] = r_max*k/(r_size-1.);

		// Check if we are at a fixed point
		if(r_array[k-1] < a and r_array[k] >= a){
			A_r = k-1; // stores radial index of fixed point A
		}
		if(r_array[k-1] < b and r_array[k] >= b){
			B_r = k-1; // stores radial index of fixed point B
		}
	}
	if(A_r == 0 or B_r == 0){
		std::cout << "ERROR! Failed to find radial index of fixed point!" << std::endl;
	}

	// Initialize Density
	for(int i=0;i<phi_size;i++){
		phi = phi_array[i];
		for(int j=0;j<mu_size;j++){
			mu = mu_array[j];
			for(int k=0;k<r_size;k++){
				r = r_array[k];

				// Populate only cells within ellipsoid
				x = r*sqrt(1-pow(mu,2))*cos(phi);
				y = r*sqrt(1-pow(mu,2))*sin(phi);
				z = r*mu;
				f = pow(x/a,2)+pow(y/b,2)+pow(z/c,2);
				if(f<=1){
					density[i][j][k] = rho_0; //dimfull
				}
				else{
					density[i][j][k] = 0.;
				}

			}
		}
	}
}


/* Define function to calculate legendre polynomials P[l][mu_index] */
void Planet::init_P(void){
	P_array = std::vector<std::vector<double>>(L_max,std::vector<double>(mu_size,0.0));
	double P_minus_2,P_minus_1,P;
	double P0,P1;
	for(int j=0;j<mu_size;j++){
		double x = mu_array[j];
		for(int l=0;l<L_max;l++){
			P0 = 1;
			P1 = x;
			if(l==0){
				P = P0;
			}else if(l == 1){
				P = P1;
			}else{
				P_minus_2 = P0;
				P_minus_1 = P1;
				for(int a=2;a<=l;a++){
					P = ((2*(a-1)+1)*x*P_minus_1 - (a-1)*P_minus_2)/a;
					P_minus_2 = P_minus_1;
					P_minus_1 = P;				
				}
			}	
			P_array[l][j] = P;
		}
	}
}


/* Define a function to calculate the gravitational potential of the planet */
/* This is from Hachisu 1986b Eq. 36 */
void Planet::get_potential(void){
	double phi,mu,r;
	double sum;  // holds the partial sum of phi at i,j,k
	double G = 6.67e-8;
	double e_m;
	for(int i=0;i<phi_size;i++){
		phi = phi_array[i];
		for(int j=0;j<mu_size;j++){
			mu = mu_array[j];
			for(int k=0;k<r_size;k++){
				r = r_array[k];
				sum = 0.;
				for(int l=0;l<L_max;l=l+2){
					for(int m=0;m<=l;m=m+2){
						if(m==0){
							e_m = 1.; 
						}else{
							e_m = 2.;
						}
						sum = sum - e_m*std::tgamma(l-m+1)/std::tgamma(l+m+1)*D3_array[l][m][k]*P_array[l][j]*cos(m*phi);
					}
				}
				g_potential[i][j][k] = G*sum;
			}
		}
	}
}
void Planet::init_D3(void){
	D3_array = 	std::vector<std::vector<std::vector<double>>>
	(L_max, std::vector< std::vector<double>>(L_max,std::vector<double>(r_size,0.0)));
	for(int k=0;k<r_size;k++){
		for(int m=0;m<L_max;m=m+2){
			for(int l=0;l<L_max;l=l+2){	
				double value=0;
				for(int s=0;s<r_size-2;s=s+2){
					value = value+1/6.*(r_array[s+2]-r_array[s])*
						(f_l(r_array[s],r_array[k],l)*D2_array[s][l][m]+4*f_l(r_array[s+1],r_array[k],l)*
						D2_array[s+1][l][m]+f_l(r_array[s+2],r_array[k],l)*D2_array[s+2][l][m]);
				}
				D3_array[l][m][k] = value;
			}
		}
	}

}
void Planet::init_D2(void){
	D2_array = 	std::vector<std::vector<std::vector<double>>>	
	(r_size, std::vector< std::vector<double>>(L_max,std::vector<double>(L_max,0.0)));
	for(int m=0;m<L_max;m=m+2){
		for(int l=0;l<L_max;l=l+2){
			for(int s=0;s<r_size;s++){
				double value = 0.;	
				for(int t=0;t<mu_size-2;t=t+2){
					value = value + 1/3.*(mu_array[t+2]-mu_array[t])*(P_array[l][t]*
						D1_array[t][s][m]+4*P_array[l][t+1]*D1_array[t+1][s][m]+
						P_array[l][t+2]*D1_array[t+2][s][m]);
				}
				D2_array[s][l][m] = value;
			}
		}
	}
}
void Planet::init_D1(void){
	D1_array = 	std::vector<std::vector<std::vector<double>>>	
	(mu_size, std::vector< std::vector<double>>(r_size,std::vector<double>(L_max,0.0)));
	for(int m=0;m<L_max;m=m+2){
		for(int s=0;s<r_size;s++){
			for(int t=0;t<mu_size;t++){
				double value = 0.;
				for(int u=0;u<phi_size-2;u=u+2){
					value = value + 2/3.*(phi_array[u+2]-phi_array[u])*(
						cos(m*phi_array[u])*density[u][t][s]+
						4*cos(m*phi_array[u+1])*density[u+1][t][s]+
						cos(m*phi_array[u+2])*density[u+2][t][s]);
				}
				D1_array[t][s][m] = value;
			}	
		}
	}
}
double Planet::f_l(double r1, double r2, int l){
	double value = 0.;
	if(r1==0 and r2==0){
		return 0.;
	}
	if(r1 <= r2){
		value = pow(r1,double(l)+2)/pow(r2,double(l)+1);
	}else{
		value = pow(r2,double(l))/pow(r1,double(l)-1.);
	}
	return value;
}


/* Define function to calculate enthalpy */
/* Hachisu 1986b Eq.8 */
void Planet::get_enthalpy(double K0, double K0_prime,double rho_0){
	double phi,mu,r;
	double minthalpy;
	double maxthalpy;

	// First get omega squared and C
	omega_sq = 2*(g_potential[A_phi][A_mu][A_r]-g_potential[B_phi][B_mu][B_r])/
		(pow(r_array[A_r],2)*(1-pow(mu_array[A_mu],2))-pow(r_array[B_r],2)*(1-pow(mu_array[B_mu],2)));
	C = g_potential[A_phi][A_mu][A_r]-1/2.*omega_sq*pow(r_array[A_r],2)*(1-pow(mu_array[A_mu],2));
	
	minthalpy =  C - g_potential[0][0][0]; // initialize min enthalpy
	maxthalpy =  C - g_potential[0][0][0]; // initialize max enthalpy

	// Now get enthalpy at each grid point
	for(int i=0;i<phi_size;i++){
		phi = phi_array[i];
		for(int j=0;j<mu_size;j++){
			mu = mu_array[j];
			for(int k=0;k<r_size;k++){
				r = r_array[k];
				enthalpy[i][j][k] = C - g_potential[i][j][k] + omega_sq*pow(r,2)*(1-pow(mu,2))/2.;
				if(enthalpy[i][j][k] < minthalpy){
					minthalpy = enthalpy[i][j][k];
				}
				if(enthalpy[i][j][k] > maxthalpy){
					maxthalpy = enthalpy[i][j][k];
				}
			}
		}
	}
	// Make sure enthalpy is 0 at fixed points
	if(enthalpy[A_phi][A_mu][A_r] > 1e-6 or enthalpy[B_phi][B_mu][B_r] > 1e-6){
		std::cout << "ERROR! Enthalpy at fixed point is non-zero!" << std::endl;
	}
}


/* Define function to use enthalpy[i][j][k] to update density[i][j][k]*/
/* Hachisu 1986b (Section 2) but for an Birch–Murnaghan equation of state */
void Planet::evolve(double K0, double K0_prime,double rho_0){

	//std::cout << "Density at A:\t" << density[A_phi][A_mu][A_r] << std::endl;
	//std::cout << "Density at B:\t" << density[B_phi][B_mu][B_r] << std::endl;
	if(std::abs(density[A_phi][A_mu][A_r]-rho_0)>0.01 or std::abs(density[B_phi][B_mu][B_r]-rho_0)>0.01){
		std::cout << "ERROR! Density at fixed point has changed!" << std::endl;
	}

	// Update Density at each point using enthalpy at each point
	for(int i=0;i<phi_size;i++){
		for(int j=0;j<mu_size;j++){
			for(int k=0;k<r_size;k++){
				density[i][j][k] = newton_raphson(enthalpy[i][j][k],density[i][j][k]/rho_0,K0,K0_prime,rho_0);
			}
		}
	}
	Mass = get_mass();
	get_potential();
	get_enthalpy(K0,K0_prime,rho_0);
}
double Planet::newton_raphson(double H, double x_guess, double K0, double K0_prime,double rho_0){
	int max_iter = 1000;
	int i;
	double tol = 1e-6;
	double x_old,x_new,value;
	if(x_guess == 0){
		x_old = 1e-6; // don't use 0 as a guess
	}else{
		x_old = x_guess;
	}
	for(i=0;i<=max_iter;i++){
		value = root_func(H,x_old,K0,K0_prime,rho_0);
		if(std::abs(value)<tol){
			x_new = x_old;
			break; // root found
		}else{
			x_new = x_old-root_func(H,x_old,K0,K0_prime,rho_0)/root_func_deriv(H,x_old,K0,K0_prime,rho_0);
			x_old= x_new;
		}
		if(i==max_iter){
			std::cout << "ERROR - UNABLE TO LOCATE ROOT!\t" << std::endl;
			x_new = -1;
			break;
		}
	}
	return x_new*rho_0;
}
/*Birch–Murnaghan equation of state*/
double Planet::root_func(double H, double x,double K0,double K0_prime,double rho_0){
	return 2*rho_0*H/(3*K0)-3/4.*(K0_prime-4)*(3/2.*pow(x,2.)-7/2.*pow(x,(4/3.))+5/2.*pow(x,(2/3.))-1/2.)-(7/4.*pow(x,(4/3.))-5/2.*pow(x,(2/3.))+3/4.);
	//return 16*rho_0*H/(3*K0)+3*K0_prime-18-pow(x,4./3)*(98-21*K0_prime)
	//		-pow(x,2./3)*(15*K0_prime-80)-pow(x,2)*(9*K0_prime-36);
}
double Planet::root_func_deriv(double H,double x,double K0,double K0_prime,double rho_0){
	return -3/4.*(K0_prime-4)*(3*pow(x,2.)-14/3.*pow(x,(1/3.))+5/3.*pow(x,(-1/3.)))-(7/3.*pow(x,(1/3.))-5/3.*pow(x,(-1/3.)));
	// return -4./3*pow(x,1./3)*(98-21*K0_prime)-2./3*pow(x,-1./3)*(15*K0_prime-80)
	// 		-2*x*(9*K0_prime-36);
}


/* Define Functions to integrate density and find the total mass */
double Planet::get_mass(void){
	double M=0;
	for(int k=0;k<r_size-2;k=k+2){
		M = M + 1/6.*(r_array[k+2]-r_array[k])*(pow(r_array[k],2)*Q2_array[k]+4*pow(r_array[k+1],2)*
			Q2_array[k+1]+pow(r_array[k+2],2)*Q2_array[k+2]);
	}
	return M;
}
void Planet::init_Q2(void){
	Q2_array = std::vector<double>(r_size,0.0);
	for(int k=0;k<r_size;k++){	
		double value = 0.;
		for(int j=0;j<mu_size-2;j=j+2){
			value = value + 1/3.*(mu_array[j+2]-mu_array[j])*(Q1_array[j][k]+4*Q1_array[j+1][k]+Q1_array[j+2][k]);
		}
		Q2_array[k] = value;
	}
}
void Planet::init_Q1(void){
	Q1_array = std::vector<std::vector<double>>(mu_size,std::vector<double>(r_size,0.0));
	for(int k=0;k<r_size;k++){	
		for(int j=0;j<mu_size;j++){
			double value = 0.;
			for(int i=0;i<phi_size-2;i=i+2){
				value = value + 4/6.*(phi_array[i+2]-phi_array[i])*(density[i][j][k]+4*density[i+1][j][k]+density[i+2][j][k]);
			}
			Q1_array[j][k] = value;
		}
	}
}


/* Define function to caculate the gravitational potential of the planet */
double Planet::get_W(void){
	double W=0;
	for(int k=0;k<r_size-2;k=k+2){
		W = W - 1/12.*(r_array[k+2]-r_array[k])*(pow(r_array[k],2)*
			S2_array[k]+4*pow(r_array[k+1],2)*S2_array[k+1]+pow(r_array[k+2],2)*S2_array[k+2]);
	}
	return W;
}
void Planet::init_S2(void){
	S2_array = std::vector<double>(r_size,0.0);
	for(int k=0;k<r_size;k++){
		double value = 0.;
		for(int j=0;j<mu_size-2;j=j+2){
			value = value + 1/3.*(mu_array[j+2]-mu_array[j])*(S1_array[j][k]+4*S1_array[j+1][k]+S1_array[j+2][k]);
		}
		S2_array[k] = value;
	}
}
void Planet::init_S1(void){
	S1_array = std::vector<std::vector<double>>(mu_size,std::vector<double>(r_size,0.0));
	for(int k=0;k<r_size;k++){
		for(int j=0;j<mu_size;j++){
			double value = 0.;
			for(int i=0;i<phi_size-2;i=i+2){
				value = value + 2/3.*(phi_array[i+2]-phi_array[i])*(density[i][j][k]*
					g_potential[i][j][k]+4*density[i+1][j][k]*g_potential[i+1][j][k]+density[i+2][j][k]*g_potential[i+2][j][k]);
			}
			S1_array[j][k] = value;
		}
	}
}


