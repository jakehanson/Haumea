#include "header.h"
#include <stdio.h>
#include <cmath>
//#include <stdexcept>


/* Define our constructor. Pass Args to this to create an instance of Planet Class */
Planet::Planet(int n_r,int n_mu, int n_phi, int L)
{
	/* Store information as member of structure */
	r_size = n_r;
	mu_size = n_mu;
	phi_size = n_phi;
	L_max = L;
	r_array = std::vector<double>(r_size);
	mu_array = std::vector<double>(mu_size);
	phi_array = std::vector<double>(phi_size);

	/* Generate Grid */	
	density = std::vector<std::vector<std::vector<double>>>
	(phi_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(r_size,0.0)));

	enthalpy = std::vector<std::vector<std::vector<double>>>
	(phi_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(r_size,0.0)));

	g_potential = std::vector<std::vector<std::vector<double>>>
	(phi_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(r_size,0.0)));

}


/* Define function to initialize values inside grid */
void Planet::init_density(double a,double b,double c,double rho_rock,double r_max)
{
	double r,mu,phi;
	double x,y,z;
	double f;

	/* Initialize Grid */
	for(int i=0;i<phi_size;i++){
		phi_array[i] = M_PI*i/(2*(phi_size-1.));
	}
	for(int j=0;j<mu_size;j++){
		mu_array[j] = j/(mu_size-1.);
	}
	for(int k=0;k<r_size;k++){
		r_array[k] = r_max*k/(r_size-1.);
	}

	/* Initialize Density */
	for(int i=0;i<phi_size;i++){
		phi = phi_array[i];
		for(int j=0;j<mu_size;j++){
			mu = mu_array[j];
			for(int k=0;k<r_size;k++){
				r = r_array[k];

				x = r*sqrt(1-pow(mu,2))*sin(phi);
				y = r*sqrt(1-pow(mu,2))*cos(phi);
				z = r*mu;

				f = pow(x/a,2)+pow(y/b,2)+pow(z/c,2);
				
				if(f<=1){
					density[i][j][k] = rho_rock; //dimfull
				}
				else{
					density[i][j][k] = 0.;
				}
			}
		}
	}

}

/* Define function to calculate legendre polynomials */
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
void Planet::init_potential(void){
	double phi,mu,r;
	double sum;  // holds the partial sum of phi at i,j,k
	double G = 6.67e-8;
	for(int i=0;i<phi_size;i++){
		phi = phi_array[i];
		for(int j=0;j<mu_size;j++){
			mu = mu_array[j];
			for(int k=0;k<r_size;k++){
				r = r_array[k];
				sum = -1*D3_array[0][0][k]; // initialize with first term (l=0 m=0)
				for(int l=2;l<L_max;l=l+2){
					sum = sum - (D3_array[l][0][k]*P_array[l][j]+1/12.*D3_array[l][2][k]*P_array[l][j]*cos(2.*phi));
				}
				g_potential[i][j][k] = G*sum;
				//std::cout << "Grav Potential:\t" << sum << std::endl;

			}
		}
	}
}

void Planet::init_D3(void){
	D3_array = 	std::vector<std::vector<std::vector<double>>>
	(L_max, std::vector< std::vector<double>>(3,std::vector<double>(r_size,0.0)));
	for(int k=0;k<r_size;k++){
		for(int m=0;m<=2;m=m+2){
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
	(r_size, std::vector< std::vector<double>>(L_max,std::vector<double>(3,0.0)));
	for(int m=0;m<=2;m=m+2){
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
	(mu_size, std::vector< std::vector<double>>(r_size,std::vector<double>(3,0.0)));
	for(int m=0;m<=2;m=m+2){
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


/* Define a function to check the gravitational potential of the planet */
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

