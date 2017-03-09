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

	/* Generate Grid */	
	density = std::vector<std::vector<std::vector<double>>>
	(r_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(phi_size,0.0)));

	enthalpy = std::vector<std::vector<std::vector<double>>>
	(r_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(phi_size,0.0)));

	g_potential = std::vector<std::vector<std::vector<double>>>
	(r_size, std::vector< std::vector<double>>(mu_size,std::vector<double>(phi_size,0.0)));

}

/* Define function to initialize values inside grid */
void Planet::init_density(double a,double b,double c,double rho_rock, double rho_ice,double r_max)
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
					density[i][j][k] = 0.;
				}
			}
		}
	}

}

/* Define function to calculate legendre polynomials */
void Planet::init_P(int L_max){

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
void Planet::init_potential(int L_max){

	double phi,mu,r;
	double sum;  // holds the partial sum of phi at i,j,k

	for(int k=0;k<phi_size;k++){
		phi = phi_array[k];
		for(int j=0;j<mu_size;j++){
			mu = mu_array[j];
			for(int i=0;i<r_size;i++){
				r = r_array[i];
				sum = 0.;
				for(int l=0;l<L_max;l=l+2){
					sum = sum - (D_3(l,0,k)*P_array[l][j]+1/12.*D_3(l,2,k)*P_array[l][j]*cos(2.*phi));
				}
				g_potential[i][j][k] = sum;
			}
		}
	}
}

double Planet::D_3(int l, int m, int k){
	double value=0;
	for(int s=0;s<r_size-2;s=s+2){
		value = value+1/6.*(r_array[s+2]-r_array[s])*(f_l(r_array[s],r_array[k],l)*D_2(s,l,m)+4*f_l(r_array[s+1],r_array[k],l)*D_2(s+1,l,m)+f_l(r_array[s+2],r_array[k],l)*D_2(s+2,l,m));
	}
	return value;
}
double Planet::D_2(int s, int l, int m){
	double value=0.;
	for(int t=0;t<mu_size-2;t=t+2){
		value = value + 2/3.*(mu_array[t+2]-mu_array[t])*(P_array[l][t]*
			D_1(t,s,m)+4*P_array[l][t+1]*D_1(t+1,s,m)+P_array[l][t+2]*D_1(t+2,s,m));
	}
	return value;
}
double Planet::D_1(int t, int s, int m){
	double value=0.;
	for(int u=0;u<phi_size-2;u=u+2){
		value = value + 1/3.*(phi_array[u+2]-phi_array[u])*(
			cos(m*phi_array[u])*density[s][t][u]+
			4*cos(m*phi_array[u+1])*density[s][t][u+1]+
			cos(m*phi_array[u+2])*density[s][t][u+2]);
	}
	return value;
}
double Planet::f_l(double r1, double r2, int l){
	double value=0.;
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
		M = M + 1/6.*(r_array[k+2]-r_array[k])*(pow(r_array[k],2)*Q_2(k)+4*pow(r_array[k+1],2)*Q_2(k+1)+pow(r_array[k+2],2)*Q_2(k+2));
	}
	return M;	
}
double Planet::Q_2(int k){
	double value = 0.;
	for(int j=0;j<mu_size-2;j=j+2){
		value = value + 1/3.*(mu_array[j+2]-mu_array[j])*(Q_1(j,k)+4*Q_1(j+1,k)+Q_1(j+2,k));
	}
	return value;
}
double Planet::Q_1(int j,int k){
	double value = 0.;
	for(int i=0;i<phi_size-2;i=i+2){
		value = value + 4/6.*(phi_array[i+2]-phi_array[i])*(density[k][j][i]+4*density[k][j][i+1]+density[k][j][i+2]);
	}
	return value;
}
/* Define a function to check the gravitational potential of the planet */
double Planet::get_W(void){
	double W=0;
	for(int k=0;k<r_size-2;k=k+2){
		W = W - 1/12.*(r_array[k+2]-r_array[k])*(pow(r_array[k],2)*
			S_2(k)+4*pow(r_array[k+1],2)*S_2(k+1)+pow(r_array[k+2],2)*S_2(k+2));
	}
	return W;
}
double Planet::S_2(int k){
	double value = 0.;
	for(int j=0;j<mu_size-2;j=j+2){
		value = value + 1/3.*(mu_array[j+2]-mu_array[j])*(S_1(j,k)+4*S_1(j+1,k)+S_1(j+2,k));
	}
	return value;
}
double Planet::S_1(int j, int k){
	double value = 0.;
	for(int i=0;i<phi_size-2;i=i+2){
		value = value + 1/3.*(phi_array[i+2]-phi_array[i])*(density[k][j][i]*
			g_potential[k][j][i]+4*density[k][j][i+1]*g_potential[k][j][i+1]+density[k][j][i+2]*g_potential[k][j][i+2]);
	}
	return value;
}

