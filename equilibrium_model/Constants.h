#pragma once
namespace Partcile_Particle_model {
	//Constants:
	//======================================================================
	//const double G = 6.67 * 10e-11;
	const double G = 1;
	const double r_c = 0.005*0.005; // в данном случае используется r_c*r_c
	const double sigma_0=0.03;
	const double r_alpha = 2.0/6.0;
	//initial conditions:
	//======================================================================
	const double t_0 = 0, t_1 = 1;
	double dt = 0.001;
	double R_max = 1;
	double M = 1;
	double r_0 = 0;
	double dr = R_max / 30;


}