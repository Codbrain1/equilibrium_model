#include<fstream>
#include<vector>
#include<random>
#include<iostream>
#include<cmath>
#include"Particle.h"
#include <iomanip>
#include<thread>
double random(double beg, double end)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(beg, end);
	return dist(gen);
}

namespace PPm = Partcile_Particle_model;

std::vector<double> sistem_E;
std::vector<double> sistem_P;
std::vector<double> sistem_M;
std::vector<double> sistem_t;
std::vector<double> sistem_E_k;
std::vector<double> sistem_E_p;
std::vector<double> sistem_r;

void calculating(std::vector<PPm::Particle>& particles, int k, int i_0, int i_1)
{
	std::vector<vec> u_i(i_1 - i_0, vec(0, 0, 0));
	for (size_t i = i_0; i < i_1; i++)
	{
		u_i[i - i_0] = particles[i].v;
		particles[i].v = particles[i].v + particles[i].F * PPm::dt;
		particles[i].r = particles[i].r + (particles[i].v + u_i[i - i_0]) * PPm::dt * 0.5;
	}


	for (int i = i_0; i < i_1; i++)
	{
		particles[i].F = vec(0, 0, 0);
	}

	for (size_t i = i_0; i < i_1; i++)
	{
		for (size_t j = 0; j < particles.size(); j++)
		{
			if (i != j)
			{
				particles[i].F = particles[i].F + PPm::F(particles[i], particles[j]);
			}
		}
	}

	for (size_t i = i_0; i < i_1; i++)
	{
		particles[i].v = (particles[i].v + u_i[i - i_0]) * 0.5 + particles[i].F * 0.5 * PPm::dt;
	}
}

void calculate_conversation_laws(std::vector<PPm::Particle>& particles, int k, int i_0, int i_1)
{
	for (size_t i = i_0; i < i_1; i++)
	{
		particles[i].E = 0;
		particles[i].P = vec(0, 0, 0);
		particles[i].M = vec(0, 0, 0);
	}

	for (int i = i_0; i < i_1; i++)
	{
		particles[i].E = particles[i].m * particles[i].v.module_2() * 0.5;	// set kinetical energy
		particles[i].P = particles[i].v * particles[i].m;					// set impulse
		particles[i].M = particles[i].r * particles[i].P;					// set moment impulse

		sistem_E_k[k] += particles[i].m * particles[i].v.module_2() * 0.5;

		for (int j = 0; j < particles.size(); j++)
		{
			if (i != j)
			{
				vec r_ij = particles[j].r - particles[i].r;
				particles[i].E -= PPm::G * particles[i].m * particles[j].m / r_ij.module();	//set potential energy
				sistem_E_p[k] -= PPm::G * particles[i].m * particles[j].m / r_ij.module();

			}
		}
		sistem_E[k] += particles[i].E;
		sistem_P[k] = sistem_P[k] + particles[i].P.module();
		sistem_M[k] = sistem_M[k] + particles[i].M.module();
	}
}

void set_initial_conditions(std::vector<PPm::Particle>&ps)
{
	//setting the initial position
		// 
		//====================================================
	ps[0].r.x = 0;
	ps[0].r.y = 0;
	ps[0].r.z = 0;
	ps[0].m = 1;
	ps[1].r.x = (PPm::R_max)*cos(0);
	ps[1].r.y = (PPm::R_max)*sin(0);
	ps[1].r.z = 0;
	ps[1].m = 1.0 / 333000.0;

	ps[2].r.x = (PPm::R_max*0.5)*cos(PPm::PI);
	ps[2].r.y = (PPm::R_max*0.5)*sin(PPm::PI);
	ps[2].r.z = 0;
	ps[2].m = 1.0 / 333000.0;

	//setting the initial velocity
	// 
	//==========================================================
	ps[1].F = F(ps[1], ps[0]) + ps[1].F;									//calculate force for litle particle
	ps[1].F = F(ps[1], ps[2]) + ps[1].F;

	ps[2].F = F(ps[2], ps[0]) + ps[1].F;									//
	ps[2].F = F(ps[2], ps[1]) + ps[1].F;									//

	double v_asimutal = sqrt(ps[1].r.module() * ps[1].F.module());
	double phi = atan2(ps[1].r.y, ps[1].r.x);								// calculate velocity for litle particle
	double r = ps[1].r.module();											// (cylindrical coordinates)

	double v_asimutal1 = sqrt(ps[2].r.module() * ps[2].F.module());
	double phi1 = atan2(ps[2].r.y, ps[2].r.x);
	double r1 = ps[2].r.module();

	//setting the initial velocity
	// 
	//==========================================================
	ps[1].v.x = -r * v_asimutal * sin(phi);									//
	ps[1].v.y = r * v_asimutal * cos(phi);									// calculate velocity for litle particle
	ps[1].v.z = 0;															// (kartesian coordinates)

	ps[2].v.x = -r1 * v_asimutal1 * sin(phi1);
	ps[2].v.y = r1 * v_asimutal1 * cos(phi1);
	ps[2].v.z = 0;
}

int main()
{


	std::vector<PPm::Particle> particles(PPm::N); //array particals in model
	//setting the initial conditions
	// 
	//====================================================================================================
	set_initial_conditions(particles);

	std::cout << "set initial conditions\n";
	
	sistem_E.push_back(0);
	sistem_P.push_back(0);
	sistem_M.push_back(0);
	sistem_t.push_back(0);
	sistem_E_k.push_back(0);
	sistem_E_p.push_back(0);
	sistem_r.push_back(0);
	sistem_r[0] = particles[1].r.module();
	calculate_conversation_laws(particles, 0, 1, particles.size());
	std::cout << "set intitial conversation laws\n";
	//==========================================================

// integration of differential equations
// 
//========================================================================================================================
	
	std::cout << "write into file\n";

	std::ofstream positions("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация_измерений\\positions.txt");
	positions << PPm::N << std::endl;
	positions << 0 << std::endl;
	for (int i = 0; i < particles.size(); i++)
	{
		positions << std::setprecision(10) << particles[i].r.x << " " << particles[i].r.y << " " << particles[i].r.z << " " << particles[i].v.x << " " << particles[i].v.y << " " << particles[i].v.z << std::endl;
	}

	int k = 1;
	int b = 1;
	std::cout << "start calculating\n";

	for (double t = PPm::t_0 + PPm::dt; t <= PPm::t_1; t += PPm::dt)
	{
		
		calculating(particles, k, 1, particles.size());

		if (b % PPm::div == 0) {
			sistem_E.push_back(0);
			sistem_P.push_back(0);
			sistem_M.push_back(0);
			sistem_t.push_back(t);
			sistem_E_k.push_back(0);
			sistem_E_p.push_back(0);
			sistem_r.push_back(0);
			sistem_r[k] = particles[1].r.module();
			calculate_conversation_laws(particles, sistem_E.size()-1, 1, particles.size());
			
		}
		if (b % PPm::div == 0) {
			positions << t << std::endl;
			for (int i = 0; i < particles.size(); i++)
			{
				positions << std::setprecision(10) << particles[i].r.x << " " << particles[i].r.y << " " << particles[i].r.z << " " << particles[i].v.x << " " << particles[i].v.y << " " << particles[i].v.z << std::endl;
			}
			k++;
		}
		b++;

	}
	std::cout << "end calculating\n";
	std::cout << "write conwersation laws into file\n";
	std::ofstream conversation_laws("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация_измерений\\measurements.txt");
	for (int i = 0; i < sistem_E.size(); i++)
	{
		conversation_laws <<std::fixed<< std::setprecision(15) << sistem_t[i] << " " << sistem_E[i]-sistem_E[0] << " " << sistem_P[i] << " " << sistem_M[i] -sistem_M[0]<< " " << sistem_E_k[i]<< " " << sistem_E_p[i]<< " " << sistem_r[i]-sistem_r[0]<< std::endl;
	}
	//====================================================================================================
}