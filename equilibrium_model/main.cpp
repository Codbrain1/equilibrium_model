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

	for (size_t i = i_0; i < i_1; i++)
	{
		particles[i].E = particles[i].m * particles[i].v.module_2() * 0.5;
		particles[i].P = particles[i].v * particles[i].m;
		particles[i].M = particles[i].r * particles[i].P;
		for (size_t j = 0; j < particles.size(); j++)
		{
			if (i != j)
			{
				vec r_ij = particles[j].r - particles[i].r;
				particles[i].E -= 0.5 * PPm::G * particles[i].m * particles[j].m / sqrt(r_ij.module_2() + PPm::r_c);

			}
		}
		sistem_E[k] += particles[i].E;
		sistem_P[k] = sistem_P[k] + particles[i].P.module();
		sistem_M[k] = sistem_M[k] + particles[i].M.module();
	}

}

int main()
{


	std::vector<PPm::Particle> particles(PPm::N); //array particals in model
	//setting the initial conditions
	// 
	//====================================================================================================

		//setting the initial position
		// 
		//====================================================
	particles[0].r.x = 0;
	particles[0].r.y = 0;
	particles[0].r.z = 0;
	particles[0].m = 1;
	particles[1].r.x = (PPm::R_max) * cos(PPm::PI);
	particles[1].r.y = (PPm::R_max) * sin(PPm::PI);
	particles[1].r.z = 0;
	particles[1].m = 1.0 / 333000.0;
	sistem_E.push_back(0);
	sistem_P.push_back(0);
	sistem_M.push_back(0);
	sistem_t.push_back(0);

	//setting the initial velocity
	// 
	//==========================================================


	particles[1].F = F(particles[1], particles[0]) + particles[1].F;	//calculate force for litle particle

	double v_asimutal = sqrt(particles[1].r.module() * particles[1].F.module());	//
	double phi = atan2(particles[1].r.y, particles[1].r.x);							// calculate velocity for litle particle
	double r = particles[1].r.module();												// (cylindrical coordinates)

	particles[1].v.x = -r * v_asimutal * sin(phi);										//
	particles[1].v.y = r * v_asimutal * cos(phi);									// calculate velocity for litle particle
	particles[1].v.z = 0;															// (kartesian coordinates)

	double r_module = (particles[0].r - particles[1].r).module();
	particles[1].E = 0.5 * particles[1].m * particles[1].v.module_2() - 0.5 * (PPm::G * particles[0].m * particles[1].m) / (sqrt(r_module * r_module + PPm::r_c));
	particles[1].P = particles[1].v * particles[1].m;
	particles[1].M = particles[1].r * particles[1].P;
	sistem_E[0] += particles[1].E;
	sistem_P[0] += particles[1].P.module();
	sistem_M[0] += particles[1].M.module();
	std::cout << std::setprecision(15) << " E= " << sistem_E[0] << " P= " << sistem_P[0] << " M= " << sistem_M[0] << std::endl;
	//==========================================================

// integration of differential equations
// 
//========================================================================================================================
	std::ofstream positions("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация_измерений\\positions.txt");
	positions << PPm::N << std::endl;
	positions << 0 << std::endl;
	for (int i = 0; i < particles.size(); i++)
	{
		positions << std::setprecision(15) << particles[i].r.x << " " << particles[i].r.y << " " << particles[i].r.z << std::endl;
	}

	int k = 1;
	double dphi = PPm::PI / 200;
	double _phi = PPm::PI / 200;
	for (double t = PPm::t_0 + PPm::dt; t <= PPm::t_1; t += PPm::dt)
	{
		sistem_E.push_back(0);
		sistem_P.push_back(0);
		sistem_M.push_back(0);
		sistem_t.push_back(t);

		double alpha = PPm::G * particles[0].m * particles[1].m;
		double eccentricity = sqrt(1 + (particles[1].E * particles[1].M.module_2()) / (particles[1].m * alpha * alpha));
		double p = particles[1].M.module_2() / (particles[1].m * alpha);

		double _r = p / (1 + eccentricity * cos(_phi));

		particles[1].r.x = _r * cos(_phi);
		particles[1].r.y = _r * cos(_phi);
		particles[1].r.z = 0;
		_phi += dphi;

		for (size_t i = 0; i < particles.size(); i++)
		{
			particles[i].E = 0;
			particles[i].P = vec(0, 0, 0);
			particles[i].M = vec(0, 0, 0);
		}

		/*std::vector<vec> u_i(particles.size(), vec(0, 0, 0));
		for (size_t i = 1; i < particles.size(); i++)
		{
			u_i[i] = particles[i].v;
			particles[i].v = particles[i].v + particles[i].F * PPm::dt;
			particles[i].r = particles[i].r + (particles[i].v + u_i[i]) * PPm::dt * 0.5;
		}


		for (int i = 1; i < particles.size(); i++)
		{
			particles[i].F = vec(0, 0, 0);
		}

		for (size_t i = 1; i < particles.size(); i++)
		{
			for (size_t j = 0; j < particles.size(); j++)
			{
				if (i != j)
				{
					particles[i].F = particles[i].F + PPm::F(particles[i], particles[j]);
				}
			}
		}

		for (size_t i = 1; i < particles.size(); i++)
		{
			particles[i].v = (particles[i].v + u_i[i]) * 0.5 + particles[i].F * 0.5 * PPm::dt;
		}*/

		double r = (particles[0].r - particles[1].r).module();
		particles[1].E = 0.5 * particles[1].m * particles[1].v.module_2() - 0.5 * (PPm::G * particles[0].m * particles[1].m) / (sqrt(r * r + PPm::r_c));
		particles[1].P = particles[1].v * particles[1].m;
		particles[1].M = particles[1].r * particles[1].P;
		sistem_E[k] += particles[1].E;
		sistem_P[k] = sistem_P[k] + particles[1].P.module();
		sistem_M[k] = sistem_M[k] + particles[1].M.module();


		std::cout << k;
		std::cout << std::setprecision(15) << " E= " << sistem_E[k] << " P= " << sistem_P[k] << " M= " << sistem_M[k] << std::endl;
		positions << t << std::endl;
		for (int i = 0; i < particles.size(); i++)
		{
			positions << std::setprecision(15) << particles[i].r.x << " " << particles[i].r.y << " " << particles[i].r.z << std::endl;
		}
		k++;

	}
	std::ofstream conversation_laws("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация_измерений\\measurements.txt");
	for (int i = 0; i < sistem_E.size(); i++)
	{
		conversation_laws << std::setprecision(15) << sistem_t[i] << " " << sistem_E[i] << " " << sistem_P[i] << " " << sistem_M[i] << std::endl;
	}
	//====================================================================================================
}