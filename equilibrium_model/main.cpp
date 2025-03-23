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

std::vector<double> sistem_E((int)(PPm::t_1 - PPm::t_0) / PPm::dt + 1, 0);
std::vector<double> sistem_P((int)(PPm::t_1 - PPm::t_0) / PPm::dt + 1, 0);
std::vector<double> sistem_M((int)(PPm::t_1 - PPm::t_0) / PPm::dt + 1, 0);

void calculating(std::vector<PPm::Particle>& particles, int k, int i_0, int i_1)
{

	for (size_t i = i_0; i < i_1; i++)
	{
		particles[i].r = particles[i].r + particles[i].v * PPm::dt + particles[i].F * 0.5 * PPm::dt * PPm::dt;
	}

	std::vector<vec> F_i(i_1 - i_0, vec(0, 0, 0));
	for (int i = 0; i < F_i.size(); i++)
	{
		F_i[i] = particles[i + i_0].F;
		particles[i+i_0].F = vec(0, 0, 0);
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
		particles[i].v = particles[i].v + (F_i[i - i_0] + particles[i].F) * 0.5 * PPm::dt;
	}

	for (size_t i = i_0; i < i_1; i++)
	{
		for (size_t j = 0; j < particles.size(); j++)
		{
			if (i != j)
			{
				particles[i].E -= PPm::G * particles[i].m * particles[j].m / (particles[j].r - particles[i].r).module();
				
			}
		}
		particles[i].E += particles[i].m * particles[i].v.module_2()*0.5;
		particles[i].P = particles[i].v * particles[i].m;
		particles[i].M = particles[i].r * particles[i].P;
		sistem_E[k] += particles[i].E;
		sistem_P[k] = sistem_P[k] + particles[i].P.module();
		sistem_M[k] = sistem_M[k] + particles[i].M.module();
	}


}
int main()
{


	std::vector<PPm::Particle> particles(PPm::N); //array particals in model
	//TODO: add initial conditions for Particals



//setting the initial conditions
// 
//====================================================================================================

	//setting the initial position
	// 
	//====================================================
	double index = 0;
	double r_k = PPm::r_0 + PPm::dr * 0.5;
	int sum = 0;
	double Summ_sigma_r = 0;
	while (r_k <= PPm::R_max)
	{

		double sigma_r = PPm::sigma(r_k);
		int N_k = round(sigma_r * PPm::PI * r_k * PPm::dr * PPm::N);
		sum += N_k;
		double phi = 0;
		double d_phi = 2 * PPm::PI / N_k;
		double ds = PPm::PI * r_k * PPm::dr;
		Summ_sigma_r += sigma_r * ds;
		for (int i = index; i < index + N_k; i++)
		{
			particles[i].r.x = (r_k + random(-PPm::dr * 0.5, PPm::dr * 0.5)) * cos(phi);
			particles[i].r.y = (r_k + random(-PPm::dr * 0.5, PPm::dr * 0.5)) * sin(phi);
			particles[i].r.z = 0;
			phi += d_phi;
		}
		r_k += PPm::dr;
		index += N_k;
	}
	std::cout << "\n" << sum << " " << Summ_sigma_r << std::endl;

	if (sum == PPm::N)
	{
		//setting the initial velocity
		// 
		//==========================================================
		for (size_t i = 0; i < particles.size(); i++)
		{
			for (size_t j = 0; j < particles.size(); j++)
			{

				if (i != j)
				{
					// TODO: realize impuls of moment
					particles[i].F = F(particles[i], particles[j]) + particles[i].F;
					particles[i].E -= PPm::G * particles[i].m * particles[j].m / (particles[j].r - particles[i].r).module();
				}
			}
		}
		for (size_t i = 0; i < particles.size(); i++)
		{
			double v_radial = sqrt(particles[i].r.module() * particles[i].F.module());
			particles[i].v.x = v_radial*cos(atan(particles[i].F.y / particles[i].F.x));
			particles[i].v.y = v_radial*sin(atan(particles[i].F.y / particles[i].F.x));
			particles[i].v.z = 0;
			particles[i].E += particles[i].m * particles[i].v.module_2() / 2.0;
			particles[i].P = particles[i].v * particles[i].m;
			particles[i].M = particles[i].r * particles[i].P;
			sistem_E[0] += particles[i].E;
			sistem_P[0] += particles[i].P.module();
			sistem_M[0] += particles[i].M.module();
		}
		std::cout << std::setprecision(20) << " E= " << sistem_E[0] << " P= " << sistem_P[0] << " M= " << sistem_M[0] << std::endl;

// integration of differential equations
// 
//========================================================================================================================
		int k = 1;
		for (double t = PPm::t_0; t <= PPm::t_1; t += PPm::dt)
		{
			for (size_t i = 0; i < particles.size(); i++)
			{
				particles[i].E = 0;
				particles[i].P = vec(0, 0, 0);
				particles[i].M = vec(0, 0, 0);
			}

			std::thread th1(calculating, std::ref(particles), k, 0, particles.size() / 4);
			std::thread th2(calculating, std::ref(particles), k, particles.size() / 4, particles.size() / 2);
			std::thread th3(calculating, std::ref(particles), k, particles.size() / 2, 3 * particles.size() / 4);
			std::thread th4(calculating, std::ref(particles), k, 3 * particles.size() / 4, particles.size());
			th1.join();
			th2.join();
			th3.join();
			th4.join();
			
			std::cout << k;
			std::cout << std::setprecision(10) << " E= " << sistem_E[k] << " P= " << sistem_P[k] << " M= " << sistem_M[k] << std::endl;
			k++;
		}
		std::ofstream conversation_laws("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация измерений\\measurements.txt");
		for (int i = 0; i < sistem_E.size(); i++)
		{
			conversation_laws << std::setprecision(20) << i << " " << sistem_E[i] << " " << sistem_P[i] << " " << sistem_M[i] << std::endl;
		}
//====================================================================================================
	}
}