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
void calculating(std::vector<PPm::Particle>& particles, int i, int j_0, int j_1)
{
	for (size_t j = j_0; j <j_1; j++)
	{
		if (i != j)
		{
			// TODO: realize impuls of moment
			particles[i].v = F(particles[i], particles[j]) * PPm::dt + particles[i].v;
			particles[i].r = particles[i].r + particles[i].v * PPm::dt;
			particles[i].E += particles[i].m * particles[i].v.module_2() / 2.0 - PPm::G * particles[i].m * particles[j].m / (particles[j].r - particles[i].r).module();
			particles[i].P = particles[i].P + particles[i].v * particles[i].m;
			particles[i].M = particles[i].M + particles[i].r * particles[i].P;
		}
	}
}
int main()
{

	
	std::vector<PPm::Particle> particles(PPm::N); //array particals in model
	//TODO: add initial conditions for Particals



	//setting the initial conditions
	// 
	//====================================================================================================
	double index = 0;
	double r_k = PPm::r_0 + PPm::dr / 2.0;
	int sum = 0;
	while (r_k <= PPm::R_max)
	{
		double sigma_r = PPm::sigma(r_k);
		int N_k = round(sigma_r * PPm::PI * r_k * PPm::dr * PPm::N);
		sum += N_k;
		double phi = 0;
		double d_phi = 2 * PPm::PI / N_k;
		for (int i = index; i < index + N_k; i++)
		{
			particles[i].r.x = (r_k+random(r_k-PPm::dr,r_k+PPm::dr)) * cos(phi);
			particles[i].r.y = (r_k + random(r_k - PPm::dr, r_k + PPm::dr)) * sin(phi);
			particles[i].r.z = 0;
			phi += d_phi;
		}
		r_k += PPm::dr;
		index += N_k;
	}
	std::cout << "\n" << sum << std::endl;

	for (size_t i = 0; i < particles.size(); i++)
	{
		vec F_ij;
		for (size_t j = 0; j < particles.size(); j++)
		{
			
			if (i != j)
			{
				// TODO: realize impuls of moment
				F_ij = F(particles[i], particles[j]) + F_ij;
				particles[i].E -= PPm::G * particles[i].m * particles[j].m / (particles[j].r - particles[i].r).module();
			}
		}
		double v_radial = sqrt(particles[i].r.module() * F_ij.module());
		particles[i].v.x = v_radial * cos(atan(F_ij.y / F_ij.x));
		particles[i].v.y = v_radial * sin(atan(F_ij.y / F_ij.x));
		particles[i].v.z = 0;
		particles[i].E += particles[i].m * particles[i].v.module_2() / 2.0;
		particles[i].P = particles[i].v * particles[i].m;
		particles[i].M = particles[i].r * particles[i].P;
	}

	// integration of differential equations
	// 
	//====================================================================================================
	int k = 0;
	std::vector<double> sistem_E((int)(PPm::t_1 - PPm::t_0) / PPm::dt, 0);
	std::vector<double> sistem_P((int)(PPm::t_1 - PPm::t_0) / PPm::dt, 0);
	std::vector<double> sistem_M((int)(PPm::t_1 - PPm::t_0) / PPm::dt, 0);

	for (double t = PPm::t_0; t < PPm::t_1; t += PPm::dt)
	{
		for (size_t i = 0; i < particles.size(); i++)
		{
			
			/*std::thread th1(calculating,particles,i,0, particles.size()/4);
			std::thread th2(calculating, particles, i, particles.size() / 4, particles.size() / 2);
			std::thread th3(calculating, particles, i, particles.size() / 2, particles.size() / 4+ particles.size() / 2);
			std::thread th4(calculating, particles, i, particles.size() / 4 + particles.size() / 2, particles.size());
			th1.join();
			th2.join();
			th3.join();
			th4.join();*/
			sistem_E[k] += particles[i].E;
			sistem_P[k] = sistem_P[k] + particles[i].P.module();
			sistem_M[k] = sistem_M[k] + particles[i].M.module();

		}
		for (size_t i = 0; i < particles.size(); i++)
		{
			particles[i].E = 0;
			particles[i].P = vec(0, 0, 0);
			particles[i].M = vec(0, 0, 0);
		}
		std::cout << k;
		std::cout << std::setprecision(10) << " E= " << sistem_E[k] <<" P= "<< sistem_P[k] <<" M= "<< sistem_M[k] << std::endl;
		k++;
	}
	std::ofstream conversation_laws("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация измерений\\measurements.txt");
	for (int i = 0; i < sistem_E.size(); i++)
	{
		conversation_laws<< std::setprecision(10) << i <<" "<< sistem_E[i] << " " << sistem_P[i] << " " << sistem_M[i] << std::endl;
	}
	//====================================================================================================

}