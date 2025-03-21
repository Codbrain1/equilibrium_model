#include<fstream>
#include<vector>
#include<random>
#include<iostream>
#include"Particle.h"

double random(double beg, double end)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(beg,end);
	return dist(gen);
}
int main()
{
	namespace PPm = Partcile_Particle_model;
	double t_0 = 0, t_1 = 100;
	double dt = 1;
	unsigned N = 200000;
	std::vector<PPm::Particle> PP_model(N); //array particals in model
	//TODO: add initial conditions for Particals
	 
	 
	
	//setting the initial conditions
	// 
	//====================================================================================================
	

	
	// integration of differential equations
	// 
	//====================================================================================================
	std::vector<PPm::Particle> sistem_parametrs(100, PPm::Particle());
	for (int t = 0; t < 100; t++)
	{
		for (size_t i = 0; i < PP_model.size(); i++)
		{
			for (size_t j = 0; j < PP_model.size(); j++)
			{
				if (i != j)
				{
					// TODO: realize impuls of moment
					PP_model[i].v = F(PP_model[i], PP_model[j]) * dt + PP_model[i].v;
					PP_model[i].r = PP_model[i].r + PP_model[i].v * dt;
					PP_model[i].E += PP_model[i].m * PP_model[i].v.module_2() / 2.0 - PPm::G * PP_model[i].m * PP_model[j].m / (PP_model[j].r - PP_model[i].r).module();
					PP_model[i].P = PP_model[i].P + PP_model[i].m * PP_model[i].v;
				}
			}
			sistem_parametrs[t].E += PP_model[i].E;
			sistem_parametrs[t].P = sistem_parametrs[t].P + PP_model[i].P;
			//sistem_parametrs[t].M = sistem_parametrs[t].M + PP_model[i].M;
		}
	}
	for (int t=0;t<100;t++)
	{
		std::cout << "t" << t << ": E = " << sistem_parametrs[t].E << " P = " <<sistem_parametrs[t].P.module()<< std::endl;
	}
	//====================================================================================================

}