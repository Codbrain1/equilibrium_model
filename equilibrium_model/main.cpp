#include<fstream>
#include<vector>
#include"gravity_differential_equation.cpp"
int main()
{
	namespace ppg = PP_gravity;
	double t_0 = 0, t_1 = 100;
	double dt = 1;
	unsigned N = 200000;
	std::vector<Particle> PP_model(N); //array particals in model
	//TODO: add initial conditions for Particals

	std::vector<Particle> sistem_parametrs(100, Particle());
	for (int t = t_0; t <= t_1; t += dt)
	{
		for (size_t i = 0; i < PP_model.size(); i++)
		{
			for (size_t j = 0; j < PP_model.size(); j++)
			{
				if (i != j)
				{
					//PP_model[i].F = PP_model[i].F + ppg::F(PP_model[i], PP_model[j]); 
					// TODO: realize law of conservation of moment
					PP_model[i].v = ppg::F(PP_model[i], PP_model[j]) * dt + PP_model[i].v;
					PP_model[i].r = PP_model[i].r + PP_model[i].v * dt;
					PP_model[i].E += +PP_model[i].m * PP_model[i].v.module_2() / 2.0 - ppg::G * PP_model[i].m * PP_model[j].m / (PP_model[j].r - PP_model[i].r).module();
					PP_model[i].P = PP_model[i].P + PP_model[i].m * PP_model[i].v;
				}
			}
			sistem_parametrs[t].E += PP_model[i].E;
			sistem_parametrs[t].P = sistem_parametrs[t].P + PP_model[i].P;
			//sistem_parametrs[t].M = sistem_parametrs[t].M + PP_model[i].M;
		}
	}

}