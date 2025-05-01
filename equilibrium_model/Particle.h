#pragma once
#include"Constants.h"
#include"vec.h"
namespace Partcile_Particle_model {
	
	struct Particle
	{
	public:

		double m, E;
		vec r, v, F, M, P;
		// TODO: add function
		Particle() :r(0, 0, 0), v(0, 0, 0), F(0, 0, 0), M(0, 0, 0), P(0, 0, 0)
		{
			m = Partcile_Particle_model::M/N;
			E = 0;

		}

	};

	vec F(const Particle& particle, const Particle& particle1)
	{
		vec r_ij = particle1.r - particle.r;
		double r =sqrt(r_ij.module_2()+r_c);
		return  r_ij*(G * (particle1.m) / (r*r*r));
	}
	/// <summary>
	/// функция вычисляющая поверхностную плотность в тонком диске, определена в полярных координатах
	/// </summary>
	/// <param name="r"></param>
	/// <returns></returns>
	double sigma(double r) { return sigma_0 *r_alpha* r_alpha /r; }

	//double sigma_exp(double r) { return sigma_0 * exp(-r / r_alpha); }
}