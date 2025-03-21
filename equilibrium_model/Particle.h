#pragma once
#include<math.h>
#include"Constants.h"
#include"vec.h"

namespace Partcile_Particle_model {

	struct Particle
	{
	public:

		double m, E;
		vec r, v, F, a, M, P;
		// TODO: add function
		Particle() :r(0, 0, 0), v(0, 0, 0), F(0, 0, 0), a(0, 0, 0), M(0, 0, 0), P(0, 0, 0)
		{
			m = 0;
			E = 0;
		}

	};

	vec F(const Particle& particle, const Particle& particle1)
	{
		vec r_ij = particle1.r - particle.r;
		return G * (particle1.m) / pow((r_ij.module_2() + r_c), 3.0 / 2.0) * r_ij;
	}
}