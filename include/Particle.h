#pragma once
#include"Constants.h"
#include"vec.h"
namespace Particle_Particle_model {
	
	struct  Particle
	{
	public:

		long double m, E;
		vec r, v, F, L, P, a;
		// TODO: add function
		Particle():
			r(0, 0, 0), 
			v(0, 0, 0), 
			F(0, 0, 0), 
			L(0, 0, 0), 
			P(0, 0, 0),
			a(0, 0, 0),
			m(M/N),
			E(0)
		{}
		void Evaluate_a() { a = F / m; }

	};

	inline vec F(const Particle& p1, const Particle& p2)
	{
		const vec r_ij = p2.r - p1.r;
		const long double r =sqrt(r_ij.module_2()+r_c);
		return  r_ij*(G * (p2.m*p1.m) / (r*r*r));
	}

	inline vec calc_f(const Particle& p1, const Particle& p2)
	{
		const vec r_ij = p2.r - p1.r;
		const long double r = sqrt(r_ij.module_2() + r_c);
		return  r_ij * (G * (p2.m) / (r * r * r));
	}	
}