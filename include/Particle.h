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
			r(0.0L, 0.0L, 0.0L), 
			v(0.0L, 0.0L, 0.0L), 
			F(0.0L, 0.0L, 0.0L), 
			L(0.0L, 0.0L, 0.0L), 
			P(0.0L, 0.0L, 0.0L),
			a(0.0L, 0.0L, 0.0L),
			m(M/N),
			E(0.0L)
		{}
		void Evaluate_a() { a = F / m; }

	};

	inline vec F(const Particle& p1, const Particle& p2)
	{
		const vec r_ij = p2.r - p1.r;
		const long double r =sqrtl(r_ij.module_2()+r_c);
		return  r_ij*(G * (p2.m*p1.m) / (r*r*r));
	}

	inline vec calc_f(const Particle& p1, const Particle& p2)
	{
		const vec r_ij = p2.r - p1.r;
		const long double r = sqrtl(r_ij.module_2() + r_c);
		return  r_ij * (G * (p2.m) / (r * r * r));
	}	
}