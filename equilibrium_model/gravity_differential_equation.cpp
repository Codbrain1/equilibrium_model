#include"Constants.h"
#include"Particle.h"
namespace PP_gravity {//Partical-Partical gravity interaction

	vec F(const Particle& particle, const Particle& particle1)
	{
		vec r_ij = particle1.r - particle.r;
		return G * (particle1.m) / pow((r_ij.module_2() + particle1.r.module_2()),3.0/2.0) * r_ij;
	}

}