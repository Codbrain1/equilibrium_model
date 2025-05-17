#pragma once
namespace Particle_Particle_model {
	//Constants:
	//======================================================================
	//const long double G = 6.67 * 10e-11;
	const long double G = 1.0L;
	const long double r_c = 0.05L*0.05L; // в данном случае используется r_c*r_c
	const long double PI = 3.1415926535897932384626433832795L;
	//initial conditions:
	//======================================================================
	const int N = 512; //число частиц
	const long double t_0 = 0.0L, t_1 = 1000.0L; //начальное и конеченое врем
	inline int div = 1; //отвечает за частоту записи данных в файл 1=каждая итерация
	inline long double dt = 0.001L; //шаг времени
	inline long double mindt = 1e-15L;
	inline long double maxdt = 0.01L;
	inline long double R_max = 1L;//максимальный радиус
	inline long double M = 1.0L;//общая масса частиц
	inline long double r_0 = 0.0L; //начальный радиус
	inline long double dr = R_max; //шаг по радиальной компоненте
	const long double r_alpha =0.25L;
	const long double sigma_0=0.0L;

}