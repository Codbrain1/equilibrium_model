#pragma once
namespace Particle_Particle_model {
	//Constants:
	//======================================================================
	//const long double G = 6.67 * 10e-11;
	const long double G = 1;
	const long double r_c = 0.05*0.0; // в данном случае используется r_c*r_c
	const long double PI = 3.141592653589793;
	//initial conditions:
	//======================================================================
	const int N = 3; //число частиц
	const long double t_0 = 0, t_1 = 1000; //начальное и конеченое время
	inline int div = 1; //отвечает за частоту записи данных в файл 1=каждая итерация
	inline long double dt = 0.00001; //шаг времени
	inline long double mindt = 1e-15;
	inline long double maxdt = 0.01;
	inline long double R_max = 1;//максимальный радиус
	inline long double M = 1;//общая масса частиц
	inline long double r_0 = 0; //начальный радиус
	inline long double dr = R_max; //шаг по радиальной компоненте
	const long double r_alpha =0.25;
	const long double sigma_0=0;

}