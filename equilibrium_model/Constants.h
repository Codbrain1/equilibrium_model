#pragma once
namespace Partcile_Particle_model {
	//Constants:
	//======================================================================
	//const double G = 6.67 * 10e-11;
	const double G = 1;
	const double r_c = 0.005*0.005; // в данном случае используется r_c*r_c
	const double PI = 3.141592653589793238462643383279502884;
	//initial conditions:
	//======================================================================
	const double N = 2; //число частиц
	const double t_0 = 0, t_1 = 10; //начальное и конеченое время
	const int div = 1; //отвечает за частоту записи данных в файл 1=каждая итерация
	double dt = 0.1; //шаг времени
	double R_max = 1;//максимальный радиус
	double M = 1;//общая масса частиц
	double r_0 = 0; //начальбный радиус
	double dr = R_max; //шаг по радиальной компоненте
	const double r_alpha =R_max/sqrt(N);
	const double sigma_0=N/PI;

}