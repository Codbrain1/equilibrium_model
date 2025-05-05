#pragma once
namespace Particle_Particle_model {
	//Constants:
	//======================================================================
	//const double G = 6.67 * 10e-11;
	const double G = 1;
	const double r_c = 0.05*0.05; // в данном случае используется r_c*r_c
	const double PI = 3.141592653589793;
	//initial conditions:
	//======================================================================
	const int N = 512; //число частиц
	const double t_0 = 0, t_1 = 10; //начальное и конеченое время
	int div = 1; //отвечает за частоту записи данных в файл 1=каждая итерация
	double dt = 0.001; //шаг времени
	double mindt = 1e-10;
	double maxdt = 0.01;
	double R_max = 1;//максимальный радиус
	double M = 1;//общая масса частиц
	double r_0 = 0; //начальбный радиус
	double dr = R_max; //шаг по радиальной компоненте
	const double r_alpha =0.25;
	const double sigma_0=0;

}