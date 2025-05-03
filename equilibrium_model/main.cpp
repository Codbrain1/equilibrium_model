#include<fstream>
#include<vector>
#include<random>
#include<iostream>
#include<cmath>
#include"Particle.h"
#include <iomanip>
#include<chrono>
#include<omp.h>

double random(double beg, double end)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(beg, end);
	return dist(gen);
}

namespace PPm = Particle_Particle_model;

double sistem_E;
double sistem_P;
double sistem_M;
double sistem_t;
double sistem_E_k;
double sistem_E_p;

void calc_forces(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
	for (int i = i_0; i < i_1; i++)
	{
		ps[i].F = vec(0, 0, 0);
	}

	for (int i = i_0; i < i_1; i++)
	{
		for (int j = i + 1; j < ps.size(); j++)
		{
			if (i != j)
			{
				auto F_ij = PPm::F(ps[i], ps[j]);
				ps[i].F = ps[i].F + F_ij;
				ps[j].F = ps[j].F - F_ij;
			}
		}
	}
}

void calc_acceleration_asinc(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
#pragma omp parallel for
	for (int i = i_0; i < i_1; i++)
	{
		ps[i].a = vec(0, 0, 0);

		
		for (int j = 0; j < ps.size(); j++)
		{
			if (i != j)
			{
				auto f_ij = PPm::calc_f(ps[i], ps[j]);
				ps[i].a = ps[i].a + f_ij;
			}
		}
	}
}

void KDK_parallel(std::vector<PPm::Particle>& ps) {
	// Проверка на корректность количества потоков

	const int N = ps.size();
	std::vector<vec> u_i(N);  // Сохраняем начальные скорости для всех частиц

	// Функция для обновления координат и скоростей (Drift + Kick)
	#pragma omp parallel for 
	for (int i = 0; i < N; ++i) {
		u_i[i] = ps[i].v + ps[i].a * PPm::dt;                     // Kick (обновление скорости)
		ps[i].r = ps[i].r + (ps[i].v + u_i[i]) * 0.5 * PPm::dt;  // Drift (обновление координат)
	}

	calc_acceleration_asinc(ps, 0, ps.size());

	// Финальное обновление скоростей (Kick)
	#pragma omp parallel for 
	for (int i = 0; i < N; i++) {
		ps[i].v = (ps[i].v + u_i[i]) * 0.5 + ps[i].a * PPm::dt * 0.5;
	}
}

void KDK(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
	std::vector<vec> u_i(ps.size(), vec(0, 0, 0));
	for (int i = i_0; i < i_1; i++)
	{
		u_i[i - i_0] = ps[i].v;
		ps[i].v = ps[i].v + ps[i].a * PPm::dt;
		ps[i].r = ps[i].r + (ps[i].v + u_i[i - i_0]) * PPm::dt * 0.5;
	}
	calc_forces(ps, 0, ps.size());
	for (auto& p : ps)
		p.Evaluate_a();

	for (int i = 0; i < ps.size(); i++)
	{
		ps[i].v = (ps[i].v + u_i[i]) * 0.5 + ps[i].a * 0.5 * PPm::dt;
	}
}

void RungeKutta4(std::vector<PPm::Particle>& particles, int i_0, int i_1) {
	std::vector<vec> k1(PPm::N, vec(0, 0, 0)), k2(PPm::N, vec(0, 0, 0)), k3(PPm::N, vec(0, 0, 0)), k4(PPm::N, vec(0, 0, 0));
	std::vector<vec> k1v(PPm::N, vec(0, 0, 0)), k2v(PPm::N, vec(0, 0, 0)), k3v(PPm::N, vec(0, 0, 0)), k4v(PPm::N, vec(0, 0, 0));
	// Рунге-Кутта 4 порядка
	// Шаг 1
	calc_forces(particles, i_0, i_1);
	for (int i = i_0; i < i_1; i++)
		particles[i].Evaluate_a();

	for (size_t i = i_0; i < i_1; ++i) {
		k1v[i] = particles[i].a * PPm::dt;
		k1[i] = particles[i].v * PPm::dt;
	}

	// Шаг 2
	for (size_t i = i_0; i < i_1; i++) {
		k2[i] = particles[i].r + k1[i] * 0.5;
		k2v[i] = particles[i].v + k1v[i] * 0.5;
	}

	calc_forces(particles, i_0, i_1);
	for (int i = i_0; i < i_1; i++)
		particles[i].Evaluate_a();

	for (size_t i = 0; i < i_1; i++) {
		k2v[i] = particles[i].a * PPm::dt;
		k2[i] = (particles[i].v + k1[i] * 0.5) * PPm::dt;
	}

	// Шаг 3
	for (size_t i = 0; i < i_1; ++i) {
		k3[i] = particles[i].r + k2[i] * 0.5;
		k3v[i] = particles[i].v + k2[i] * 0.5;
	}
	calc_forces(particles, i_0, i_1);
	for (int i = i_0; i < i_1; i++)
		particles[i].Evaluate_a();

	for (size_t i = i_0; i < PPm::N; ++i) {
		k3v[i] = particles[i].a * PPm::dt;
		k3[i] = (particles[i].v + k2v[i] * 0.5) * PPm::dt;
	}

	// Шаг 4
	for (size_t i = i_0; i < i_1; i++) {
		k4[i] = particles[i].r + k3[i];
		k4v[i] = particles[i].v + k3v[i];
	}
	calc_forces(particles, i_0, i_1);
	for (int i = i_0; i < i_1; i++)
		particles[i].Evaluate_a();

	for (size_t i = i_0; i < i_1; i++) {
		k4v[i] = particles[i].a * PPm::dt;
		k4[i] = (particles[i].v + k3v[i]) * PPm::dt;
	}


	for (size_t i = i_0; i < i_1; i++) {
		particles[i].r = particles[i].r + (k1[i] + k2[i] * 2 + k3[i] * 2 + k4[i]) / 6.0;
		particles[i].v = particles[i].v + (k1v[i] + k2v[i] * 2 + k3v[i] * 2 + k4v[i]) / 6.0;
	}
}

void calculate_conversation_laws_parallel(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
	#pragma omp parallel for
	for (size_t i = i_0; i < i_1; i++)
	{
		ps[i].E = 0;
		ps[i].P = vec(0, 0, 0);
		ps[i].L = vec(0, 0, 0);
	}

	double local_E_k = 0, local_E_p = 0, local_E = 0, local_P = 0, local_M = 0;

	#pragma omp parallel for reduction(+:local_E_k, local_E_p, local_E, local_P, local_M)
	for (int i = i_0; i < i_1; i++)
	{
		auto E_k = 0.5 * ps[i].m * ps[i].v.module_2();
		ps[i].E = E_k;	// set kinetical energy
		ps[i].P = ps[i].v * ps[i].m;					// set impulse
		ps[i].L = ps[i].r * ps[i].P;					// set moment impulse

		local_E_k += E_k;

		for (int j = i + 1; j < ps.size(); j++)
		{
			vec r_ij = ps[j].r - ps[i].r;
			double r = sqrt(r_ij.module_2() + PPm::r_c);
			auto E_p = -PPm::G * ps[i].m * ps[j].m / r;//set potential energy
			ps[i].E += E_p;
			local_E_p += E_p;
		}
		local_E += ps[i].E;
		local_P += ps[i].P.module();
		local_M += ps[i].L.module();
	}

	#pragma omp critical
	{
		sistem_E_k += local_E_k;
		sistem_E_p += local_E_p;
		sistem_E += local_E;
		sistem_P += local_P;
		sistem_M += local_M;
	}
}

void calculate_conversation_laws(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
	for (size_t i = i_0; i < i_1; i++)
	{
		ps[i].E = 0;
		ps[i].P = vec(0, 0, 0);
		ps[i].L = vec(0, 0, 0);
	}

	for (int i = i_0; i < i_1; i++)
	{
		auto E_k = 0.5 * ps[i].m * ps[i].v.module_2();
		ps[i].E = E_k;	// set kinetical energy
		ps[i].P = ps[i].v * ps[i].m;					// set impulse
		ps[i].L = ps[i].r * ps[i].P;					// set moment impulse

		sistem_E_k += E_k;

		for (int j = i + 1; j < ps.size(); j++)
		{
			if (i != j)
			{
				vec r_ij = ps[j].r - ps[i].r;
				double r = sqrt(r_ij.module_2() + PPm::r_c);
				auto E_p = -PPm::G * ps[i].m * ps[j].m / r;//set potential energy
				ps[i].E += E_p;
				sistem_E_p += E_p;
			}
		}
		sistem_E += ps[i].E;
		sistem_P += ps[i].P.module();
		sistem_M += ps[i].L.module();
	}
}


void set_dinamic_step(std::vector<PPm::Particle>& ps)
{
	double min_dist = (ps[0].r - ps[1].r).module();
	double max_a = ps[0].a.module();
	for (int i = 0; i < ps.size(); i++)
	{
		for (int j = i + 1; j < ps.size(); j++)
		{
			double dist = (ps[i].r - ps[j].r).module();
			if (dist < min_dist)
				min_dist = dist;
		}
		double a = ps[i].a.module();
		if (max_a < a)
			max_a = a;
	}
	double new_dt = (std::min)(0.01 * min_dist, 0.1 / sqrt(max_a));
	new_dt = (std::max)(0.5 * PPm::dt, (std::min)(2.0 * PPm::dt, new_dt));
	const double min_dt = PPm::mindt;
	const double max_dt = PPm::maxdt;
	PPm::dt = (std::max)(min_dt, (std::min)(new_dt, max_dt));
}

void set_dinamic_step_parallel(std::vector<PPm::Particle>& ps)
{
	double min_dist = (ps[0].r - ps[1].r).module();
	double max_a = ps[0].a.module();
	#pragma omp parallel for reduction(min: min_dist) reduction(max: max_a)
	for (int i = 0; i < ps.size(); i++)
	{
		for (int j = i + 1; j < ps.size(); j++)
		{
			double dist = (ps[i].r - ps[j].r).module();
			if (dist < min_dist)
				min_dist = dist;
		}
		double a = ps[i].a.module();	
		if (max_a < a)
			max_a = a;
	}
	double new_dt = (std::min)(0.01 * min_dist, 0.1 / sqrt(max_a));
	new_dt = (std::max)(0.5 * PPm::dt, (std::min)(2.0 * PPm::dt, new_dt));
	const double min_dt = PPm::mindt;
	const double max_dt = PPm::maxdt;
	PPm::dt = (std::max)(min_dt, (std::min)(new_dt, max_dt));
}

void set_exponential_disk(std::vector<PPm::Particle>& ps)
{
	using namespace Particle_Particle_model;
	//setting the initial position
		// 
		//====================================================
	std::mt19937 gen(42);
	std::uniform_real_distribution<double> angle_dist(0, 2 * PI);

	const double r_core = 0.05 * r_alpha;//ядро смягчения потенциала

	// Экспоненциальное распределение с ядром
	for (int i = 0; i < N; ++i) {
		double u = (i + 0.5) / N;
		double r = -r_alpha * log(1 - u);

		// Плавный переход в центре
		double r_eff = sqrt(r * r + r_core * r_core);
		double theta = angle_dist(gen);

		// Кривая вращения с поправкой на ядро
		double v_circ = sqrt(M * pow(r_eff, 3) / pow(r_eff + r_alpha, 3));

		// Дисперсия скоростей (20% радиальная, 10% тангенциальная)
		std::normal_distribution<double> rad_vel(0, 0.2 * v_circ);
		std::normal_distribution<double> tang_vel(v_circ, 0.1 * v_circ);

		ps[i].r.x = r * cos(theta);
		ps[i].r.y = r * sin(theta);
		ps[i].r.z = 0;

		ps[i].v.x = -tang_vel(gen) * sin(theta) + rad_vel(gen) * cos(theta);
		ps[i].v.y = tang_vel(gen) * cos(theta) + rad_vel(gen) * sin(theta);
		ps[i].v.z = 0;

		ps[i].m = M / N;
	}
	calc_forces(ps, 0, ps.size());
	for (auto& i : ps)
		i.Evaluate_a();
}

void set_initial_circle(std::vector<PPm::Particle>& ps)
{
	using namespace Particle_Particle_model;
	//setting the initial position
		// 
		//====================================================
	for (int i = 0; i < N; i++)
	{
		ps[i].r.x = R_max * cos(2 * PI / (N)*i);
		ps[i].r.y = R_max * sin(2 * PI / (N)*i);
		ps[i].r.z = 0;

	}

	calc_forces(ps, 0, ps.size());
	for (auto& i : ps)
		i.Evaluate_a();
	//setting the initial velocity
	// 
	//==========================================================
	for (auto& i : ps) {

		double v_asimutal = sqrt(i.r.module() * i.a.module());		// (cylindrical coordinates)
		double phi = atan2(i.r.y, i.r.x);
		double r = i.r.module();

		i.v.x = -r * v_asimutal * sin(phi);		// (kartesian coordinates)
		i.v.y = r * v_asimutal * cos(phi);
		i.v.z = 0;
	}
}

int main()
{
	// omp_set_num_threads(std::thread::hardware_concurrency());
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<PPm::Particle> particles(PPm::N); //array particals in model

	//setting the initial conditions
	// 
	//==========================================================
	
	//set_exponential_disk(particles);
	set_exponential_disk(particles);
	std::cout << "set initial conditions\n";

	sistem_E = 0;
	sistem_P = 0;
	sistem_M = 0;
	sistem_t = 0;
	sistem_E_k = 0;
	sistem_E_p = 0;

	calculate_conversation_laws(particles, 0, particles.size());
	std::cout << "set intitial conversation laws\n";
	//==========================================================

// integration of differential equations
// 
//========================================================================================================================

	//write initial positions into file
	std::ofstream positions("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация_измерений\\positions.txt");
	positions << PPm::N << std::endl;
	positions << PPm::t_0 << std::endl;
	for (auto& p : particles)
	{
		positions << std::setprecision(13) << p.r.x << " " << p.r.y << " " << p.r.z << " " << p.v.x << " " << p.v.y << " " << p.v.z << std::endl;
	}

	//write initial conversation laws into file
	std::ofstream conversation_laws("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация_измерений\\measurements.txt");
	conversation_laws << std::fixed << std::setprecision(13) << sistem_t << " " << sistem_E << " " << sistem_P << " " << sistem_M << " " << sistem_E_k << " " << sistem_E_p << " " << std::endl;


	//write execution time for initial calculating
	std::ofstream time("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация_измерений\\time.txt");
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	time << elapsed.count()<<" "<<PPm::dt<<std::endl;

	//	===BEGIN SIMULATION===
	int b = 1;
	std::cout << "start calculating\n";
	// HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //current simulate data
	// COORD pos = { 0, 4 };

	for (double t = PPm::t_0 + PPm::dt; t <= PPm::t_1; t += PPm::dt)
	{
		
		//KDK(particles, 0, particles.size());
		KDK_parallel(particles);
		//set_dinamic_step(particles);
		if (b % PPm::div == 0) {
			sistem_E = 0;
			sistem_P = 0;
			sistem_M = 0;
			sistem_t = t;
			sistem_E_k = 0;
			sistem_E_p = 0;
			calculate_conversation_laws_parallel(particles, 0, particles.size());

			positions << t << std::endl;
			for (auto& p : particles)
			{
				positions << std::setprecision(13) << p.r.x << " " << p.r.y << " " << p.r.z << " " << p.v.x << " " << p.v.y << " " << p.v.z << std::endl;
			}
			conversation_laws << std::fixed << std::setprecision(13) << sistem_t << " " << sistem_E << " " << sistem_P << " " << sistem_M << " " << sistem_E_k << " " << sistem_E_p << " " << std::endl;
			
			auto stop = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> el = stop - start;
			time << el.count() << " " << PPm::dt << std::endl;
		}
		b++;
		// SetConsoleCursorPosition(hConsole, pos);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end - start;
		std::cout << "\033[4;0H";
		std::cout << std::setprecision(5) << "[(t_n,t_1): " << t << "/" << PPm::t_1 << ", dt = " << PPm::dt << ", Время выполнения: " << elapsed.count() << " секунд" << "]";
	}
	std::cout << "\nend calculating\n";
	
	//	===END SIMULATION===

	//python visualisation
	//====================================================================================================

// comment	
	return 0;
}
