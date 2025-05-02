#include<fstream>
#include<vector>
#include<random>
#include<iostream>
#include<cmath>
#include"Particle.h"
#include <iomanip>
#include<thread>
#include"wrapper.h"
#include<mutex>
#include<Windows.h>

std::mutex mut;
double random(double beg, double end)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(beg, end);
	return dist(gen);
}

namespace PPm = Partcile_Particle_model;

std::vector<double> sistem_E;
std::vector<double> sistem_P;
std::vector<double> sistem_M;
std::vector<double> sistem_t;
std::vector<double> sistem_E_k;
std::vector<double> sistem_E_p;

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
				mut.lock();
				auto F_ij = PPm::F(ps[i], ps[j]);
				ps[i].F = ps[i].F + F_ij;
				ps[j].F = ps[j].F - F_ij;
				mut.unlock();
			}
		}
	}
}
void calc_forces_asinc(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
	for (int i = i_0; i < i_1; i++)
	{
		for (int j = i + 1; j < ps.size(); j++)
		{
			if (i != j)
			{
				mut.lock();
				auto F_ij = PPm::F(ps[i], ps[j]);
				ps[i].F = ps[i].F + F_ij;
				ps[j].F = ps[j].F - F_ij;
				mut.unlock();
			}
		}
	}
}
#include <vector>
#include <thread>
#include <functional>

int KDK(std::vector<PPm::Particle>& ps, int th_n) {
	// Проверка на корректность количества потоков
	if (th_n <= 0 || th_n > std::thread::hardware_concurrency()) {
		th_n = std::thread::hardware_concurrency();
	}

	const int N = ps.size();
	std::vector<vec> u_i(N);  // Сохраняем начальные скорости для всех частиц
	std::vector<std::thread> threads;

	// Функция для обновления координат и скоростей (Drift + Kick)
	for (int i = 0; i < ps.size(); ++i) {
		u_i[i] = ps[i].v;                          // Сохраняем начальную скорость
		ps[i].v = ps[i].v + ps[i].a * PPm::dt;              // Kick (обновление скорости)
		ps[i].r = ps[i].r + (ps[i].v + u_i[i]) * 0.5 * PPm::dt; // Drift (обновление координат)
	}
	for (auto& p : ps)
	{
		p.F = vec(0, 0, 0);
	}
	// Запуск потоков для обновления частиц
	int chunk_size = N / th_n;
	// Параллельное вычисление сил (Kick)
	for (int i = 0; i < th_n; ++i) {
		int start = i * chunk_size;
		int end = (i == th_n - 1) ? N : start + chunk_size;
		threads.emplace_back(calc_forces_asinc, std::ref(ps), start, end);
	}

	for (auto& t : threads) {
		t.join();
	}
	threads.clear();

	// Обновление ускорений (после вычисления сил)
	for (PPm::Particle& p : ps) {
		p.Evaluate_a();  // Предполагается, что это метод класса Particle
	}

	// Финальное обновление скоростей (Kick)
	for (int i = 0; i < ps.size(); i++) {
		ps[i].v = (ps[i].v + u_i[i]) * 0.5 + ps[i].a * PPm::dt * 0.5;
	}
	return th_n;
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

void Leapfrog(std::vector<PPm::Particle>& particles, int i_0, int i_1)
{
	calc_forces(particles, i_0, i_1);
	for (int i = i_0; i < i_1; i++)
		particles[i].Evaluate_a();

	std::vector<vec> dv(particles.size());
	for (int i = 0; i < particles.size(); i++)
	{
		dv[i] = particles[i].a * PPm::dt;
		particles[i].r = particles[i].r + (particles[i].v + dv[i] * 0.5) * PPm::dt;
		particles[i].v = particles[i].v + dv[i];

	}
}

void calculate_conversation_laws(std::vector<PPm::Particle>& ps, int k, int i_0, int i_1)
{
	for (size_t i = i_0; i < i_1; i++)
	{
		ps[i].E = 0;
		ps[i].P = vec(0, 0, 0);
		ps[i].M = vec(0, 0, 0);
	}

	for (int i = i_0; i < i_1; i++)
	{
		auto E_k = 0.5 * ps[i].m * ps[i].v.module_2();
		ps[i].E = E_k;	// set kinetical energy
		ps[i].P = ps[i].v * ps[i].m;					// set impulse
		ps[i].M = ps[i].r * ps[i].P;					// set moment impulse

		sistem_E_k[k] += E_k;

		for (int j = i + 1; j < ps.size(); j++)
		{
			if (i != j)
			{
				vec r_ij = ps[j].r - ps[i].r;
				double r = sqrt(r_ij.module_2() + PPm::r_c);

				auto E_p = -PPm::G * ps[i].m * ps[j].m / r;//set potential energy

				ps[i].E += E_p;
				sistem_E_p[k] += E_p;
			}
		}
		sistem_E[k] += ps[i].E;
		sistem_P[k] = sistem_P[k] + ps[i].P.module();
		sistem_M[k] = sistem_M[k] + ps[i].M.module();
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

void set_initial_conditions(std::vector<PPm::Particle>& ps)
{
	using namespace Partcile_Particle_model;
	//setting the initial position
		// 
		//====================================================
	double dr = 0.1;
	double r_0 = 0.2;
	for (int j = 1; j * dr + r_0 <= R_max; j++)
	{
		for (int i = (j - 1) * N / 8; i < j * N / 8; i++)
		{
			ps[i].r.x = (r_0 + j * dr * 0.5 + random(-dr * 0.5, dr * 0.5)) * cos((2 * PI * i) / (N / 8) + i * PI/6);
			ps[i].r.y = (r_0 + j * dr * 0.5 + random(-dr * 0.5, dr * 0.5)) * sin((2 * PI * i) / (N / 8) + i * PI/6);
			ps[i].r.z = 0;
			ps[i].m = M / N;
		}
	}

	calc_forces(ps, 0, N);
	for (auto& p : ps)
		p.Evaluate_a();
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
	std::vector<PPm::Particle> particles(PPm::N); //array particals in model
	//setting the initial conditions
	// 
	//==========================================================
	set_initial_conditions(particles);

	std::cout << "set initial conditions\n";

	sistem_E.push_back(0);
	sistem_P.push_back(0);
	sistem_M.push_back(0);
	sistem_t.push_back(0);
	sistem_E_k.push_back(0);
	sistem_E_p.push_back(0);

	calculate_conversation_laws(particles, 0, 0, particles.size());
	std::cout << "set intitial conversation laws\n";
	//==========================================================

// integration of differential equations
// 
//========================================================================================================================

	std::cout << "write into file\n";

	std::ofstream positions("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация_измерений\\positions.txt");
	positions << PPm::N << std::endl;
	positions << PPm::t_0 << std::endl;

	for (auto& p : particles)
	{
		positions << std::setprecision(13) << p.r.x << " " << p.r.y << " " << p.r.z << " " << p.v.x << " " << p.v.y << " " << p.v.z << std::endl;
	}

	int k = 1;
	int b = 1;
	std::cout << "start calculating\n";

	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	COORD pos = { 0, 4 };

	for (double t = PPm::t_0 + PPm::dt; t <= PPm::t_1; t += PPm::dt)
	{
		auto th_n = KDK(particles, 4);
		std::thread th(set_dinamic_step,ref(particles));
		if (b % PPm::div == 0) {
			//std::cout << "t = " << t << "dt = " << PPm:: dt << std::endl;
			sistem_E.push_back(0);
			sistem_P.push_back(0);
			sistem_M.push_back(0);
			sistem_t.push_back(t);
			sistem_E_k.push_back(0);
			sistem_E_p.push_back(0);
			calculate_conversation_laws(particles, sistem_E.size() - 1, 0, particles.size());

			positions << t << std::endl;
			for (auto& p : particles)
			{
				positions << std::setprecision(13) << p.r.x << " " << p.r.y << " " << p.r.z << " " << p.v.x << " " << p.v.y << " " << p.v.z << std::endl;
			}
			k++;
		}
		th.join();
		b++;
		SetConsoleCursorPosition(hConsole, pos);
		std::cout << "\n";
		SetConsoleCursorPosition(hConsole, pos);
		std::cout <<std::setprecision(5)<< "[(t_n,t_1): " << t << "/" << PPm::t_1 << ", dt = " << PPm::dt << ", threads: " << th_n << "]";
	}
	std::cout << "\nend calculating\n";
	std::cout << "write conwersation laws into file\n";
	std::ofstream conversation_laws("C:\\Users\\mesho\\Desktop\\научка_2025_весна\\программная_реализация_Равновесная_Модель\\визуальзация_измерений\\measurements.txt");
	for (int i = 0; i < sistem_E.size(); i++)
	{
		conversation_laws << std::fixed << std::setprecision(13) << sistem_t[i] << " " << sistem_E[i] << " " << sistem_P[i] << " " << sistem_M[i] << " " << sistem_E_k[i] << " " << sistem_E_p[i] << " " << std::endl;
	}
	//====================================================================================================
	PythonWrapper py;
	py.vcl();
	py.vt();
	py.ptc();

}