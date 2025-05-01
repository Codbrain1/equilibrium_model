#include<fstream>
#include<vector>
#include<random>
#include<iostream>
#include<cmath>
#include"Particle.h"
#include <iomanip>
#include<thread>
#include"wrapper.h"
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
				ps[i].F = ps[i].F + PPm::F(ps[i], ps[j]);
				ps[j].F = ps[i].F * (-1) + ps[j].F;
			}
		}
	}
}
void KDK(std::vector<PPm::Particle>& particles, int k, int i_0, int i_1)
{
	std::vector<vec> u_i(i_1 - i_0, vec(0, 0, 0));
	for (int i = i_0; i < i_1; i++)
	{
		u_i[i - i_0] = particles[i].v;
		particles[i].v = particles[i].v + particles[i].F * PPm::dt;
		particles[i].r = particles[i].r + (particles[i].v + u_i[i - i_0]) * PPm::dt * 0.5;
	}


	calc_forces(particles, i_0, i_1);

	for (int i = i_0; i < i_1; i++)
	{
		particles[i].v = (particles[i].v + u_i[i - i_0]) * 0.5 + particles[i].F * 0.5 * PPm::dt;
	}
}

std::vector<vec> k1(PPm::N,vec(0,0,0)), k2(PPm::N, vec(0, 0, 0)), k3(PPm::N, vec(0, 0, 0)), k4(PPm::N, vec(0, 0, 0));
std::vector<vec> k1v(PPm::N, vec(0, 0, 0)), k2v(PPm::N, vec(0, 0, 0)), k3v(PPm::N, vec(0, 0, 0)), k4v(PPm::N, vec(0, 0, 0));
// Рунге-Кутта 4 порядка
void RungeKutta4(std::vector<PPm::Particle>& particles, int k, int i_0, int i_1) {
	
	// Шаг 1
	calc_forces(particles, i_0, i_1);

	for (size_t i = i_0; i < i_1; ++i) {
		k1v[i] = particles[i].F * PPm::dt;
		k1[i] = particles[i].v * PPm::dt;
	}

	// Шаг 2
	for (size_t i = i_0; i < i_1; i++) {
		k2[i] = particles[i].r + k1[i] * 0.5;
		k2v[i] = particles[i].v + k1v[i] * 0.5;
	}

	calc_forces(particles, i_0, i_1);

	for (size_t i = 0; i < i_1; i++) {
		k2v[i] = particles[i].F * PPm::dt;
		k2[i] = (particles[i].v + k1[i] * 0.5) * PPm::dt;
	}

	// Шаг 3
	for (size_t i = 0; i < i_1; ++i) {
		k3[i] = particles[i].r + k2[i] * 0.5;
		k3v[i] = particles[i].v + k2[i] * 0.5;
	}
	for (int i = i_0; i < i_1; i++)
	{
		particles[i].F = vec(0, 0, 0);
	}

	for (size_t i = i_0; i < i_1; i++)
	{
		for (size_t j = 0; j < particles.size(); j++)
		{
			if (i != j)
			{
				particles[i].F = particles[i].F + PPm::F(particles[i], particles[j]);
			}
		}
	}

	for (size_t i = i_0; i < PPm::N; ++i) {
		k3v[i] = particles[i].F * PPm::dt;
		k3[i] = (particles[i].v + k2v[i] * 0.5) * PPm::dt;
	}

	// Шаг 4
	for (size_t i = i_0; i < i_1; i++) {
		k4[i] = particles[i].r + k3[i];
		k4v[i] = particles[i].v + k3v[i];
	}
	for (int i = i_0; i < i_1; i++)
	{
		particles[i].F = vec(0, 0, 0);
	}

	for (size_t i = i_0; i < i_1; i++)
	{
		for (size_t j = 0; j < particles.size(); j++)
		{
			if (i != j)
			{
				particles[i].F = particles[i].F + PPm::F(particles[i], particles[j]);
			}
		}
	}

	for (size_t i = i_0; i < i_1; i++) {
		k4v[i] = particles[i].F * PPm::dt;
		k4[i] = (particles[i].v + k3v[i]) * PPm::dt;
	}


	for (size_t i = i_0; i < i_1; i++) {
		particles[i].r = particles[i].r + (k1[i] + k2[i] * 2 + k3[i] * 2 + k4[i]) / 6.0;
		particles[i].v = particles[i].v + (k1v[i] + k2v[i] * 2 + k3v[i] * 2 + k4v[i]) / 6.0;
	}
}

void calculate_conversation_laws(std::vector<PPm::Particle>& particles, int k, int i_0, int i_1)
{
	for (size_t i = i_0; i < i_1; i++)
	{
		particles[i].E = 0;
		particles[i].P = vec(0, 0, 0);
		particles[i].M = vec(0, 0, 0);
	}

	for (int i = i_0; i < i_1; i++)
	{
		particles[i].E = 0.5 * particles[i].m * particles[i].v.module_2();	// set kinetical energy
		particles[i].P = particles[i].v * particles[i].m;					// set impulse
		particles[i].M = particles[i].r * particles[i].P;					// set moment impulse

		sistem_E_k[k] += particles[i].m * particles[i].v.module_2() * 0.5;

		for (int j = i+1; j < particles.size(); j++)
		{
			if (i != j)
			{
				vec r_ij = particles[j].r - particles[i].r;
				double r = r_ij.module();
				particles[i].E -= PPm::G * particles[i].m * particles[j].m / r;	//set potential energy
				sistem_E_p[k] -= PPm::G * particles[i].m * particles[j].m / r;
			}
		}
		sistem_E[k] += particles[i].E;
		sistem_P[k] = sistem_P[k] + particles[i].P.module();
		sistem_M[k] = sistem_M[k] + particles[i].M.module();
	}
}

double set_dinamic_step(std::vector<PPm::Particle>& ps)
{
	double min_dist = (ps[0].r - ps[1].r).module();
	double max_a = ps[0].F.module();
	for (int i = 0; i < ps.size(); i++)
	{
		for (int j = i + 1; j < ps.size(); j++)
		{
			double dist = (ps[i].r - ps[j].r).module();
			if (dist < min_dist)
				min_dist = dist;
		}
		double a = ps[i].F.module();
		if (max_a < a)
			max_a = a;
	}
	double new_dt = std::min(0.01*min_dist,0.1/sqrt(max_a));
	new_dt = std::max(0.5 * PPm::dt, std::min(2.0 * PPm::dt, new_dt));
	const double min_dt = 1e-10;
	const double max_dt = 0.01;
	return std::max(min_dt,std::min(new_dt,max_dt));
}
void set_initial_conditions(std::vector<PPm::Particle>& ps)
{
	using namespace Partcile_Particle_model;
	//setting the initial position
		// 
		//====================================================
	std::vector<double> radii = { 1, 1, 1, 1, 1, 1, 1};
	for (int i = 0; i < N; i++)
	{
		ps[i].r.x = radii[i] * cos((2 * PI * i) / (N));
		ps[i].r.y = radii[i] * sin((2 * PI * i) / (N));
		ps[i].r.z = 0;
		ps[i].m = M / N;
	}
	calc_forces(ps, 0, ps.size());
	//setting the initial velocity
	// 
	//==========================================================
	for (auto& i : ps) {

		double v_asimutal = sqrt(i.r.module() * i.F.module());		// (cylindrical coordinates)
		double phi = atan2(i.r.y, i.r.x);
		double r = i.r.module();

		i.v.x = -r * v_asimutal * sin(phi);		// (kartesian coordinates)
		i.v.y = r * v_asimutal * cos(phi);
		i.v.z = 0;
	}
}
void diagnostic(std::vector<PPm::Particle>& particles)
{
	// ===== ДИАГНОСТИКА НАЧАЛЬНЫХ УСЛОВИЙ =====
	std::cout << "\n=== Детальная проверка начальных условий ===\n";
	// 1. Вывод позиций и скоростей
	std::cout << "\nКоординаты и скорости частиц:\n";
	std::cout << std::setprecision(15);
	for (int i = 0; i < PPm::N; i++) {
		std::cout << "Частица " << i + 1 << ":\n";
		std::cout << "  Позиция: (" << particles[i].r.x << ", "
			<< particles[i].r.y << ", " << particles[i].r.z << ")\n";
		std::cout << "  Скорость: (" << particles[i].v.x << ", "
			<< particles[i].v.y << ", " << particles[i].v.z << ")\n";
		std::cout << "  Модуль скорости: " << particles[i].v.module() << "\n";
	}

	// 2. Вывод масс
	std::cout << "\nМассы частиц:\n";
	for (int i = 0; i < PPm::N; i++) {
		std::cout << "m" << i + 1 << " = " << particles[i].m << "\n";
	}

	// 3. Проверка сил
	std::cout << "\nСилы, действующие на частицы:\n";
	calc_forces(particles, 0, particles.size());
	for (int i = 0; i < PPm::N; i++) {
		std::cout << "F" << i + 1 << " = (" << particles[i].F.x << ", "
			<< particles[i].F.y << ", " << particles[i].F.z << ")\n";
		std::cout << "  |F| = " << particles[i].F.module() << "\n";
	}

	// 4. Детальный расчет энергии
	std::cout << "\nДетальный расчет энергии:\n";
	double total_K = 0, total_U = 0;

	// Кинетическая энергия
	std::cout << "Кинетическая энергия:\n";
	for (int i = 0; i < PPm::N; i++) {
		double K_i = 0.5 * particles[i].m * particles[i].v.module_2();
		std::cout << "  K" << i + 1 << " = " << K_i << "\n";
		total_K += K_i;
	}
	std::cout << "  Суммарная K = " << total_K << "\n";

	// Потенциальная энергия
	std::cout << "\nПотенциальная энергия (попарно):\n";
	for (int i = 0; i < PPm::N; i++) {
		for (int j = i + 1; j < PPm::N; j++) {
			vec r_ij = particles[j].r - particles[i].r;
			double dist = r_ij.module();
			double U_ij = -PPm::G * particles[i].m * particles[j].m / dist;
			std::cout << "  U" << i + 1 << j + 1 << " = " << U_ij
				<< " (r_" << i + 1 << j + 1 << " = " << dist << ")\n";
			total_U += U_ij;
		}
	}
	std::cout << "  Суммарная U = " << total_U << "\n";

	// Полная энергия
	double total_E = total_K + total_U;
	std::cout << "\nПолная энергия системы: E = K + U = "
		<< total_K << " + " << total_U << " = " << total_E << "\n";

	// 5. Проверка сохранения момента импульса
	std::cout << "\nМомент импульса:\n";
	vec total_L(0, 0, 0);
	for (int i = 0; i < PPm::N; i++) {
		vec L_i = particles[i].r * (particles[i].v * particles[i].m);
		std::cout << "  L" << i + 1 << " = (" << L_i.x << ", "
			<< L_i.y << ", " << L_i.z << ")\n";
		total_L = total_L + L_i;
	}
	std::cout << "  Суммарный L = (" << total_L.x << ", "
		<< total_L.y << ", " << total_L.z << ")\n";

	std::cout << "\n=== Проверка завершена ===\n\n";

	// ===== КОНЕЦ ДИАГНОСТИКИ =====
	// ===== КОНЕЦ ДИАГНОСТИКИ =====
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

	for (double t = PPm::t_0 + PPm::dt; t <= PPm::t_1; t += PPm::dt)
	{

		KDK(particles, k, 0, particles.size());
		PPm::dt = set_dinamic_step(particles);
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
		b++;

	}
	std::cout << "end calculating\n";
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
	//py.ptc();

}