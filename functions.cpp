#include"header.h"

//вычисление поверхностной плотности диска
//long double sigma(long double r) { return sigma_0 *r_alpha* r_alpha /r; }
long double sigma_exp(long double r) 
{ 
    using namespace PPm;
    return sigma_0 * exp(-r / r_alpha); 
}
//случайное числа в диапазоне [a,b] с равномерной вероятностью
long double random(long double beg, long double end)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(beg, end);
	return dist(gen);
}
//вычисление сил для всех частиц
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
//паралельное вычисление ускорений для всех частиц на процессоре
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
//метод интегрирования Kick-drift-kikck паралельно на процессоре
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

	calc_acceleration_asinc(ps, 0, N);

	// Финальное обновление скоростей (Kick)
	#pragma omp parallel for 
	for (int i = 0; i < N; i++) {
		ps[i].v = (ps[i].v + u_i[i]) * 0.5 + ps[i].a * PPm::dt * 0.5;
	}
}
//метод интегрирования Kick-drift-kikck
void KDK(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
	std::vector<vec> u_i(ps.size(), vec(0, 0, 0));
	for (int i = i_0; i < i_1; i++) 
	{
		u_i[i - i_0] = ps[i].v;
		ps[i].v = ps[i].v + ps[i].a * PPm::dt; //Kick (обновление скорости)
		ps[i].r = ps[i].r + (ps[i].v + u_i[i - i_0]) * PPm::dt * 0.5;//drift (обновление позиции)
	}
	calc_forces(ps, 0, ps.size());
	for (auto& p : ps)
		p.Evaluate_a(); //пересчет ускорений

	for (int i = 0; i < ps.size(); i++)
	{
		ps[i].v = (ps[i].v + u_i[i]) * 0.5 + ps[i].a * 0.5 * PPm::dt; //kick финальное обновление скорости
	}
}
//TODO: Rungekutta4 находится в разработке
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
//вычисление сохраняющихся величин(паралельно на процессоре)
void calculate_conversation_laws_parallel(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
#pragma omp single
	{
		sistem_E_k = 0; //кинетическая энергия
		sistem_E_p = 0; //потенциальная энергия
		sistem_E = 0; //полная энергия
		sistem_P = vec(0, 0, 0); //импульс
		sistem_L = vec(0, 0, 0); //момент импульса
	}
	long double local_E_k = 0, local_E_p = 0;
	vec local_L = vec(0, 0, 0);
	vec local_P = vec(0, 0, 0);

#pragma omp parallel for reduction(+:local_E_k, local_E_p, local_P, local_L)
	for (int i = i_0; i < i_1; i++)
	{
		auto E_k = 0.5 * ps[i].m * ps[i].v.module_2();

		ps[i].P = ps[i].v * ps[i].m;					// вычисляем импульс
		ps[i].L = ps[i].r * ps[i].P;					// вычисляем момент импульса

		local_E_k += E_k;

		for (int j = i + 1; j < ps.size(); j++) { //вычисляем потенциальную энергию
			vec r_ij = ps[j].r - ps[i].r;
			long double r = sqrt(r_ij.module_2() + PPm::r_c);
			auto E_p = -PPm::G * ps[i].m * ps[j].m / r;
			local_E_p += E_p;
		}
		local_P = local_P + ps[i].P;
		local_L = local_L + ps[i].L;
	}
//обновление глобальных переменных
#pragma omp critical
	{
		sistem_E_k += local_E_k;
		sistem_E_p += local_E_p;
		sistem_E += local_E_k + local_E_p;
		sistem_P = sistem_P + local_P;
		sistem_L = sistem_L + local_L;
	}
}

//рассчет центра масс
vec calculate_centr_mass(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
	vec R(0, 0, 0);
	for (int i = i_0; i < i_1; i++)
	{
		R = R + ps[i].r * ps[i].m;
	}
	R = R / PPm::M;
	return R;
}
//последовательное вычисление законов сохранения
void calculate_conversation_laws(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
	#pragma omp single
	{
		sistem_E_k = 0;
		sistem_E_p = 0;
		sistem_E = 0;
		sistem_P = vec(0, 0, 0);
		sistem_L = vec(0, 0, 0);
	}

	for (int i = i_0; i < i_1; i++)
	{
		long double E_k = 0.5 * ps[i].m * ps[i].v.module_2();

		ps[i].P = ps[i].v * ps[i].m;					// set impulse
		ps[i].L = ps[i].r * ps[i].P;					// set moment impulse

		sistem_E_k += E_k;

		for (int j = i + 1; j < ps.size(); j++)
		{
			vec r_ij = ps[j].r - ps[i].r;
			long double r = sqrt(r_ij.module_2() + PPm::r_c);
			auto E_p = -PPm::G * ps[i].m * ps[j].m / r;//set potential energy
			sistem_E_p += E_p;
		}
		sistem_P = sistem_P + ps[i].P;
		sistem_L = sistem_L + ps[i].L;
	}
	sistem_E = sistem_E_k + sistem_E_p;
}

//вычисление шага интегрирования (без ограничений, последовательно)
//TODO: проверит на корректность работы
void set_dinamic_step_temp(std::vector<PPm::Particle>& ps)
{
	long double min_dist = (ps[0].r - ps[1].r).module();
	long double max_a = ps[0].a.module();

	// Находим минимальное расстояние между частицами и максимальное ускорение
	for (int i = 0; i < ps.size(); i++)
	{
		for (int j = i + 1; j < ps.size(); j++)
		{
			long double dist = (ps[i].r - ps[j].r).module();
			if (dist < min_dist)
				min_dist = dist;
		}
		long double a = ps[i].a.module();
		if (max_a < a)
			max_a = a;
	}

	// Вычисляем новый шаг интегрирования без ограничений на изменение
	long double new_dt = std::min(0.01 * min_dist, 0.1 / sqrt(max_a));

	// Применяем только глобальные ограничения (min_dt и max_dt)
	const long double min_dt = PPm::mindt;
	const long double max_dt = PPm::maxdt;
	PPm::dt = std::max(min_dt, std::min(new_dt, max_dt));
}
//вычисление шага интегрированая(ограничение - изменение не более чем в 2 раза за шаг, последовательно)
//TODO: проверить на корректность работу
void set_dinamic_step(std::vector<PPm::Particle>& ps)
{
	long double min_dist = (ps[0].r - ps[1].r).module();
	long double max_a = ps[0].a.module();
	for (int i = 0; i < ps.size(); i++)
	{
		for (int j = i + 1; j < ps.size(); j++)
		{
			long double dist = (ps[i].r - ps[j].r).module();
			if (dist < min_dist)
				min_dist = dist;
		}
		long double a = ps[i].a.module();
		if (max_a < a)
			max_a = a;
	}
	long double new_dt = std::min(0.01 * min_dist, 0.1 / sqrt(max_a));
	new_dt = std::max(0.5 * PPm::dt, std::min(2.0 * PPm::dt, new_dt));
	const long double min_dt = PPm::mindt;
	const long double max_dt = PPm::maxdt;
	PPm::dt = (std::max)(min_dt, (std::min)(new_dt, max_dt));
}
//вычисление шага интегрированая(ограничение - изменение не более чем в 10 раз за шаг, паралельно)
//TODO: проверить на корректность работу
void set_dinamic_step_parallel(std::vector<PPm::Particle>& ps)
{
	long double min_dist = (ps[0].r - ps[1].r).module();
	long double max_a = ps[0].a.module();
	#pragma omp parallel for reduction(min: min_dist) reduction(max: max_a)
	for (int i = 0; i < ps.size(); i++)
	{
		for (int j = i + 1; j < ps.size(); j++)
		{
			long double dist = (ps[i].r - ps[j].r).module();
			if (dist < min_dist)
				min_dist = dist;
		}
		long double a = ps[i].a.module();
		if (max_a < a)
			max_a = a;
	}
	long double new_dt = (std::min)(0.01 * min_dist, 0.1 / sqrt(max_a));
	new_dt = (std::max)(PPm::dt/10.0, (std::min)(10.0 * PPm::dt, new_dt));
	const long double min_dt = PPm::mindt;
	const long double max_dt = PPm::maxdt;
	PPm::dt = (std::max)(min_dt, (std::min)(new_dt, max_dt));
}
//функция задает экспонециальное распределение диска
void set_exponential_disk(std::vector<PPm::Particle>& ps, int i_0, int i_1, long double _x_, long double _y_)
{
	using namespace Particle_Particle_model;
	std::mt19937 gen(42);
	std::uniform_real_distribution<long double> angle_dist(0, 2 * PI);
	const long double r_core = 0.05 * r_alpha; // ядро смягчения потенциала
	int _N = i_1 - i_0;
	// Экспонциальное распределение с ядром
	for (int i = i_0; i < i_1; ++i) {
		long double u = (i - i_0 + 0.5) / _N;
		long double r = -r_alpha * log(1 - u);

		// Плавный переход в центре
		long double r_eff = sqrt(r * r + r_core * r_core);
		long double theta = angle_dist(gen);

		// Кривая вращения с поправкой на ядро
		long double v_circ = sqrt(M * pow(r_eff, 3) / pow(r_eff + r_alpha, 3));

		// Положение частицы
		ps[i].r.x = r * cos(theta) + _x_;
		ps[i].r.y = r * sin(theta) + _y_;
		ps[i].r.z = 0;

		// Чисто круговая скорость (тангенциальная)
		ps[i].v.x = -v_circ * sin(theta);  // v_φ = -v_circ * sin(θ)
		ps[i].v.y = v_circ * cos(theta);   // v_φ = v_circ * cos(θ)
		ps[i].v.z = 0;

		// Добавляем случайные компоненты (радиальную и тангенциальную)
		std::normal_distribution<long double> rad_vel(0, 0.2 * v_circ);
		std::normal_distribution<long double> tang_vel(0, 0.1 * v_circ);

		long double v_rad = rad_vel(gen);
		long double v_tang = tang_vel(gen);

		// Преобразование случайных компонент в декартовы координаты
		ps[i].v.x += v_rad * cos(theta) - v_tang * sin(theta);
		ps[i].v.y += v_rad * sin(theta) + v_tang * cos(theta);

		ps[i].m = M / N;
	}

	calc_forces(ps, i_0, i_1);
	for (int i = i_0; i < i_1; ++i) {
		ps[i].Evaluate_a();
	}
}
//задает распределение частиц по окружности
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
		ps[i].m = M / N;
	}

	calc_forces(ps, 0, ps.size());
	for (auto& i : ps)
		i.Evaluate_a();

	//setting the initial velocity
	// 
	//==========================================================
	for (auto& i : ps) {

		long double v_asimutal = sqrt(i.r.module() * i.a.module());		// (cylindrical coordinates)
		long double phi = atan2(i.r.y, i.r.x);
		long double r = i.r.module();

		i.v.x = -r * v_asimutal * sin(phi);		// (kartesian coordinates)
		i.v.y = r * v_asimutal * cos(phi);
		i.v.z = 0;
	}
}