#include"header.h"
#include<algorithm>

//вычисление поверхностной плотности диска
//long double sigma(long double r) { return sigma_0 *r_alpha* r_alpha /r; }
long double sigma_exp(long double r) 
{ 
    using namespace PPm;
    return sigma_0 * expl(-r / r_alpha); 
}
//случайное числа в диапазоне [a,b] с равномерной вероятностью
long double random(long double beg, long double end)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<long double> dist(beg, end);
	return dist(gen);
}
//вычисление сил для всех частиц
void calc_forces(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
	for (int i = i_0; i < i_1; i++)
	{
		ps[i].F = vec(0.0L, 0.0L, 0.0L);
	}

	for (int i = i_0; i < i_1; i++)
	{
		for (int j = i + 1; j < ps.size(); j++)
		{
			if (i != j)
			{
				vec F_ij = PPm::F(ps[i], ps[j]);
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
		ps[i].a = vec(0.0L, 0.0L, 0.0L);

		for (int j = 0; j < ps.size(); j++)
		{
			if (i != j)
			{
				vec f_ij = PPm::calc_f(ps[i], ps[j]);
				ps[i].a = ps[i].a + f_ij;
			}
		}
	}
}
//метод интегрирования Kick-drift-kikck паралельно на процессоре
void KDK_parallel(std::vector<PPm::Particle>& ps) {
	const int N = ps.size();
	std::vector<vec> u_i(N);  // Сохраняем начальные скорости для всех частиц

	// Функция для обновления координат и скоростей (Drift + Kick)
	#pragma omp parallel for 
	for (int i = 0; i < N; ++i) {
		u_i[i] = ps[i].v + ps[i].a * PPm::dt;                     // Kick (обновление скорости)
		ps[i].r = ps[i].r + (ps[i].v + u_i[i]) * 0.5L * PPm::dt;  // Drift (обновление координат)
	}

	calc_acceleration_asinc(ps, 0, N);

	// Финальное обновление скоростей (Kick)
	#pragma omp parallel for 
	for (int i = 0; i < N; i++) {
		ps[i].v = (ps[i].v + u_i[i]) * 0.5L + ps[i].a * PPm::dt * 0.5L;
	}
}
//метод интегрирования Kick-drift-kikck
void KDK(std::vector<PPm::Particle>& ps, int i_0, int i_1)
{
	std::vector<vec> u_i(ps.size(), vec(0.0L, 0.0L, 0.0L));
	for (int i = i_0; i < i_1; i++) 
	{
		u_i[i - i_0] = ps[i].v;
		ps[i].v = ps[i].v + ps[i].a * PPm::dt; //Kick (обновление скорости)
		ps[i].r = ps[i].r + (ps[i].v + u_i[i - i_0]) * PPm::dt * 0.5L;//drift (обновление позиции)
	}
	calc_forces(ps, 0, ps.size());
	for (auto& p : ps)
		p.Evaluate_a(); //пересчет ускорений

	for (int i = 0; i < ps.size(); i++)
	{
		ps[i].v = (ps[i].v + u_i[i]) * 0.5L + ps[i].a * 0.5L * PPm::dt; //kick финальное обновление скорости
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
		sistem_E_k = 0.0L; //кинетическая энергия
		sistem_E_p = 0L; //потенциальная энергия
		sistem_E = 0.0L; //полная энергия
		sistem_P = vec(0.0L, 0.0L, 0.0L); //импульс
		sistem_L = vec(0.0L, 0.0L, 0.0L); //момент импульса
	}
	long double local_E_k = 0.0L, local_E_p = 0.0L;
	long double local_Px = 0.0L, local_Py = 0.0L, local_Pz = 0.0L;
    long double local_Lx = 0.0L, local_Ly = 0.0L, local_Lz = 0.0L;

	#pragma omp parallel for reduction(+:local_E_k, local_E_p, local_Px, local_Py, local_Pz, local_Lx, local_Ly, local_Lz)
	for (int i = i_0; i < i_1; i++)
	{
		long double E_k = 0.5L * ps[i].m * ps[i].v.module_2();

		ps[i].P = ps[i].v * ps[i].m;					// вычисляем импульс
		ps[i].L = ps[i].r * ps[i].P;					// вычисляем момент импульса

		local_E_k += E_k;
		
		for (int j = i + 1; j < ps.size(); j++) { //вычисляем потенциальную энергию
			vec r_ij = ps[j].r - ps[i].r;
			long double r = sqrtl(r_ij.module_2() + PPm::r_c);
			long double E_p = -PPm::G * ps[i].m * ps[j].m / r;
			local_E_p += E_p;
		}
		local_Px += ps[i].P.x;
        local_Py += ps[i].P.y;
        local_Pz += ps[i].P.z;
        local_Lx += ps[i].L.x;
        local_Ly += ps[i].L.y;
        local_Lz += ps[i].L.z;
	}
	vec local_P(local_Px, local_Py, local_Pz);
    vec local_L(local_Lx, local_Ly, local_Lz);
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
	vec R(0.0L, 0.0L, 0.0L);
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
		sistem_E_k = 0.0L;
		sistem_E_p = 0.0L;
		sistem_E = 0.0L;
		sistem_P = vec(0.0L, 0.0L, 0.0L);
		sistem_L = vec(0.0L, 0.0L, 0.0L);
	}

	for (int i = i_0; i < i_1; i++)
	{
		long double E_k = 0.5L * ps[i].m * ps[i].v.module_2();

		ps[i].P = ps[i].v * ps[i].m;					// set impulse
		ps[i].L = ps[i].r * ps[i].P;					// set moment impulse

		sistem_E_k += E_k;

		for (int j = i + 1; j < ps.size(); j++)
		{
			vec r_ij = ps[j].r - ps[i].r;
			long double r = sqrtl(r_ij.module_2() + PPm::r_c);
			long double E_p = -PPm::G * ps[i].m * ps[j].m / r;//set potential energy
			sistem_E_p += E_p;
		}
		sistem_P = sistem_P + ps[i].P;
		sistem_L = sistem_L + ps[i].L;
	}
	sistem_E = sistem_E_k + sistem_E_p;
}

//вычисление шага интегрирования (без ограничений, последовательно)
//TODO: проверит на корректность работы
void set_dynamic_step_parallel_temp(std::vector<PPm::Particle>& ps) {
    if (ps.size() < 2) return;

    // Параметры системы
    constexpr long double r_soft = 0.05L;  // Радиус гравитационного обрезания
    constexpr long double SAFETY_DIST = 0.001L;
    constexpr long double SAFETY_ACCEL = 0.005L;
    constexpr long double MIN_VELOCITY = 1.0e-50L;
    constexpr long double GRAVITY_FACTOR = 0.0001L;
    constexpr long double MAX_CHANGE_FACTOR = 2.0L;
    constexpr long double MIN_CHANGE_FACTOR = 0.5L;
	constexpr long double ADAPT_SMOOTH = 0.3L;
	
    // Инициализация
    long double min_dist = std::numeric_limits<long double>::max();
    long double max_a = 0.0L;
    long double max_v = 0.0L;
    long double min_orbital_period = std::numeric_limits<long double>::max();
    #pragma omp parallel for reduction(min: min_dist, min_orbital_period) \
                             reduction(max: max_a, max_v)
    for (int i = 0; i < ps.size(); i++) {
        const long double mi = static_cast<long double>(ps[i].m);
        const vec r_i = ps[i].r;
        
        for (int j = i + 1; j < ps.size(); j++) {
            // Расстояние с учётом мягкого радиуса
            const vec dr = r_i - ps[j].r;
            const long double dist = std::max(dr.module(), r_soft);
            min_dist = std::min(min_dist, dist);

            // Орбитальный период с поправкой на мягкий радиус
            const long double mu = static_cast<long double>(PPm::G) * (mi + static_cast<long double>(ps[j].m));
            const long double period = 2.0L * M_PIl * sqrtl(dist*dist*dist / mu);
            min_orbital_period = std::min(min_orbital_period, period);
        }

        // Максимальные ускорение и скорость
        max_a = std::max(max_a, static_cast<long double>(ps[i].a.module()));
        max_v = std::max(max_v, static_cast<long double>(ps[i].v.module()));
    }

    // Критерии с учётом обрезания
    const long double dt_dist = min_dist / (max_v + MIN_VELOCITY) * SAFETY_DIST;
    
    // Максимальное возможное ускорение из-за обрезания
    const long double max_effective_a = PPm::G * ps[0].m * ps.size() / (r_soft * r_soft);
    const long double dt_accel = (max_a > 0.0L) 
                               ? SAFETY_ACCEL / sqrtl(std::min(max_a, max_effective_a)) 
                               : LDBL_MAX;

    const long double dt_orbital = min_orbital_period * GRAVITY_FACTOR;

    // Выбираем минимальный шаг (исключаем dt_freefall, так как он учтён в dt_accel)
    const long double new_dt = std::min({dt_dist, dt_accel, dt_orbital});

    // Плавная адаптация
 
    long double adaptive_dt = ADAPT_SMOOTH * new_dt + 
                            (1.0L - ADAPT_SMOOTH) * static_cast<long double>(PPm::dt);

    // Финальные ограничения
    PPm::dt = std::clamp(adaptive_dt,
                        std::max(static_cast<long double>(PPm::mindt), 1.0e-30L),
                        std::min(static_cast<long double>(PPm::maxdt), 
                                min_orbital_period * 0.001L));
}
//вычисление шага интегрированая(ограничение - изменение не более чем в 2 раза за шаг, последовательно)
//TODO: проверить на корректность работу
void set_dynamic_step(std::vector<PPm::Particle>& ps)
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
	long double new_dt = std::min(0.001 * min_dist, 0.01 / sqrt(max_a));
	new_dt = std::max(0.1 * PPm::dt, std::min(10 * PPm::dt, new_dt));
	const long double min_dt = PPm::mindt;
	const long double max_dt = PPm::maxdt;
	PPm::dt = (std::max)(min_dt, (std::min)(new_dt, max_dt));
}
//вычисление шага интегрированая(ограничение - изменение не более чем в 10 раз за шаг, паралельно)
//TODO: проверить на корректность работу
void set_dynamic_step_parallel(std::vector<PPm::Particle>& ps)
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
	long double new_dt = (std::min)(0.001L * min_dist, 0.01L / sqrt(max_a));
	new_dt = std::max(PPm::dt/10.0L, std::min(10.0L* PPm::dt, new_dt));
	const long double min_dt = PPm::mindt;
	const long double max_dt = PPm::maxdt;
	PPm::dt = (std::max)(min_dt, (std::min)(new_dt, max_dt));
}
//функция задает экспонециальное распределение диска
void set_exponential_disk(std::vector<PPm::Particle>& ps, int i_0, int i_1, long double _x_, long double _y_)
{
	using namespace Particle_Particle_model;
	std::mt19937 gen(42.0L);
	std::uniform_real_distribution<long double> angle_dist(0.0L, 2.0L * PI);
	const long double r_core = 0.05L * r_alpha; // ядро смягчения потенциала
	int _N = i_1 - i_0;
	// Экспонциальное распределение с ядром
	for (int i = i_0; i < i_1; ++i) {
		long double u = (i - i_0 + 0.5L) / _N;
		long double r = -r_alpha * logl(1.0L - u);

		// Плавный переход в центре
		long double r_eff = sqrtl(r * r + r_core * r_core);
		long double theta = angle_dist(gen);

		// Кривая вращения с поправкой на ядро
		long double v_circ = sqrtl(M * powl(r_eff, 3.0L) / powl(r_eff + r_alpha, 3.0L));

		// Положение частицы
		ps[i].r.x = r * cosl(theta) + _x_;
		ps[i].r.y = r * sinl(theta) + _y_;
		ps[i].r.z = 0.0L;

		// Чисто круговая скорость (тангенциальная)
		ps[i].v.x = -v_circ * sinl(theta);  // v_φ = -v_circ * sin(θ)
		ps[i].v.y = v_circ * cosl(theta);   // v_φ = v_circ * cos(θ)
		ps[i].v.z = 0.0L;

		// Добавляем случайные компоненты (радиальную и тангенциальную)
		std::normal_distribution<long double> rad_vel(0.0L, 0.2L * v_circ);
		std::normal_distribution<long double> tang_vel(0.0L, 0.1L * v_circ);

		long double v_rad = rad_vel(gen);
		long double v_tang = tang_vel(gen);

		// Преобразование случайных компонент в декартовы координаты
		ps[i].v.x += v_rad * cosl(theta) - v_tang * sinl(theta);
		ps[i].v.y += v_rad * sinl(theta) + v_tang * cosl(theta);

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
		ps[i].r.x = R_max * cosl(2.0L * PI / (N)*i);
		ps[i].r.y = R_max * sinl(2.0L * PI / (N)*i);
		ps[i].r.z = 0.0L;
		ps[i].m = M /N;
	}

	calc_forces(ps, 0, ps.size());
	for (auto& i : ps)
		i.Evaluate_a();

	//setting the initial velocity
	// 
	//==========================================================
	for (auto& i : ps) {

		long double v_asimutal = sqrtl(i.r.module() * i.a.module());		// (cylindrical coordinates)
		long double phi = atan2l(i.r.y, i.r.x);
		long double r = i.r.module();

		i.v.x = -r * v_asimutal * sinl(phi);		// (kartesian coordinates)
		i.v.y = r * v_asimutal * cosl(phi);
		i.v.z = 0.0L;
	}
}

void set_bruk_orbit(std::vector<PPm::Particle> &ps)
{
	if(ps.size()!=3)
	{
		throw "Error size mass";
	}else{
		ps[0].r.x  = 0.8733047091;
		ps[1].r.x = -0.6254030288;
		ps[2].r.x = -0.2479016803;
		ps[0].r.y = 0;
		ps[1].r.y = 0;
		ps[1].r.y = 0;
		ps[0].v.x = 0;
		ps[1].v.x = 0;
		ps[2].v.x = 0;
		ps[0].v.y = 1.010776444;
		ps[1].v.y = -1.683353346;
		ps[2].v.y = 0.6725769022;
		ps[0].m=1;
		ps[1].m=1;
		ps[2].m=1;
		calc_forces(ps, 0, ps.size());
		for (auto& i : ps)
			i.Evaluate_a();
	}
}
