#include "header.h"

long double sistem_E=0.0L;
vec sistem_P= vec(0.0L, 0.0L, 0.0L);;
vec sistem_L= vec(0.0L, 0.0L, 0.0L);;
long double sistem_t=0.0L;
long double sistem_E_k=0.0L;
long double sistem_E_p=0.0L;
vec sistem_R=  vec(0.0L, 0.0L, 0.0L);

int main()
{
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<PPm::Particle> particles(PPm::N); //array particals in model

	//setting the initial conditions
	// 
	//==========================================================


	set_exponential_disk(particles,0,PPm::N,0,0);
	//set_exponential_disk(particles,PPm::N/2,PPm::N,1,0);
	std::cout << "set initial conditions\n";

	calculate_conversation_laws(particles, 0, particles.size());
	sistem_R = calculate_centr_mass(particles, 0, particles.size());

	std::cout << "set intitial conversation laws\n";
	//==========================================================

// integration of differential equations
// 
//========================================================================================================================
	//write initial centr mass into file
	std::ofstream centr_mass("../результат_моделирования/centr_mass.txt");

	centr_mass<<std::setprecision(30) << sistem_t << " " << sistem_R.x << " " << sistem_R.y << " " << std::endl;

	//write initial positions into file
	std::ofstream positions("../результат_моделирования/positions.txt");

	positions << PPm::N << std::endl;
	positions << PPm::t_0 << std::endl;
	for (auto& p : particles)
	{
		positions << std::setprecision(30) << p.r.x << " " << p.r.y << " " << p.r.z << " " << p.v.x << " " << p.v.y << " " << p.v.z << std::endl;
	}

	//write initial conversation laws into file
	std::ofstream conversation_laws("../результат_моделирования/measurements.txt");

	conversation_laws << std::fixed << std::setprecision(30) << sistem_t << " " << sistem_E << " " << sistem_P.module() << " " << sistem_L.module() << " " << sistem_E_k << " " << sistem_E_p << " " << std::endl;


	//write execution time for initial calculating
	std::ofstream time("../результат_моделирования/time.txt");
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<long double> elapsed = end - start;
	time << elapsed.count() << " " << PPm::dt << std::endl;

	//	===BEGIN SIMULATION===
	int b = 1;
	std::cout << "start calculating\n";

	int save = 0;
	for (long double t = PPm::t_0 + PPm::dt; t <= PPm::t_1; t += PPm::dt)
	{
		#pragma omp single
		{
			sistem_E_k = 0.0L;
			sistem_E_p = 0.0L;
			sistem_E = 0.0L;
			sistem_P = vec(0.0L, 0.0L, 0.0L);
			sistem_L = vec(0.0L, 0.0L, 0.0L);
			sistem_R = vec(0.0L, 0.0L, 0.0L);
		}
		//KDK(particles, 0, particles.size());
		KDK_parallel(particles);
		set_dynamic_step_parallel_temp(particles);
		if (PPm::dt < 1e-8L)
		{
			PPm::div = 1000000;
		}
		else if (PPm::dt < 1e-7L)
		{
			PPm::div = 100000;
		}
		else if (PPm::dt < 1e-6L)
		{
			PPm::div = 10000;
		}
		else if (PPm::dt < 1e-5L)
		{
			PPm::div = 1000;
		}
		else if (PPm::dt < 1e-4L)
		{
			PPm::div = 100;
		}
		else if (PPm::dt < 1e-3L)
		{
			PPm::div = 1000;
		}
		else if (PPm::dt < 1e-2L)
		{
			PPm::div = 100;
		}
		else {
			PPm::div = 1;
		}
		if (b % PPm::div == 0) {
			sistem_t = t;
			calculate_conversation_laws_parallel(particles, 0, particles.size());
			sistem_R = calculate_centr_mass(particles, 0, PPm::N);

			positions << t << std::endl;
			for (auto& p : particles)
			{
				positions << std::setprecision(30) << p.r.x << " " << p.r.y << " " << p.r.z << " " << p.v.x << " " << p.v.y << " " << p.v.z << std::endl;
			}

			conversation_laws << std::fixed << std::setprecision(30) << sistem_t << " " << sistem_E << " " << sistem_P.module() << " " << sistem_L.module() << " " << sistem_E_k << " " << sistem_E_p << " " << std::endl;
			centr_mass <<std::setprecision(30)<< sistem_t << " " << sistem_R.x << " " << sistem_R.y << " " << std::endl;

			auto stop = std::chrono::high_resolution_clock::now();
			std::chrono::duration<long double> el = stop - start;
			time << el.count() << " " << PPm::dt << std::endl;
			++save;
		}
		b++;
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<long double> elapsed = end - start;
		std::cout << "\033[1F";
		std::cout << std::setprecision(5) << "[(t_n,t_1): " << t << "/" << PPm::t_1 << ", dt = " << PPm::dt << ", Время выполнения: " << elapsed.count() << " секунд," << "saves: " << save << "]\n";
	}
	std::cout << "\nend calculating\n";

	//	===END SIMULATION===

	return 0;
}