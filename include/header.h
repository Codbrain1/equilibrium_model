#pragma once
#include"Constants.h"
#include"vec.h"
#include"Particle.h"
#include<fstream>
#include<vector>
#include<random>
#include<iostream>
#include<cmath>
#include <iomanip>
#include<chrono>
#include<omp.h>

extern long double sistem_E;
extern vec sistem_P;
extern vec sistem_L;
extern long double sistem_t;
extern long double sistem_E_k;
extern long double sistem_E_p;
extern vec sistem_R;

namespace PPm = Particle_Particle_model;

long double sigma_exp(long double r);

long double random(long double beg, long double end);

void calc_forces(std::vector<PPm::Particle>& ps, int i_0, int i_1);

void calc_acceleration_asinc(std::vector<PPm::Particle>& ps, int i_0, int i_1);

void KDK_parallel(std::vector<PPm::Particle>& ps);

void KDK(std::vector<PPm::Particle>& ps, int i_0, int i_1);

void RungeKutta4(std::vector<PPm::Particle>& particles, int i_0, int i_1);

void calculate_conversation_laws_parallel(std::vector<PPm::Particle>& ps, int i_0, int i_1);

vec calculate_centr_mass(std::vector<PPm::Particle>& ps, int i_0, int i_1);

void calculate_conversation_laws(std::vector<PPm::Particle>& ps, int i_0, int i_1);

void set_dinamic_step_temp(std::vector<PPm::Particle>& ps);

void set_dinamic_step(std::vector<PPm::Particle>& ps);

void set_dinamic_step_parallel(std::vector<PPm::Particle>& ps);

void set_exponential_disk(std::vector<PPm::Particle>& ps, int i_0, int i_1, long double _x_, long double _y_);

void set_initial_circle(std::vector<PPm::Particle>& ps);
