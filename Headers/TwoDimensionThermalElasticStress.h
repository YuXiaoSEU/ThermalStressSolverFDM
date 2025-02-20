//=====================================================================//
// Finite difference solver for two dimensional thermal elastic stress //
//=====================================================================//
// Author      : Pengxiao Wu                                           //
// Copyright   : All rights reserved.                                  //
// Description : C++                                                   //
// Last updated: 2025.02.18 by Pengxiao Wu                             //
//=====================================================================//

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <Windows.h>
#include <cstdlib>
#include <direct.h>

class TwoDimensionThermalElasticStress
{
protected:
	int _Mx, _My;
	int max_loops, max_frames;
	double delta_x, delta_t;
	double L, H, I;

	double T_ini, rho, c_p, lambda, alpha_heat_dimless;
	double E, mu, alpha;

	double* u_new, * u_old, * v_new, * v_old, * T_new, * T_old, * T_0;
	double* sigma_x, * sigma_y, * tau_xy, * von_mises;

protected:
	void malloc_arrays();
	void delete_arrays();

public:
	void read_parameters();
	void initialize_domain();

	void temperature_evolution();
	void temperature_boundary();
	void temperature_update();

	void displacement_evolution();
	void displacement_boundary();
	void displacement_boundary_force();
	void displacement_update();
	void calculate_stress();

	double calculate_force(double F, double y);

	void check_variables(int loops, int frame, double error);
	double calculate_global_error();

	void create_workspace();
	void create_workspace_directory();
	void write_plt_file(int loops);

	int get_max_loops();
	int get_max_frames();

public:
	TwoDimensionThermalElasticStress();
	~TwoDimensionThermalElasticStress();
};