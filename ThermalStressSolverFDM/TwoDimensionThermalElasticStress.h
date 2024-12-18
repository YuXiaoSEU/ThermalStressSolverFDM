//=====================================================================//
// Finite difference solver for two dimensional thermal elastic stress //
//=====================================================================//
// Author      : Pengxiao Wu                                           //
// Copyright   : All rights reserved.                                  //
// Description : C++, Ansi-style                                       //
// Last updated: 2024.12.18 by Pengxiao Wu                             //
//=====================================================================//

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <string>
#include <Windows.h>
#include <cstdlib> 
#include <direct.h>
#include <algorithm>

class TwoDimensionThermalElasticStress
{
protected:
	int _Mx, _My;
	int max_loops, max_frames;
	double delta_x, delta_t;

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
	void displacement_upadate();
	void compute_stress();

	void check_variables(int loops, int frame);

	void create_workspace_directory();
	void write_plt_file(int frame);

	int get_max_loops();
	int get_max_frames();

public:
	TwoDimensionThermalElasticStress();
	~TwoDimensionThermalElasticStress();
};