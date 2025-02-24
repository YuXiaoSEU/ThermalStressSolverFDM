﻿//=============================================================//
// Visualization of the eight directions in a rectangular grid //
//                       nw----n----ne                         //
//                       |     |     |                         //
//                       w-----c-----e                         //
//                       |     |     |                         //
//                       sw----s----se                         //
//=============================================================//
#include "TwoDimensionThermalElasticStress.h"

void TwoDimensionThermalElasticStress::malloc_arrays()
{
	int size = _Mx * _My;

	u_new = new double[size];
	u_old = new double[size];
	v_new = new double[size];
	v_old = new double[size];
	T_new = new double[size];
	T_old = new double[size];
	T_0   = new double[size];

	sigma_x   = new double[size];
	sigma_y   = new double[size];
	tau_xy    = new double[size];
	von_mises = new double[size];
}

void TwoDimensionThermalElasticStress::delete_arrays()
{
	delete[] u_new;
	delete[] u_old;
	delete[] v_new;
	delete[] v_old;
	delete[] T_new;
	delete[] T_old;
	delete[] T_0;

	delete[] sigma_x;
	delete[] sigma_y;
	delete[] tau_xy;
	delete[] von_mises;
}

void TwoDimensionThermalElasticStress::read_parameters()
{
	std::string line;
	std::istringstream stream_temp;
	std::string fname("Parameters.par");
	std::ifstream infile;

	infile.open(fname.c_str());

	if (infile.fail())
	{
		std::cout << "The file Parameters.par does not exist!" << std::endl;
		exit(0);
	}
	printf("================== Parameters read from file Parameters.par ===================\n");

	// Mesh parameters, delta x, delta t
	getline(infile, line);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> _Mx;
	printf("%40s = %-16d \n", "_Mx", _Mx);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> _My;
	printf("%40s = %-16d \n", "_My", _My);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> delta_x;
	printf("%40s = %-16.8f \n", "delta_x", delta_x);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> delta_t;
	printf("%40s = %-16.8f \n", "delta_t", delta_t);

	// Output settings
	getline(infile, line);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> max_loops;
	printf("%40s = %-16d \n", "maximum loops", max_loops);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> max_frames;
	printf("%40s = %-16d \n", "maximum frames", max_frames);

	printf("-------------------------------------------------------------------------------\n");
	// Temperature field settings
	getline(infile, line);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> T_ini;
	printf("%40s = %-16.8f \n", "initial temperature", T_ini);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> rho;
	printf("%40s = %-16.8f \n", "density", rho);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> c_p;
	printf("%40s = %-16.8f \n", "specific heat", c_p);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> lambda;
	printf("%40s = %-16.8f \n", "coefficient of thermal conductivity", lambda);

	// Solid field settings
	getline(infile, line);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> E;
	printf("%40s = %-16e \n", "modulus of elasticity", E);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> mu;
	printf("%40s = %-16.8f \n", "poisson ratio", mu);
	getline(infile, line);
	stream_temp.clear();
	stream_temp.str(line);
	stream_temp >> alpha;
	printf("%40s = %-16.8f \n", "coefficient of linear expansion", alpha);

	printf("============================ Parameters read end ==============================\n\n");

	malloc_arrays();
	L = delta_x * _Mx;
	H = delta_x * _My;
	I = H * H * H / 12;
}

void TwoDimensionThermalElasticStress::initialize_domain()
{
	alpha_heat_dimless = lambda / (rho * c_p);
	alpha_heat_dimless /= (delta_x * delta_x / delta_t);
	for (int i = 0; i < _Mx; i++)
	{
		for (int j = 0; j < _My; j++)
		{
			int index = _My * i + j;

			u_new[index] = 0.0;
			u_old[index] = 0.0;
			v_new[index] = 0.0;
			v_old[index] = 0.0;
			T_new[index] = T_ini;
			T_old[index] = T_ini;
			T_0[index]   = T_ini;

			sigma_x[index]   = 0.0;
			sigma_y[index]   = 0.0;
			tau_xy[index]    = 0.0;
			von_mises[index] = 0.0;
		}
	}
}

void TwoDimensionThermalElasticStress::temperature_evolution()
{
#pragma omp parallel for
	for (int i = 1; i < _Mx - 1; i++)
	{
		for (int j = 1; j < _My - 1; j++)
		{
			int index   = _My * i + j;
			int index_e = _My * (i + 1) + j;
			int index_w = _My * (i - 1) + j;
			int index_n = _My * i + j + 1;
			int index_s = _My * i + j - 1;

			T_new[index] = T_old[index] \
					+ alpha_heat_dimless * (T_old[index_e] - T_old[index]) \
					+ alpha_heat_dimless * (T_old[index_w] - T_old[index]) \
					+ alpha_heat_dimless * (T_old[index_n] - T_old[index]) \
					+ alpha_heat_dimless * (T_old[index_s] - T_old[index]);
			}
		}
}

void TwoDimensionThermalElasticStress::temperature_boundary()
{
#pragma omp parallel for
	for (int i = 0; i < _Mx; i++)
	{
		for (int j = 0; j < _My; j++)
		{
			int index   = _My * i + j;
			int index_n = _My * i + j + 1;
			int index_s = _My * i + j - 1;
			int index_e = _My * (i + 1) + j;
			int index_w = _My * (i - 1) + j;

			if (i == 0)
			{
				T_new[index] = T_new[index_e];
			}
			else if (i == _Mx - 1)
			{
				T_new[index] = T_new[index_w];
			}
			else if (j == 0)
			{
				T_new[index] = 173;
			}
			else if (j == _My - 1)
			{
				T_new[index] = 373;
			}
		}
	}
}

void TwoDimensionThermalElasticStress::temperature_update()
{
#pragma omp parallel for
	for (int i = 0; i < _Mx; i++)
	{
		for (int j = 0; j < _My; j++)
		{
			int index = _My * i + j;

			T_old[index] = T_new[index];
		}
	}
}

void TwoDimensionThermalElasticStress::displacement_evolution()
{
#pragma omp parallel for
	for (int i = 1; i < _Mx - 1; i++)
	{
		for (int j = 1; j < _My - 1; j++)
		{
			int index    = _My * i + j;
			int index_w  = _My * (i - 1) + j;
			int index_e  = _My * (i + 1) + j;
			int index_n  = _My * i + j + 1;
			int index_s  = _My * i + j - 1;
			int index_ne = _My * (i + 1) + j + 1;
			int index_nw = _My * (i - 1) + j + 1;
			int index_se = _My * (i + 1) + j - 1;
			int index_sw = _My * (i - 1) + j - 1;

			u_new[index] = ((u_old[index_e] + u_old[index_w]) + (1 - mu) * 0.5 * (u_old[index_n] + u_old[index_s]) +
							(1 + mu) * 0.5 / 4 * (v_old[index_ne] + v_old[index_sw] - v_old[index_nw] - v_old[index_se])
							- alpha * (1 + mu) * (T_old[index_e] - T_old[index_w]) * 0.5 * delta_x) / (3 - mu);

			v_new[index] = ((v_old[index_n] + v_old[index_s]) + (1 - mu) * 0.5 * (v_old[index_e] + v_old[index_w]) +
							(1 + mu) * 0.5 / 4 * (u_old[index_ne] + u_old[index_sw] - u_old[index_nw] - u_old[index_se])
							- alpha * (1 + mu) * (T_old[index_n] - T_old[index_s]) * 0.5 * delta_x) / (3 - mu);
		}
	}
}

void TwoDimensionThermalElasticStress::displacement_boundary()
{
#pragma omp parallel for
	for (int i = 0; i < _Mx; i++)
	{
		for (int j = 0; j < _My; j++)
		{
			int index    = _My * i + j;
			int index_w  = _My * (i - 1) + j;
			int index_e  = _My * (i + 1) + j;
			int index_n  = _My * i + j + 1;
			int index_s  = _My * i + j - 1;
			int index_ne = _My * (i + 1) + j + 1;
			int index_nw = _My * (i - 1) + j + 1;
			int index_se = _My * (i + 1) + j - 1;
			int index_sw = _My * (i - 1) + j - 1;

			//left
			if (i == 0)
			{
				u_new[index] = 0.0;
				v_new[index] = 0.0;
			}
			//right
			else if (i == _Mx - 1 && j > 0 && j < _My - 1)
			{
				u_new[index] = u_old[index_w] - 0.5 * mu * (v_old[index_n] - v_old[index_s]) + alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
				v_new[index] = v_old[index_w] - 0.5 * (u_old[index_n] - u_old[index_s]);
			}
			//top
			else if (j == _My - 1 && i > 0 && i < _Mx - 1)
			{
				u_new[index] = u_old[index_s] - 0.5 * (v_old[index_e] - v_old[index_w]);
				v_new[index] = v_old[index_s] - 0.5 * mu * (u_old[index_e] - u_old[index_w]) + alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			//bottom
			else if (j == 0 && i > 0 && i < _Mx - 1)
			{
				u_new[index] = u_old[index_n] + 0.5 * (v_old[index_e] - v_old[index_w]);
				v_new[index] = v_old[index_n] + 0.5 * mu * (u_old[index_e] - u_old[index_w]) - alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			//upper left corner
			//else if (j == _My - 1 && i == 0)
			//{
				//u_new[index] = u_old[index_s] - (v_old[index_e] - v_old[index]);
				//v_new[index] = v_old[index_s] - mu * (u_old[index_e] - u_old[index]) + alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			//}
			//upper right corner
			else if (j == _My - 1 && i == _Mx - 1)
			{
				u_new[index] = 0.5 * (u_old[index_w] + u_old[index_s] + v_old[index_w] - v_old[index_s]);
				v_new[index] = 0.5 * (- u_old[index_w] + u_old[index_s] + v_old[index_w] + v_old[index_s]);
				//u_new[index] = u_old[index_s] - (v_old[index] - v_old[index_w]);
				//v_new[index] = v_old[index_s] - mu * (u_old[index] - u_old[index_w]) + alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			//lower left corner
			//else if (j == 0 && i == 0)
			//{
				//u_new[index] = u_old[index_n] + v_old[index_e] - v_old[index];
				//v_new[index] = v_old[index_n] + mu * (u_old[index_e] - u_old[index]) - alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			//}
			//lower right corner
			else if (j == 0 && i == _Mx - 1)
			{
				u_new[index] = 0.5 * (u_old[index_w] + u_old[index_n] + v_old[index_n] - v_old[index_w]);
				v_new[index] = 0.5 * (u_old[index_w] - u_old[index_n] + v_old[index_n] + v_old[index_w]);
				//u_new[index] = u_old[index_n] + v_old[index] - v_old[index_w];
				//v_new[index] = v_old[index_n] + mu * (u_old[index] - u_old[index_w]) - alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
		}
	}
}

void TwoDimensionThermalElasticStress::displacement_boundary_force()
{
	//vertical downward external combined force exerted on the upper right corner
	double F = 750.0;
#pragma omp parallel for
	for (int i = 0; i < _Mx; i++)
	{
		for (int j = 0; j < _My; j++)
		{
			int index = _My * i + j;
			int index_w = _My * (i - 1) + j;
			int index_e = _My * (i + 1) + j;
			int index_n = _My * i + j + 1;
			int index_s = _My * i + j - 1;
			int index_ne = _My * (i + 1) + j + 1;
			int index_nw = _My * (i - 1) + j + 1;
			int index_se = _My * (i + 1) + j - 1;
			int index_sw = _My * (i - 1) + j - 1;

			//left
			if (i == 0 && (j == _My / 2 || j == _My / 2 - 1))
			//if(i == 0)
			{
				u_new[index] = 0.0;
				v_new[index] = 0.0;
			}
			else if (i == 0 && j != _My / 2 && j != _My / 2 - 1)
			{
				double y = j * delta_x - H / 2;
				double local_force = calculate_force(F, y);
				double G = 0.5 * E / (1 + mu);

				u_new[index] = 0.0;
				v_new[index] = v_old[index_e] - local_force * delta_x / G;
			}
			//right
			else if (i == _Mx - 1 && j > 0 && j < _My - 1)
			{
				double y = j * delta_x - H / 2;
				double local_force = calculate_force(F, y);
				double G = 0.5 * E / (1 + mu);

				u_new[index] = u_old[index_w] - 0.5 * mu * (v_old[index_n] - v_old[index_s]) + alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
				v_new[index] = v_old[index_w] - 0.5 * (u_old[index_n] - u_old[index_s]) + delta_x * local_force / G;
			}
			//top
			else if (j == _My - 1 && i > 0 && i < _Mx - 1)
			{
				u_new[index] = u_old[index_s] - 0.5 * (v_old[index_e] - v_old[index_w]);
				v_new[index] = v_old[index_s] - 0.5 * mu * (u_old[index_e] - u_old[index_w]) + alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			//bottom
			else if (j == 0 && i > 0 && i < _Mx - 1)
			{
				u_new[index] = u_old[index_n] + 0.5 * (v_old[index_e] - v_old[index_w]);
				v_new[index] = v_old[index_n] + 0.5 * mu * (u_old[index_e] - u_old[index_w]) - alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			//upper right corner
			else if (j == _My - 1 && i == _Mx - 1)
			{
				double y = H / 2;
				double local_force = calculate_force(F, y);
				double G = 0.5 * E / (1 + mu);

				u_new[index] = (u_old[index_w] - mu * u_old[index_s] + mu * (v_old[index_s] - v_old[index_w]) + (alpha * (T_old[index] - T_ini) * (1 + mu) - mu * F / G) * delta_x) / (1 - mu);
				v_new[index] = (u_old[index_w] - u_old[index_s] + mu * v_old[index_s] - v_old[index_w] + (alpha * (T_old[index] - T_ini) * (1 + mu) - F / G) * delta_x) / (mu - 1);
			}
			//lower right corner
			else if (j == 0 && i == _Mx - 1)
			{
				double y = -H / 2;
				double local_force = calculate_force(F, y);

				double G = 0.5 * E / (1 + mu);

				u_new[index] = (u_old[index_w] - mu * u_old[index_n] + mu * (v_old[index_w] - v_old[index_n]) + (alpha * (T_old[index] - T_ini) * (1 + mu) + mu * F / G) * delta_x) / (1 - mu);
				v_new[index] = (u_old[index_w] - u_old[index_n] + v_old[index_w] - mu * v_old[index_n] + (alpha * (T_old[index] - T_ini) * (1 + mu) + F / G) * delta_x) / (1 - mu);
			}
		}
	}
}

void TwoDimensionThermalElasticStress::displacement_update()
{
#pragma omp parallel for
	for (int i = 0; i < _Mx; i++)
	{
		for (int j = 0; j < _My; j++)
		{
			int index = _My * i + j;
			u_old[index] = u_new[index];
			v_old[index] = v_new[index];
		}
	}
}

void TwoDimensionThermalElasticStress::calculate_stress()
{
#pragma omp parallel for
	for (int i = 0; i < _Mx; i++)
	{
		for (int j = 0; j < _My; j++)
		{
			int index = _My * i + j;
			int index_w = _My * (i - 1) + j;
			int index_e = _My * (i + 1) + j;
			int index_n = _My * i + j + 1;
			int index_s = _My * i + j - 1;
			int index_ne = _My * (i + 1) + j + 1;
			int index_nw = _My * (i - 1) + j + 1;
			int index_se = _My * (i + 1) + j - 1;
			int index_sw = _My * (i - 1) + j - 1;

			double du_dx;
			if (i == 0)
				du_dx = (u_new[index_e] - u_new[index]) / delta_x;
			else if (i == _Mx - 1) 
				du_dx = (u_new[index] - u_new[index_w]) / delta_x;
			else 
				du_dx = (u_new[index_e] - u_new[index_w]) / (2.0 * delta_x);

			double dv_dy;
			if (j == 0)
				dv_dy = (v_new[index_n] - v_new[index]) / delta_x;
			else if (j == _My - 1)
				dv_dy = (v_new[index] - v_new[index_s]) / delta_x;
			else
				dv_dy = (v_new[index_n] - v_new[index_s]) / (2.0 * delta_x);

			double du_dy;
			if (j == 0)
				du_dy = (u_new[index_n] - u_new[index]) / delta_x;
			else if (j == _My - 1)
				du_dy = (u_new[index] - u_new[index_s]) / delta_x;
			else 
				du_dy = (u_new[index_n] - u_new[index_s]) / (2.0 * delta_x);

			double dv_dx;
			if (i == 0)
				dv_dx = (v_new[index_e] - v_new[index]) / delta_x;
			else if (i == _Mx - 1)
				dv_dx = (v_new[index] - v_new[index_w]) / delta_x;
			else
				dv_dx = (v_new[index_e] - v_new[index_w]) / (2.0 * delta_x);

			double DeltaT = T_old[index] - T_ini;
			double T_term = alpha * DeltaT;

			sigma_x[index] = (E / (1.0 - mu * mu)) * (du_dx + mu * dv_dy) - (E * T_term) / (1.0 - mu);
			sigma_y[index] = (E / (1.0 - mu * mu)) * (dv_dy + mu * du_dx) - (E * T_term) / (1.0 - mu);
			tau_xy[index] = (E / (2.0 * (1.0 + mu))) * (du_dy + dv_dx);
			von_mises[index] = sqrt((sigma_x[index] + sigma_y[index]) * (sigma_x[index] + sigma_y[index]) - 3 * sigma_x[index] * sigma_y[index] + 3 * tau_xy[index] * tau_xy[index]);
		}
	}
}

double TwoDimensionThermalElasticStress::calculate_force(double F, double y)
{
	double force = -(0.5 * F / I) * (H * H / 4 - y * y);
	return force;
}

void TwoDimensionThermalElasticStress::check_variables(int loops, int frame, double error)
{
	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(10) << "*** frame = " << frame << "\t\t";
	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(10) << "loops = " << loops << "\t\t";
	std::cout << std::scientific << std::setprecision(10) << "max_error = " << error << "\t\t";
	std::cout << "\n\n";
}

double TwoDimensionThermalElasticStress::calculate_global_error()
{
	double max_error = 0.0;
	for (int i = 0; i < _Mx; i++)
	{
		for (int j = 0; j < _My; j++)
		{
			int index = _My * i + j;
			if (max_error < std::abs(u_new[index] - u_old[index]) || max_error < std::abs(v_new[index] - v_old[index]))
			{
				max_error = max(std::abs(u_new[index] - u_old[index]), std::abs(v_new[index] - v_old[index]));
			}
		}
	}
	return max_error;
}

void TwoDimensionThermalElasticStress::create_workspace()
{
	std::string directory = "./workspace";
	struct _stat info;
	if (_stat(directory.c_str(), &info) == 0 && (info.st_mode & _S_IFDIR))
		printf("************************ The workspace already exists *************************\n\n");
	else
	{
		if (_mkdir(directory.c_str()) != 0)
		{
			std::cerr << "Error: Unable to create workspace." << std::endl;
			exit(EXIT_FAILURE);
		}
		printf("******************* Workspace created in the root directory *******************\n\m");
	}
}

void TwoDimensionThermalElasticStress::create_workspace_directory()
{
	std::string directory = "./workspace/TwoDimensionThermalElasticStress";
	struct _stat info;
	if (_stat(directory.c_str(), &info) == 0 && (info.st_mode & _S_IFDIR))
		printf("******************* The workspace directory already exists ********************\n");
	else
	{
		if (_mkdir(directory.c_str()) != 0)
		{
			std::cerr << "Error: Unable to create workspace directory." << std::endl;
			exit(EXIT_FAILURE);
		}
		printf("************** Workspace directory created in the root directory **************\n");
	}
}

void TwoDimensionThermalElasticStress::write_plt_file(int loops)
{
	std::stringstream output_filename;
	std::ofstream     output_file;
	output_filename.clear();
	output_filename.str("");
	output_filename << "./workspace/TwoDimensionThermalElasticStress/stress_" << loops << ".plt";
	output_file.open(output_filename.str().c_str(), std::ios::trunc);

	output_file << "TITLE     = \" stress and displacement \"\n";
	output_file << "VARIABLES = \"X\", \"Y\", \"sigma_x\",\"sigma_y\",\"tau_xy \",\"von_mises\",\"u\",\"v\",\"temperature\"";
	output_file << "\n";
	output_file << "ZONE I=" << _Mx << ", J=" << _My << ", F=POINT\n";

	for (int i = 0; i < _Mx; i++)
	{
		for (int j = 0; j < _My; j++)
		{
			int index = _My * i + j;
			output_file << i << "\t" << j << "\t" << sigma_x[index] << "\t" << sigma_y[index] << "\t" << tau_xy[index] << "\t"
				<< von_mises[index] << "\t" << u_new[index]
				<< "\t" << v_new[index] << "\t"
				<< T_new[index];
			output_file << "\n";
		}
	}
	output_file.close();
}

int TwoDimensionThermalElasticStress::get_max_loops()
{
	return max_loops;
}

int TwoDimensionThermalElasticStress::get_max_frames()
{
	return max_frames;
}

TwoDimensionThermalElasticStress::TwoDimensionThermalElasticStress()
{
}

TwoDimensionThermalElasticStress::~TwoDimensionThermalElasticStress()
{
	delete_arrays();
}