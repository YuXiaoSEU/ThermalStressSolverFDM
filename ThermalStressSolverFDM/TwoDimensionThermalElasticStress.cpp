//=============================================================//
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

			if (j == 0 && i > 0 && i < _Mx - 1)
			{
				u_new[index] = u_old[index_n] + 0.5 * (v_old[index_e] - v_old[index_w]);
				v_new[index] = v_old[index_n] + 0.5 * mu * (u_old[index_e] - u_old[index_w]) - alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			else if (j == 0 && i == 0)
			{
				//u_new[index] = u_old[index_n] + v_old[index_e] - v_old[index];
				//v_new[index] = v_old[index_n] + mu * (u_old[index_e] - u_old[index]) - alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			else if (j == 0 && i == _Mx - 1)
			{
				u_new[index] = u_old[index_n] + v_old[index] - v_old[index_w];
				v_new[index] = v_old[index_n] + mu * (u_old[index] - u_old[index_w]) - alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			else if (j == _My - 1 && i > 0 && i < _Mx - 1)
			{
				u_new[index] = u_old[index_s] - 0.5 * (v_old[index_e] - v_old[index_w]);
				v_new[index] = v_old[index_s] - 0.5 * mu * (u_old[index_e] - u_old[index_w]) + alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			else if (j == _My - 1 && i == 0)
			{
				//u_new[index] = u_old[index_s] - (v_old[index_e] - v_old[index]);
				//v_new[index] = v_old[index_s] - mu * (u_old[index_e] - u_old[index]) + alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			else if (j == _My - 1 && i == _Mx - 1)
			{
				u_new[index] = u_old[index_s] - (v_old[index] - v_old[index_w]);
				v_new[index] = v_old[index_s] - mu * (u_old[index] - u_old[index_w]) + alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
			}
			else if (i == _Mx - 1 && j > 0 && j < _My - 1)
			{
				u_new[index] = u_old[index_w] - 0.5 * mu * (v_old[index_n] - v_old[index_s]) + alpha * (T_old[index] - T_ini) * (1 + mu) * delta_x;
				v_new[index] = v_old[index_w] - 0.5 * (u_old[index_n] - u_old[index_s]);

			}
			else if (i == 0)
			{
				u_new[index] = 1e-18;
				v_new[index] = 1e-18;
			}
		}
	}
}

void TwoDimensionThermalElasticStress::displacement_upadate()
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

void TwoDimensionThermalElasticStress::compute_stress()
{

}

void TwoDimensionThermalElasticStress::check_variables(int loops, int frame)
{
	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(10) << "*** frame = " << frame << "\t\t";
	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(10) << "loops = " << loops << "\t\t";
	std::cout << "\n\n";
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
		printf("************** Worspace directory created in the root directory ***************\n");;
	}
}

void TwoDimensionThermalElasticStress::write_plt_file(int frame)
{
	std::stringstream output_filename;
	std::ofstream     output_file;
	output_filename.clear();
	output_filename.str("");
	output_filename << "./workspace/TwoDimensionThermalElasticStress/stress_" << frame << ".plt";
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