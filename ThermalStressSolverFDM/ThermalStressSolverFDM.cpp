#include "TwoDimensionThermalElasticStress.h"

int main()
{
	TwoDimensionThermalElasticStress mystress;

	mystress.read_parameters();
	mystress.initialize_domain();
	mystress.create_workspace_directory();

	int loops = 0, frame = 0;
	int max_loops = mystress.get_max_loops();
	int max_frames = mystress.get_max_frames();

	printf("\n============================ Simulation begins ================================\n");
	mystress.write_plt_file(loops);
	while (loops < max_loops)
	{
		loops++;

		mystress.temperature_evolution();
		mystress.temperature_boundary();
		mystress.temperature_update();

		mystress.displacement_evolution();
		mystress.displacement_boundary();
		mystress.displacement_upadate();

		if (loops % (max_loops / max_frames) == 0)
		{
			frame = loops / (max_loops / max_frames);
			mystress.check_variables(loops, frame);
			mystress.write_plt_file(loops);
		}
	}
	printf("=============================== Simulation ends ===============================\n");
	std::cin.get();
	exit(0);
}