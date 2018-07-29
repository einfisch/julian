#include "vtu_write.cpp"
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <typeinfo>
class vtu_converter{

	public:
		void heart_to_vtu (std::string, std::string, std::string);
		void fibers_to_vtu (std::string, std::string, double);

};

void vtu_converter::heart_to_vtu(std::string connectivity, std::string coordinates, std::string out_file_name){
	std::vector<int> num;
	num = heart_to_file(connectivity, coordinates);
	int num_points = num[0];
	int num_cells = num[1] + num[0];
	int region = num[2];
	create_file(coordinates, num_points, num_cells, region, out_file_name, false);

}

void vtu_converter::fibers_to_vtu(std::string coordinates, std::string out_file_name, double d){
	
	std::vector<int> num;
	num = fiber_directions(coordinates, d);
	std::cout << num.size() << std::endl;
	int num_points = num[0];
	int num_cells = num[1] + num[0];
	int region = num[2];

	create_file("vertices_fibers", num_points, num_cells, region, out_file_name, true);

}
