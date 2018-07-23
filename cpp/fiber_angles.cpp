#include <iostream>
#include <cmath>
#include "zinc/context.hpp"
#include "zinc/element.hpp"
#include "zinc/field.hpp"
#include "zinc/fieldcache.hpp"
#include "zinc/fieldmodule.hpp"
#include "zinc/region.hpp"

using namespace OpenCMISS::Zinc;
double distance(double coord1[], double coord2[]){

	double focus = 35.25;
	double x1 = focus * std::cosh(coord1[0]) * std::cos(coord1[1]);
	double y1 = focus * std::sinh(coord1[0]) * std::sin(coord1[1]) * std::cos(coord1[2]);
	double z1 = focus * std::sinh(coord1[0]) * std::sin(coord1[1]) * std::sin(coord1[2]);
	double x2 = focus * std::cosh(coord2[0]) * std::cos(coord2[1]);
	double y2 = focus * std::sinh(coord2[0]) * std::sin(coord2[1]) * std::cos(coord2[2]);
	double z2 = focus * std::sinh(coord2[0]) * std::sin(coord2[1]) * std::sin(coord2[2]);
	double d = std::pow((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2), 0.5);
	return d;
}
void to_euclidian(double coords[], double out_values[]){

	double focus = 35.25;
	out_values[0] = focus * std::cosh(coords[0]) * std::cos(coords[1]);
	out_values[1] = focus * std::sinh(coords[0]) * std::sin(coords[1]) * std::cos(coords[2]);
	out_values[2] = focus * std::sinh(coords[0]) * std::sin(coords[1]) * std::sin(coords[2]);

}
void to_prolate_spheroidal(double coords[], double out_values[]){
	 
	double focus = 35.25;
	double x = coords[0];
	double y = coords[1];
	double z = coords[2];

	double a = std::sqrt(y*y + z*z + (x + focus)*(x + focus));
	double b = std::sqrt(y*y + z*z + (x - focus)*(x - focus));

	out_values[0] = std::acosh((1.0/(2*focus))*(a + b));
	out_values[1] = std::acos((1.0/(2*focus))*(a - b));
	out_values[2] = std::atan2(z,y);

}

void get_fiber_values(double coords[], double precision, double out_values[]){

	double d;
	//double out_values[3];
	double ps_coords[3];
	to_prolate_spheroidal(coords, ps_coords);

	Context context("heart");
	Region region = context.getDefaultRegion();
	region.readFile("heart1.exfile");
	Fieldmodule fieldmodule = region.getFieldmodule();
	Field field = fieldmodule.findFieldByName("coordinates");
	Field field_fibers = fieldmodule.findFieldByName("fibres");
	Mesh mesh = fieldmodule.findMeshByDimension(3);
	Fieldcache cache = fieldmodule.createFieldcache();
	Nodeset node_set = fieldmodule.findNodesetByName("nodes");
	Elementiterator ele_iter = mesh.createElementiterator();

	double minimum = 1000;
	int element_index = -1;
	double x = 0.5;
	double y = 0.5;
	double z = 0.5;
	double x_min = x;
	double y_min = y;
	double z_min = z;

	double local_node_coordinates[3] = {x, y, z};
	Element element = ele_iter.next();

	while(element.isValid()){

		cache.setMeshLocation(element, 3,  local_node_coordinates);
		field.evaluateReal(cache, 3, out_values);
		d = distance(out_values, ps_coords);
		if(d < minimum){
			minimum = d;
			element_index = element.getIdentifier();
		}
		element = ele_iter.next();
	}

	element = mesh.findElementByIdentifier(element_index);
	//#split the element in 8 smaller cubes, and search for the cube with the center closest
	//#to the input coordinates.
	//#continue with the closest cube until the given precision is reached.
	//#Example analogous in 2D: x marks the center of each of the four smaller squares.
	//#o marks the input coordinates.
	//#continue with the new square until the given precision is reached.
	//###########################################
	//#      |  x  |  x  |		 |  x  |  x  |#
	//#STEP1 |_____|_____|STEP2->|_____|_____|#
	//#      |  x  |  x  |		 |  x  |  x  |#
	//#      |_____|o____| 		 |o____|_____|#
	//###########################################
	int refinement = 2;
	double x1;
	double y1;
	double z1;
	while((d > precision) && (refinement < 50)){
		for (int i = -1; i <=1; i = i+2){
			for (int j = -1; j <= 1; j = j+2){
				for (int k = -1; k <= 1; k = k+2){
					x1 = x + i * std::pow(0.5,refinement);
					y1 = y + j * std::pow(0.5,refinement);
					z1 = z + k * std::pow(0.5,refinement);

					local_node_coordinates[0] = x1;
					local_node_coordinates[1] = y1;
					local_node_coordinates[2] = z1;

					cache.setMeshLocation(element, 3, local_node_coordinates);
					field.evaluateReal(cache, 3, out_values);
					d = distance(out_values, ps_coords);

					if (d < minimum){
						
						minimum = d;
						x_min = x1;
						y_min = y1;
						z_min = z1;
					}

				}
			}
		}
	x = x_min;
	y = y_min;
	z = z_min;
	std::cout << x << " " << y << " " << z << " " << std::endl;
	refinement += 1;	
	}
	local_node_coordinates[0] = x;
	local_node_coordinates[1] = y;
	local_node_coordinates[2] = z;
	cache.setMeshLocation(element, 3, local_node_coordinates);
	field_fibers.evaluateReal(cache, 3, out_values);
	field.evaluateReal(cache, 3, coords);//TEST
	std::cout << "distance: " << distance(coords, ps_coords) << std::endl;
}
int main(){

    
	double out_values[3];
	double coords[3] = {-22.762729289,11.152001135,-39.54201334};
	double precision = 0.00001;
    get_fiber_values(coords, precision, out_values);
    std::cout << "result:  " << out_values[0] << " " << out_values[1]<< " " << out_values[2] << std::endl << std::endl;

    return 0;
}
