#include <iostream>
#include <fstream>
#include <cmath>
#include "zinc/context.hpp"
#include "zinc/element.hpp"
#include "zinc/field.hpp"
#include "zinc/fieldcache.hpp"
#include "zinc/fieldmodule.hpp"
#include "zinc/region.hpp"
#include <string>
#include <vector>
#include <stdio.h>
#include <map>
#include <assert.h>
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
class exfile_reader{
		const char* in_file; int number_of_buckets;
	public:
		void set_file(const char*);
		void set_number_of_buckets(int);
		void create_node_coordinates(std::string);
		void create_node_connectivity(std::string);
		void create_directional_field(const char*, int , int , int , bool );
		void create_fiber_values(std::string);
		void create_interpolation(const char* , int , int , int );
};

void exfile_reader::set_file(const char* new_in_file){
	in_file = new_in_file;
}
void exfile_reader::set_number_of_buckets(int buckets){
	number_of_buckets = buckets;
}
void exfile_reader::create_node_coordinates(std::string out_file){
	//out_file = "/tmp/" + out_file;
	std::ofstream myfile;
	myfile.open(out_file, std::ios::trunc);

	Context context("heart");
	Region region = context.getDefaultRegion();
	region.readFile(in_file);
	Fieldmodule fieldmodule = region.getFieldmodule();
	Field field = fieldmodule.findFieldByName("coordinates");
	Fieldcache cache = fieldmodule.createFieldcache();

	Nodeset node_set = fieldmodule.findNodesetByName("nodes");
	Nodeiterator node_iter = node_set.createNodeiterator();
	int counter = 0;
	Node node = node_iter.next();
	double out_values[3];
	while(node.isValid()){
		//iterate over every node in the mesh and write the coordinates to the output file.
		cache.setNode(node);
		field.evaluateReal(cache, 3, out_values);
		myfile << out_values[0] << " " <<out_values[1] << " " << out_values[2] <<"\r\n";
		node = node_iter.next();
		counter += 1;
	}



}

void exfile_reader::create_node_connectivity(std::string out_file){
	std::ofstream myfile;
	myfile.open(out_file, std::ios:: trunc);
	Context context("heart");
	Region region = context.getDefaultRegion();
	region.readFile(in_file);
	Fieldmodule fieldmodule = region.getFieldmodule();
	Field field = fieldmodule.findFieldByName("coordinates");
	Fieldcache cache = fieldmodule.createFieldcache();

	Mesh mesh = fieldmodule.findMeshByDimension(3);
	Elementiterator ele_iter = mesh.createElementiterator();

	//int counter = 0;
	Element element = ele_iter.next();
	Elementfieldtemplate ele_template = element.getElementfieldtemplate(field, 1);
	int n = ele_template.getNumberOfLocalNodes();
	Node current_node;
	int local_index;
	while(element.isValid()){
		//iterate over every element and write the local nodes of the element to the output file
		for(int i = 1; i <= n; i++){
			current_node = element.getNode(ele_template, i);
			cache.setElement(element);
			local_index = current_node.getIdentifier();
			myfile << local_index << " ";

		}
		myfile << "\r\n";
		element = ele_iter.next();
		ele_template = element.getElementfieldtemplate(field,1);
	}
}

void exfile_reader::create_directional_field(const char* out_file, int refinement1, int refinement2, int refinement3, bool surface){
	int fine_mesh_size = 60*std::pow(2,refinement1 + refinement2 + refinement3);
	std::cout << "generating fibers..." << std::endl;
	std::ofstream myfile;
	myfile.open(out_file, std::ios::trunc);
	Context context("heart");
	Region region = context.getDefaultRegion();
	region.readFile(in_file);
	Fieldmodule fieldmodule = region.getFieldmodule();
	Field field = fieldmodule.findFieldByName("coordinates");
	Field field_fibers = fieldmodule.findFieldByName("fibres");
	Mesh mesh = fieldmodule.findMeshByDimension(3);
	Fieldcache cache = fieldmodule.createFieldcache();
	Nodeset node_set = fieldmodule.findNodesetByName("nodes");
	Nodetemplate node_template = node_set.createNodetemplate();
	node_template.defineField(field);
	node_template.defineField(field_fibers);

	Nodeiterator node_iter = node_set.createNodeiterator();
	//length of sides of the new elements in respect to local coordinate of parent element
	double l1 = 1.0/std::pow(2,refinement1);
	double l2 = 1.0/std::pow(2,refinement2);
	double l3 = 1.0/std::pow(2,refinement3);

	int counter = 0;
	
	int num_ele = mesh.getSize(); //number of elements contained in the mesh befor refinement.
	int mesh_size = num_ele;
	//create map "buckets", where each entry of buckets is a vector of Nodes.
	//each vector contains nodes with a similar sum od x,y and z coordinates.
	//if we want to create a new node and want to check if that node has already been created, we can 
	//limit our search to the bucket with similar nodes.
	std::map<int, std::vector<int>> buckets;

	Node node = node_iter.next();
	double global_coords[3];
	double out_values[3];
	double euclidian_coords[3];
	double coord_value;
	double min_buckets = 1000000;
	double max_buckets = -1000000;
	while(node.isValid()){
		counter += 1;
		cache.setNode(node);
		field.evaluateReal(cache, 3, out_values);
		to_euclidian(out_values, euclidian_coords);
		coord_value = euclidian_coords[0] + euclidian_coords[1] + euclidian_coords[2];
		if(coord_value < min_buckets){min_buckets = coord_value;}
		if(coord_value > max_buckets){max_buckets = coord_value;}
		node = node_iter.next();
	}
	//set the interval of each bucket as a uniform interval between the biggest and smallest values
	//of the initial nodes
	
	double bucket_interval = (max_buckets - min_buckets)/number_of_buckets;
	for(int i = 0; i < number_of_buckets; i++){
		buckets[i] = std::vector<int>{} ;
		
	}
	node_iter = node_set.createNodeiterator();
	node = node_iter.next();
	int c = 0;
	fieldmodule.beginChange();
	node_set.destroyAllNodes();
	int fibers_created = 0;
	Element element;
	
	Node current_node;
	Elementfieldtemplate ele_field_template;
	while(c < mesh_size){
		element = mesh.findElementByIdentifier(c + 1);
		ele_field_template = element.getElementfieldtemplate(field, 1);
		
		
		//iterate over xi_1, xi_2 and xi_3 axis, depending on number of refinements
		for(int i = 0; i < std::pow(2,refinement1); i++){
			for(int j = 0; j < std::pow(2,refinement2); j++){
				for(int k = 0; k < std::pow(2,refinement3); k++){
					//rendering of surface fibers or central fibers
					//local_node_coordinates are the coordinates at this element where one fiber is 
					//being visualized.
					double local_node_coordinates[3];
					if(surface){
						
						local_node_coordinates[0] = 0.5*l1+i*l1;
						local_node_coordinates[1] = 0.5*l2+j*l2;
						local_node_coordinates[2] = 1*l3+k*l3;
					}
					else{
						local_node_coordinates[0] = 0.5*l1+i*l1;
						local_node_coordinates[1] = 0.5*l2+j*l2;
						local_node_coordinates[2] = 0.5*l3+k*l3; }

					//fiber direction is given as a fiber angle between xi_1 axis in local coordinates
					//and fiber direction. Create a straight line that has slope-angle equal to fiber angle.
					//create additional node close to the node we want to evaluate and transform to global
					//coordinates.
					cache.setMeshLocation(element, 3, local_node_coordinates);
					field.evaluateReal(cache, 3, global_coords);
					to_euclidian(global_coords, euclidian_coords);
					field_fibers.evaluateReal(cache, 3, out_values);
					double fiber_angle = out_values[0];
					double slope = std::tan(fiber_angle);
					double d = 0.0001*std::cos(fiber_angle); 
					double xi_3 = local_node_coordinates[2];
					double xi_2 = local_node_coordinates[1] + slope*d;
					double xi_1 = local_node_coordinates[0] + d;
					double coord_value = euclidian_coords[0] + euclidian_coords[1] + euclidian_coords[2];
					
					
					int bucket_value = floor((coord_value - min_buckets)/bucket_interval);
					
					if(bucket_value < 0){bucket_value = 0;}
					if(bucket_value > number_of_buckets){bucket_value = number_of_buckets;}
					
					bool already_existing = false;
					/*if(!(bucket_value == 0)){
					 		for(std::size_t m = 0; m < buckets[bucket_value-1].size() ; m++){
					 			if((buckets[bucket_value-1].size()) == 0){break;}
					 			assert(m < buckets[bucket_value-1].size());
					 			Node node_test = node_set.findNodeByIdentifier(buckets[bucket_value-1][m]);
					 			cache.setNode(node_test);
					 			double coord_test[3];
					 			field.evaluateReal(cache, 3, coord_test);
					 			if(distance(global_coords, coord_test) < 0.0001){
					 				already_existing = true;
					 				
					 				break;
					 			}
					 		}
					 	}
						if(!(already_existing)){
							for(std::size_t m = 0; m < buckets[bucket_value].size(); m++){
								if((buckets[bucket_value].size()) == 0){break;}
					 			Node node_test = node_set.findNodeByIdentifier(buckets[bucket_value][m]);
					 			cache.setNode(node_test);
					 			double coord_test[3];
					 			field.evaluateReal(cache, 3, coord_test);
					 			if(distance(global_coords, coord_test) < 0.0001){
					 				already_existing = true;
					 				
					 				break;
					 			}
					 		}
						} 
						if(!(already_existing) && (bucket_value < number_of_buckets)){
							for(std::size_t m = 0; m < buckets[bucket_value  + 1].size(); m++){
								if((buckets[bucket_value+1].size()) == 0){break;}
					 			Node node_test = node_set.findNodeByIdentifier(buckets[bucket_value+1][m]);
					 			cache.setNode(node_test);
					 			double coord_test[3];
					 			field.evaluateReal(cache, 3, coord_test);
					 			if(distance(global_coords, coord_test) < 0.0001){
					 				already_existing = true;
					 				
					 				break;
					 			}
					 		}
						}
					*/
					if(already_existing == false){
						//if no node with the same coordinates is found, create a new one and use it.
						fibers_created += 1;
						node = node_set.createNode(counter + 1, node_template);
						cache.setNode(node);
						field.assignReal(cache,3, global_coords);
						field_fibers.assignReal(cache,3, out_values);
						buckets[bucket_value].push_back(node.getIdentifier());
						counter += 1;
						double new_coords[3] = {xi_1, xi_2, xi_3};
						cache.setMeshLocation(element, 3,new_coords);
						field.evaluateReal(cache, 3, new_coords);
						node = node_set.createNode(counter + 1, node_template);
						cache.setNode(node);
						field.assignReal(cache,3, new_coords);
						counter += 1;
					}

				}
			}
		}
		mesh.destroyElement(element);
		std::cout << std::flush << 100*fibers_created/fine_mesh_size << "% of all fibers created..." << "\r";
		c += 1;
	}
	for(int i = 1; i < 100; i++){
		Node node_dest = node_set.findNodeByIdentifier(i);
		node_set.destroyNode(node_dest);
	}	
	region.writeFile(out_file);
	fieldmodule.endChange();
	std::cout << "fibers created successfully!" << std::endl;
}

void exfile_reader::create_fiber_values(std::string out_file){
	std::ofstream myfile;
	myfile.open(out_file, std::ios::trunc);
	Context context("heart");
	Region region = context.getDefaultRegion();
	region.readFile(in_file);
	Fieldmodule fieldmodule = region.getFieldmodule();
	Field field = fieldmodule.findFieldByName("fibres");
	Fieldcache cache = fieldmodule.createFieldcache();
	Nodeset node_set = fieldmodule.findNodesetByName("nodes");
	Mesh mesh = fieldmodule.findMeshByDimension(3);
	Elementiterator ele_iter = mesh.createElementiterator();
	Nodeiterator node_iter = node_set.createNodeiterator();
	int counter = 0;
	Node node = node_iter.next();
	double out_values[3];
	std::vector<Node> vec;
	vec.push_back(node);

	while (node.isValid()){

		cache.setNode(node);
		field.evaluateReal(cache, 3, out_values);
		myfile << out_values[0] << " " <<out_values[1] << " " << out_values[2] <<"\r\n";
		node = node_iter.next();
		counter += 1;

	}
}

void exfile_reader::create_interpolation(const char* out_file, int refinement1, int refinement2, int refinement3){
	int fine_mesh_size = 60*std::pow(2,refinement3 + refinement2 + refinement1);
	std::cout << "generating mesh..." << std::endl;
	std::ofstream myfile;
	myfile.open(out_file, std::ios::trunc);
	Context context("heart");
	Region region = context.getDefaultRegion();
	region.readFile(in_file);
	Fieldmodule fieldmodule = region.getFieldmodule();
	Field field = fieldmodule.findFieldByName("coordinates");
	Field field_fibers = fieldmodule.findFieldByName("fibres");
	Mesh mesh = fieldmodule.findMeshByDimension(3);
	Fieldcache cache = fieldmodule.createFieldcache();
	Nodeset node_set = fieldmodule.findNodesetByName("nodes");
	Nodetemplate node_template = node_set.createNodetemplate();
	node_template.defineField(field);
	node_template.defineField(field_fibers);
	/*node_template.setValueNumberOfVersions(field,  1, node.VALUE_LABEL_VALUE, 1);//define, that d/ds1, d/ds2, d2/ds1ds2 should be saved in the node lamda
	node_template.setValueNumberOfVersions(field,  1, 3, 1);//mu
	node_template.setValueNumberOfVersions(field,  1, 4, 1);//theta*/
	Nodeiterator node_iter = node_set.createNodeiterator();
	//length of sides of the new elements in respect to local coordinate of parent element
	double l1 = 1.0/std::pow(2,refinement1);
	double l2 = 1.0/std::pow(2,refinement2);
	double l3 = 1.0/std::pow(2,refinement3);

	int counter = 0;
	
	int num_ele = mesh.getSize(); //number of elements contained in the mesh befor refinement.
	int mesh_size = num_ele;
	std::map<int, std::vector<int>> buckets;

	Node node = node_iter.next();
	double global_coords[3];
	double out_values[3];
	double euclidian_coords[3];
	double global_angles[3];
	double coord_value;
	double min_buckets = 1000000;
	double max_buckets = -1000000;
	while(node.isValid()){
		counter += 1;
		cache.setNode(node);
		field.evaluateReal(cache, 3, out_values);
		to_euclidian(out_values, euclidian_coords);
		coord_value = euclidian_coords[0] + euclidian_coords[1] + euclidian_coords[2];
		if(coord_value < min_buckets){min_buckets = coord_value;}
		if(coord_value > max_buckets){max_buckets = coord_value;}
		node = node_iter.next();
	}
	//set the interval of each bucket as a uniform interval between the biggest and smallest values
	//of the initial nodes
	
	double bucket_interval = (max_buckets - min_buckets)/number_of_buckets;

	for(int i = 0; i < number_of_buckets; i++){
		buckets[i] = std::vector<int>{} ;
		
	}
	int node_index;
	node_iter = node_set.createNodeiterator();
	node = node_iter.next();
	double bucket_value;
	while(node.isValid()){
		//insert initial nodes in buckets.
		//counter += 1;
		cache.setNode(node);
		field.evaluateReal(cache, 3, out_values);
		//
		to_euclidian(out_values, euclidian_coords);
		coord_value = euclidian_coords[0] + euclidian_coords[1] + euclidian_coords[2];
		bucket_value = floor((coord_value - min_buckets)/bucket_interval);
		node_index = node.getIdentifier();
		buckets[bucket_value].push_back(node_index);
		node = node_iter.next();
	}
	int c = 0;
	fieldmodule.beginChange(); //declare that the field is going to be manipulated
	Element element;
	Elementfieldtemplate ele_field_template;
	std::vector<std::vector<double>> local_node_coordinates;
	std::vector<int> node_identifiers;
	std::vector<int> node_identifiers_wedge;
	while(c < mesh_size){
		element = mesh.findElementByIdentifier(c + 1);
		ele_field_template = element.getElementfieldtemplate(field,1);
		for(int i = 0; i < std::pow(2,refinement1); i++){
			for(int j = 0; j < std::pow(2,refinement2); j++){
				for(int k = 0; k < std::pow(2,refinement3); k++){
					//we save the node indices of the new element here
					node_identifiers.clear();
					//define 8 nodes with their local element coordinates. They will form a new element.
					local_node_coordinates = {{0 + i*l1,0 + j*l2,0 + k*l3},
					 {l1 + i*l1 ,0 + j*l2,0 + k*l3},
					 {0 + i*l1,l2 + j*l2,0 + k*l3},
					 {l1 + i*l1,l2 + j*l2,0 + k*l3},
					 {0 + i*l1,0 + j*l2,l3 + k*l3},
					 {l1 + i*l1,0 + j*l2,l3 + k*l3},
					 {0 + i*l1 ,l2 + j*l2,l3 + k*l3},
					 {l1 + i*l1 ,l2 + j*l2,l3 + k*l3}};
					for(int n = 0; n < 8; n++){
					 	std::vector<double> local_coords_vec = local_node_coordinates[n];
					 	double local_coords[3]; //= {local_coords_vec[0], local_coords_vec[1], local_coords_vec[2]};
					 	local_coords[0] = local_coords_vec[0];
					 	local_coords[1] = local_coords_vec[1];
					 	local_coords[2] = local_coords_vec[2];
					 	cache.clearLocation();
					 	
					 	cache.setMeshLocation(element, 3, local_coords);
					 	field.evaluateReal(cache, 3, global_coords);
					 	//std::cout << global_coords[0]<<" " << global_coords[1] <<" " <<global_coords[2] << std::endl;
					 	field_fibers.evaluateReal(cache, 3, global_angles);
					 	to_euclidian(global_coords, euclidian_coords);

					 	coord_value = euclidian_coords[0] + euclidian_coords[1] + euclidian_coords[2];
					 	
					 	bucket_value = floor((coord_value - min_buckets)/bucket_interval);
					 	//std::cout << coord_value << " " << bucket_value << std::endl;
					 	
					 	if(bucket_value < 0){bucket_value = 0;}
					 	if(bucket_value > number_of_buckets){bucket_value = number_of_buckets;}

					 	bool already_existing = false;

					 	if(!(bucket_value == 0)){
					 		for(std::size_t m = 0; m < buckets[bucket_value-1].size() ; m++){
					 			if((buckets[bucket_value-1].size()) == 0){break;}
					 			assert(m < buckets[bucket_value-1].size());
					 			Node node_test = node_set.findNodeByIdentifier(buckets[bucket_value-1][m]);
					 			cache.setNode(node_test);
					 			double coord_test[3];
					 			field.evaluateReal(cache, 3, coord_test);
					 			if(distance(global_coords, coord_test) < 0.0001){
					 				already_existing = true;
					 				node_identifiers.push_back(node_test.getIdentifier());
					 				break;
					 			}
					 		}
					 	}
						if(!(already_existing)){
							for(std::size_t m = 0; m < buckets[bucket_value].size(); m++){
								if((buckets[bucket_value].size()) == 0){break;}
					 			Node node_test = node_set.findNodeByIdentifier(buckets[bucket_value][m]);
					 			cache.setNode(node_test);
					 			double coord_test[3];
					 			field.evaluateReal(cache, 3, coord_test);
					 			if(distance(global_coords, coord_test) < 0.0001){
					 				already_existing = true;
					 				node_identifiers.push_back(node_test.getIdentifier());
					 				break;
					 			}
					 		}
						} 
						if(!(already_existing) && (bucket_value < number_of_buckets)){
							for(std::size_t m = 0; m < buckets[bucket_value  + 1].size(); m++){
								if((buckets[bucket_value+1].size()) == 0){break;}
					 			Node node_test = node_set.findNodeByIdentifier(buckets[bucket_value+1][m]);
					 			cache.setNode(node_test);
					 			double coord_test[3];
					 			field.evaluateReal(cache, 3, coord_test);
					 			if(distance(global_coords, coord_test) < 0.0001){
					 				already_existing = true;
					 				node_identifiers.push_back(node_test.getIdentifier());
					 				break;
					 			}
					 		}
						}
						if(already_existing == false){
							node = node_set.createNode(counter + 1, node_template);
							//std::cout << node.isValid() << std::endl;
							cache.setNode(node);
							field.assignReal(cache, 3, global_coords);
							field_fibers.assignReal(cache, 3, global_angles);
							node_identifiers.push_back(counter + 1);
							buckets[bucket_value].push_back(node.getIdentifier());
							//node = buckets[bucket_value][buckets[bucket_value].size()];

							counter += 1;

							}
						}
					//Elementbasis = fieldmodule.createElementbasis(2,1);
						Elementtemplate element_template = mesh.createElementtemplate();
						if(node_identifiers[0] == node_identifiers[1]){
							node_identifiers_wedge.clear();
							node_identifiers_wedge.push_back(node_identifiers[0]);
							node_identifiers_wedge.push_back(node_identifiers[2]);
							node_identifiers_wedge.push_back(node_identifiers[3]);
							node_identifiers_wedge.push_back(node_identifiers[4]);
							node_identifiers_wedge.push_back(node_identifiers[6]);
							node_identifiers_wedge.push_back(node_identifiers[7]);
							node_identifiers_wedge.push_back(0);
							node_identifiers_wedge.push_back(0);
							} 
						
						

						
							//set shape and properties of the element to be created..
						int node_indexes[8] = {1,2,3,4,5,6,7,8};
						Elementbasis basis_xi_3;
						Elementbasis basis_xi_2;
						Elementbasis basis_xi_1;
						basis_xi_3 = fieldmodule.createElementbasis(3, basis_xi_3.FUNCTION_TYPE_CUBIC_HERMITE);
						basis_xi_2 = fieldmodule.createElementbasis(3, basis_xi_2.FUNCTION_TYPE_LINEAR_LAGRANGE);
						basis_xi_1 = fieldmodule.createElementbasis(3, basis_xi_2.FUNCTION_TYPE_LINEAR_LAGRANGE);

						element_template.setElementShapeType(element.SHAPE_TYPE_CUBE);
						element_template.setNumberOfNodes(8);

						element_template.defineFieldSimpleNodal(field, 1, basis_xi_3, 8, node_indexes);
						element_template.defineFieldSimpleNodal(field, 2, basis_xi_2, 8, node_indexes);
						element_template.defineFieldSimpleNodal(field, 3, basis_xi_1, 8, node_indexes);
						int z = 1;
						if(node_identifiers[0] == node_identifiers[1]){
							for(int identifier = 0; identifier < 8; identifier++){
								node = node_set.findNodeByIdentifier(node_identifiers_wedge[identifier]);
								element_template.setNode(z,node);
								z += 1;
								}
							}
						else{
							for(int identifier = 0; identifier < 8; identifier++){
								node = node_set.findNodeByIdentifier(node_identifiers[identifier]);
								element_template.setNode(z,node);
								z += 1;
								}
							}
						mesh.defineElement(num_ele + 1 , element_template);
						num_ele += 1;
						

					}
				}
		
			
		
			}
		std::cout << std::flush << 100*(num_ele - 60)/fine_mesh_size << "% of all elements created..." << "\r";
		mesh.destroyElement(element);
		c += 1;
		}
		region.writeFile(out_file);
		fieldmodule.endChange();
		std::cout << "generation of mesh with refinement: " << refinement1 << ", " << refinement2 << ", " << refinement3
		<< " finished..." << std::endl;
		//for(int i = 0; i < number_of_buckets; i++){
		//	std::cout << buckets[i].size() << std::endl;
		//}
}
/*int main(){

	exfile_reader ex_reader;
	ex_reader.set_number_of_buckets(10000);
	ex_reader.set_file("heart1.exfile");
	ex_reader.create_node_coordinates("test.txt");
	ex_reader.create_node_connectivity("test1.txt");
	ex_reader.create_fiber_values("test2.txt");
	ex_reader.create_directional_field("fibers.exfile", 2,2,0, true);
	ex_reader.create_interpolation("heart_444.exfile" , 4,4,4);
	return 1; 
}*/