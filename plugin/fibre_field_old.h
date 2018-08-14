/*
 * Copyright (c) 2018:  G-CSC, Goethe University Frankfurt
 * Authors: Julian Hilbert
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 *      example usage (on unit square)
 *      	elemDisc = ConvectionDiffusion("c", "Inner", disc)
 *			elemDisc:set_diffusion(CardiacFibreField(100, 0, 1, 0.01))
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__FIBRE_FIELD__
#define __H__UG__LIB_DISC__SPATIAL_DISC__FIBRE_FIELD__

// extern headers
#include "fiber_angles.cpp"
#include <vector>
#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"
#include <cmath>
#include "zinc/context.hpp"
#include "zinc/element.hpp"
#include "zinc/field.hpp"
#include "zinc/fieldcache.hpp"
#include "zinc/fieldmodule.hpp"
#include "zinc/region.hpp"
#include <map>

namespace ug{
	namespace OpenCMISS{

template <typename TData ,int dim, typename TRet = void>
class CardiacFibreField
		: public StdGlobPosData<CardiacFibreField<TData, dim, TRet >, TData, dim, TRet>
/*template <int dim>
class CardiacFibreField
		: public StdGlobPosData<CardiacFibreField<MathVector<dim>, dim, void>, MathVector<dim>, dim, void>*/
{

public:

	
	::OpenCMISS::Zinc::Region region;
	
	::OpenCMISS::Zinc::Context context;
	std::map<int, std::vector<std::vector<double>>> buckets;
	double min_buckets;
	double max_buckets;
	double bucket_interval;
	int number_of_buckets;

	CardiacFibreField(std::string filename)
	: context("heart")
	{	const char* c = filename.c_str();
		::OpenCMISS::Zinc::Field field;
		::OpenCMISS::Zinc::Field field_fibers;
		::OpenCMISS::Zinc::Mesh mesh;
		::OpenCMISS::Zinc::Fieldcache cache;
		::OpenCMISS::Zinc::Nodeset node_set;
		::OpenCMISS::Zinc::Elementiterator ele_iter;
		::OpenCMISS::Zinc::Element element;
		::OpenCMISS::Zinc::Fieldmodule fieldmodule;
		std::cerr << "context is valid: " << context.isValid() << std::endl;
		region = context.getDefaultRegion();
		std::cerr << "read file: " << region.readFile(c) << std::endl;
		fieldmodule = region.getFieldmodule();
		field = fieldmodule.findFieldByName("coordinates");
		field_fibers = fieldmodule.findFieldByName("fibres");
		
		::OpenCMISS::Zinc::Node node;
		mesh = fieldmodule.findMeshByDimension(3);
		cache = fieldmodule.createFieldcache();
		node_set = fieldmodule.findNodesetByName("nodes");
		::OpenCMISS::Zinc::Nodeiterator node_iter = node_set.createNodeiterator();
		ele_iter = mesh.createElementiterator();
		::OpenCMISS::Zinc::Nodetemplate node_template = node_set.createNodetemplate();
		node_template.defineField(field);
		node_template.defineField(field_fibers);
		int refinement = 4;
		
		int fine_mesh_size = 60*std::pow(2,3*refinement);
		int counter = 0;
		//int node_index;
		number_of_buckets = 100000;
		double l1 = 1.0/std::pow(2,refinement+1);
		double l2 = 1.0/std::pow(2,refinement+1);
		double l3 = 1.0/std::pow(2,refinement);
		min_buckets = 1000000;
		max_buckets = -1000000;
		double coord_value;
		double out_values[3];
		double euclidian_coords[3];
		double global_coords[3];
		double global_angles[3];
		node = node_iter.next();
		while(node.isValid()){
			counter += 1;
			cache.setNode(node);
			field.evaluateReal(cache, 3, out_values);
			to_euclidian(out_values, euclidian_coords);
			coord_value = euclidian_coords[0] + euclidian_coords[1] + euclidian_coords[2];
			//std::cerr << node.getIdentifier() << std::endl;
			if(coord_value < min_buckets){min_buckets = coord_value;}
			if(coord_value > max_buckets){max_buckets = coord_value;}
			node = node_iter.next();
		}
		std::vector<double> node_element;
		
		//min_buckets = min_buckets - 2;
		bucket_interval = (max_buckets - min_buckets)/number_of_buckets;
		for(int i = 0; i < number_of_buckets; i++){
			std::vector<std::vector<double>> dummy_vec;
			buckets[i] = dummy_vec;
		}
		node_iter = node_set.createNodeiterator();
		element = ele_iter.next();
		node = node_iter.next();
		double bucket_value;
		while(node.isValid()){

			cache.setNode(node);
			field.evaluateReal(cache, 3, out_values);
			//
			to_euclidian(out_values, euclidian_coords);
			coord_value = euclidian_coords[0] + euclidian_coords[1] + euclidian_coords[2];
			bucket_value = floor((coord_value - min_buckets)/bucket_interval);
			//node_index = node.getIdentifier();
			//buckets[bucket_value].push_back(node_index);
			//std::cerr << buckets[bucket_value].size() << " "<< bucket_value << std::endl;
			node = node_iter.next();
		}
		int num_ele = mesh.getSize();
		
		node_element.push_back(-1);
		node_element.push_back(-1);
		node_element.push_back(-1);
		node_element.push_back(-1);
		node_element.push_back(-1);
		std::vector<std::vector<double>> local_node_coordinates;
		while(element.isValid()){
		
			for(int i = 0; i < std::pow(2,refinement+1); i++){
				for(int j = 0; j < std::pow(2,refinement+1); j++){
					for(int k = 0; k < std::pow(2,refinement); k++){
						//we save the node indices of the new element here
						
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
						 			
						 			Node node_test = node_set.findNodeByIdentifier(buckets[bucket_value-1][m][0]);
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
						 			Node node_test = node_set.findNodeByIdentifier(buckets[bucket_value][m][0]);
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
						 			Node node_test = node_set.findNodeByIdentifier(buckets[bucket_value+1][m][0]);
						 			cache.setNode(node_test);
						 			double coord_test[3];
						 			field.evaluateReal(cache, 3, coord_test);
						 			if(distance(global_coords, coord_test) < 0.0001){
						 				already_existing = true;
						 				
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
								node_element[0] = node.getIdentifier();
								node_element[1] = element.getIdentifier();
								node_element[2] = local_coords_vec[0];
								node_element[3] = local_coords_vec[1];
								node_element[4] = local_coords_vec[2];
								buckets[bucket_value].push_back(node_element);
								//node = buckets[bucket_value][buckets[bucket_value].size()];

								counter += 1;

								}
							}
						num_ele++;
						}
					}
				}
			std::cout << std::flush << 100*(num_ele - 60)/fine_mesh_size << "% of preprocessing..." << "\r";
			element = ele_iter.next();
			}
		//int max_bucket = 0;
		/*std::cerr << min_buckets << std::endl;
		for(int i = 0; i < number_of_buckets; i++){
			if(buckets[i].size() == 0){max_bucket++; std::cerr << max_bucket << std::endl;}
		}*/
		std::cout << "PREPROCESSINGFINISHED" << std::endl;
	};



	virtual ~CardiacFibreField()
	{

	}

	/// maps nu: x -> nu(x)
	inline TRet evaluate(TData& nu, const MathVector<dim>& x_vec_in, number time, int si) const
	{	::OpenCMISS::Zinc::Field field;
		::OpenCMISS::Zinc::Field field_fibers;
		::OpenCMISS::Zinc::Mesh mesh;
		::OpenCMISS::Zinc::Fieldcache cache;
		::OpenCMISS::Zinc::Nodeset node_set;
		::OpenCMISS::Zinc::Elementiterator ele_iter;
		::OpenCMISS::Zinc::Element element;
		::OpenCMISS::Zinc::Fieldmodule fieldmodule;
		fieldmodule = region.getFieldmodule();
		field = fieldmodule.findFieldByName("coordinates");
		field_fibers = fieldmodule.findFieldByName("fibres");
		mesh = fieldmodule.findMeshByDimension(3);
		cache = fieldmodule.createFieldcache();
		node_set = fieldmodule.findNodesetByName("nodes");

		::OpenCMISS::Zinc::Node node;
		double xs[3];
		double xps[3];
		for (int i = 0; i < 3; i++){
			xs[i] = x_vec_in[i];
		}
		
		to_prolate_spheroidal(xs, xps);
		//std::cout << "coordinates of node: "<< xps[0] << " " << xps[1] << " " << xps[2] << std::endl;
		double x[3];
		x[0] = xs[0];
		x[1] = xs[1];
		x[2] = xs[2];
		double out[3];
		double x1[3];
		double d;

		
		//std::cerr << x_vec_in << std::endl;

		to_euclidian(xps, x);
		double coord_value = x[0] + x[1] + x[2];
		double bucket_value = floor((coord_value - min_buckets)/bucket_interval);
		//double precision = 0.01;
		double min = 1000;
		double local_coords[3];
		int node_index;
		int ele_index;
		int j;
		int iter;
		std::size_t i;
		for(i = 0; i < buckets.at(bucket_value).size(); i++){
				node = node_set.findNodeByIdentifier(buckets.at(bucket_value)[i][0]);
				cache.setNode(node);
				field.evaluateReal(cache, 3, out);
				//to_prolate_spheroidal(x, xps);
				d = distance(out, xps);
				
				if(d < min){
					node_index = node.getIdentifier(); 
					min = d;
					//std::cerr << "distance is "<< " " << d << std::endl;
					ele_index =  buckets.at(bucket_value)[i][1];
					local_coords[0] = buckets.at(bucket_value)[i][2];
					local_coords[1] = buckets.at(bucket_value)[i][3];
					local_coords[2] = buckets.at(bucket_value)[i][4];}
				}
				j = 0;
				iter = 0;
				//search close buckets for similair nodes if none are found in initial bucket..
				while( min > 0.5 && iter < 20){
					j++;
					iter++;
					if((bucket_value -j) > 0){
						for(i = 0; i < buckets.at(bucket_value - j).size(); i++){
							node = node_set.findNodeByIdentifier(buckets.at(bucket_value-j)[i][0]);
							cache.setNode(node);
							field.evaluateReal(cache, 3, out);
							//to_prolate_spheroidal(x, xps);
							d = distance(out, xps);
							
							if(d < min){
								node_index = node.getIdentifier(); 
								min = d;
								//std::cerr << "distance is "<< " " << d << std::endl;
								ele_index =  buckets.at(bucket_value-j)[i][1];
								local_coords[0] = buckets.at(bucket_value-j)[i][2];
								local_coords[1] = buckets.at(bucket_value-j)[i][3];
								local_coords[2] = buckets.at(bucket_value-j)[i][4];
							}
						}
					}
					if((bucket_value +j) < number_of_buckets){
						for(i = 0; i < buckets.at(bucket_value+j).size(); i++){
							node = node_set.findNodeByIdentifier(buckets.at(bucket_value+j)[i][0]);
							cache.setNode(node);
							field.evaluateReal(cache, 3, out);
							//to_prolate_spheroidal(x, xps);
							d = distance(out, xps);
							
							if(d < min){
								node_index = node.getIdentifier(); 
								min = d;
								//std::cerr << "distance is "<< " " << d << std::endl;
								ele_index =  buckets.at(bucket_value+j)[i][1];
								local_coords[0] = buckets.at(bucket_value+j)[i][2];
								local_coords[1] = buckets.at(bucket_value+j)[i][3];
								local_coords[2] = buckets.at(bucket_value+j)[i][4];
							}
						}
					}
				}
		/*if(buckets.at(bucket_value).size() != 0){
			for(i = 0; i < buckets.at(bucket_value).size(); i++){
				node = node_set.findNodeByIdentifier(buckets.at(bucket_value)[i][0]);
				cache.setNode(node);
				field.evaluateReal(cache, 3, out);
				//to_prolate_spheroidal(x, xps);
				d = distance(out, xps);
				
				if(d < min){
					node_index = node.getIdentifier(); 
					min = d;
					//std::cerr << "distance is "<< " " << d << std::endl;
					ele_index =  buckets.at(bucket_value)[i][1];
					local_coords[0] = buckets.at(bucket_value)[i][2];
					local_coords[1] = buckets.at(bucket_value)[i][3];
					local_coords[2] = buckets.at(bucket_value)[i][4];}
				}
			}

		//PROBABLY CHECK FOR EDGE CASES OF THE MAP OTHERWISE THERE MIGHT BE SEGMENTATION FAULT!
		else{
			int j;
			while((buckets.at(bucket_value-j).size() == 0) && buckets.at(bucket_value+j).size() == 0){
				j++;
			}
			
			for(i = 0; i < buckets.at(bucket_value - j).size(); i++){
				node = node_set.findNodeByIdentifier(buckets.at(bucket_value-j)[i][0]);
				cache.setNode(node);
				field.evaluateReal(cache, 3, out);
				//to_prolate_spheroidal(x, xps);
				d = distance(out, xps);
				
				if(d < min){
					node_index = node.getIdentifier(); 
					min = d;
					//std::cerr << "distance is "<< " " << d << std::endl;
					ele_index =  buckets.at(bucket_value-j)[i][1];
					local_coords[0] = buckets.at(bucket_value-j)[i][2];
					local_coords[1] = buckets.at(bucket_value-j)[i][3];
					local_coords[2] = buckets.at(bucket_value-j)[i][4];}
				}
			for(i = 0; i < buckets.at(bucket_value+j).size(); i++){
				node = node_set.findNodeByIdentifier(buckets.at(bucket_value+j)[i][0]);
				cache.setNode(node);
				field.evaluateReal(cache, 3, out);
				//to_prolate_spheroidal(x, xps);
				d = distance(out, xps);
				if(d < min){
					node_index = node.getIdentifier(); 
					min = d;
					//std::cerr << "distance is "<< " " << d << std::endl;
					ele_index =  buckets.at(bucket_value+j)[i][1];
					local_coords[0] = buckets.at(bucket_value+j)[i][2];
					local_coords[1] = buckets.at(bucket_value+j)[i][3];
					local_coords[2] = buckets.at(bucket_value+j)[i][4];
					}
				}
			}
		*/
		element = mesh.findElementByIdentifier(ele_index);
		node = node_set.findNodeByIdentifier(node_index);
		cache.setNode(node);
		field.evaluateReal(cache,3,x1);
		//std::cout << "local node coordinates" << local_coords[0] << " "  << local_coords[1] << " " << local_coords[2] << std::endl;
		//std::cout << "coordinates of found node: "<< x1[0] << " " << x1[1] << " " << x1[2] << std::endl;
		cache.setMeshLocation(element,3 , local_coords);
		field_fibers.evaluateReal(cache, 3, out);
		field.evaluateReal(cache,3,x1);
		//std::cout << "coordinates of found node with local coords: "<< x1[0] << " " << x1[1] << " " << x1[2] << std::endl;
		double fiber_angle = out[0];
		//std::cerr<< "fiber angle of found node: " << fiber_angle << std::endl;
		double slope = std::tan(fiber_angle);
		d = 0.0001*std::cos(fiber_angle);
		local_coords[0] +=  d;
		local_coords[1] += slope*d;
		
		//std::cerr << ele_index << std::endl;
		cache.setMeshLocation(element, 3, local_coords);
		field.evaluateReal(cache, 3, xps);
		//std::cout << "coordinates of new node: "<< xps[0] << " " << xps[1] << " " << xps[2] << "with distance " << distance(x1, xps) << std::endl;
		to_euclidian(xps, out);
		/*std::cerr << "found vector is: " << out[0] << " " << out[1] << " "  << out[2] << std::endl;
		std::cerr << "distance is: " << min << std::endl;
		std::cerr << "steps taken:" << iter << std::endl;
		std::cerr << "_______________________________________________" << std::endl;*/
		xps[0] = x[0] - out[0];
		xps[1] = x[1] - out[1];
		xps[2] = x[2] - out[2];
		double u = 1.0/(std::sqrt(xps[0]*xps[0] + xps[1]*xps[1] + xps[2]*xps[2]));
		out[0] = u*xps[0];
		out[1] = u*xps[1];
		out[2] = u*xps[2];

		//std::cout<< "final fiber vector: "<< out[0] << " " << out[1] << " " << out[2] << std::endl;

		//get_fiber_values(coords, precision, out_values, fieldmodule);
		for(int i = 0; i < dim; i++){
			nu[i] = out[i];
		}
		//nu = 0.0;
		//nu[1] = 1.0;
		if(nu[0] != nu[0]){std::cerr << "Koordinaten: "<<  x_vec_in << std::endl;std::cerr << nu << std::endl;}
		//std::cerr << out_values[0] <<  " "<< out_values[1] << " " << out_values[2] << std::endl;
		//std::cerr << nu << std::endl; 


	}
	

	std::string config_string() const;

private:
	

    MathVector<dim> m_sigma;

  



};


} // end OpenCMISS
} // end ug


#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__LOGNORMAL_FIELD__ */