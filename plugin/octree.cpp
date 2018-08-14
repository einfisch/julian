#include <vector>
#include <iostream>
#include <queue>
#include <stack>
#include <memory>
#include <cmath>
struct ArrayOct {
	ArrayOct() = default;
	ArrayOct(size_t parent_id_) : parent_id(parent_id_) {}

	size_t parent_id;
	std::vector<std::vector<double>> values;
	std::vector<std::vector<double>> bbox;
	std::vector<size_t> children_ids;
	double x_max;
	double x_min;
	double y_max;
	double y_min;
	double z_max;
	double z_min;

};

//size_t lookup(const std::vector<ArrayOct>& tree, coords){}
//TODO: TEMPLATE THIS VALURE VECTOR
bool insert_leaf(std::vector<ArrayOct>& tree, std::vector<double> coords, std::vector<double> val){
	double x = coords[0];
	double y = coords[1];
	double z = coords[2];

	
	int node_id = 0;
	size_t level = 0;
	bool leaf = false;

	while(!leaf){
		if(tree[node_id].children_ids.size() == 0){leaf = true; break;}
		double x_node = (tree[node_id].x_max + tree[node_id].x_min)/2.0;
		double y_node = (tree[node_id].y_max + tree[node_id].y_min)/2.0;
		double z_node = (tree[node_id].z_max + tree[node_id].z_min)/2.0;
		//since z is switching fastest and x slowest, we can easily check for our new node:
		/*std::cout << "x: " << x_node << std::endl;
		std::cout << "y: "<< y_node << std::endl;
		std::cout << "z: "<< z_node << std::endl;*/

		/*for (size_t j = 0; j < tree[node_id].children_ids.size(); j++)
			std::cout << tree[node_id].children_ids[j] << std::endl;*/

		int child = 0;
		if(x > x_node){child += 4;}
		if(y > y_node){child += 2;}
		if(z > z_node){child += 1;}
		node_id = tree[node_id].children_ids[child]  ;
		/*std::cout << "new node id: "<< node_id << std::endl;
		std::cout << "lower left : " << tree[node_id].x_min << " " << tree[node_id].y_min << " " << tree[node_id].z_min << std::endl;
		std::cout << "upper right: " << tree[node_id].x_max << " " << tree[node_id].y_max << " " << tree[node_id].z_max << std::endl;*/
	} 
	tree[node_id].values.push_back(val);
	return true;
}
std::vector<double> lookup(const std::vector<ArrayOct>& tree, std::vector<double> coords){
	double x = coords[0];
	double y = coords[1];
	double z = coords[2];

	
	int node_id = 0;
	size_t level = 0;
	bool leaf = false;

	while(!leaf){
		if(tree[node_id].children_ids.size() == 0){leaf = true; break;}
		double x_node = (tree[node_id].x_max + tree[node_id].x_min)/2.0;
		double y_node = (tree[node_id].y_max + tree[node_id].y_min)/2.0;
		double z_node = (tree[node_id].z_max + tree[node_id].z_min)/2.0;
		//since z is switching fastest and x slowest, we can easily check for our new node:
		/*std::cout << "x: " << x_node << std::endl;
		std::cout << "y: "<< y_node << std::endl;
		std::cout << "z: "<< z_node << std::endl;*/

		/*for (size_t j = 0; j < tree[node_id].children_ids.size(); j++)
			std::cout << tree[node_id].children_ids[j] << std::endl;*/

		int child = 0;
		if(x > x_node){child += 4;}
		if(y > y_node){child += 2;}
		if(z > z_node){child += 1;}
		node_id = tree[node_id].children_ids[child]  ;
		/*std::cout << "new node id: "<< node_id << std::endl;
		std::cout << "lower left : " << tree[node_id].x_min << " " << tree[node_id].y_min << " " << tree[node_id].z_min << std::endl;
		std::cout << "upper right: " << tree[node_id].x_max << " " << tree[node_id].y_max << " " << tree[node_id].z_max << std::endl;*/
	} 
	double min = 1000;
	std::size_t val_index;
	if(tree[node_id].values.size() == 0){
		//std::cout << "no value found for " << x <<" " << y << " " << z << std::endl;
		std::vector<double> out = {0,0,0,1,1,1};
		return out;
	}
	for(std::size_t i = 0; i < tree[node_id].values.size(); i++){
		double xi = tree[node_id].values[i][0];
		double yi = tree[node_id].values[i][1];
		double zi = tree[node_id].values[i][2];
		double d = std::sqrt((xi-x)*(xi-x) + (yi-y)*(yi-y) +(zi-z)*(zi-z));
		if(d < min){min = d; val_index = i;}
	}
	/*std::cout << "size of bucket: " << tree[node_id].values.size() << std::endl;*/
	std::vector<double> out;
	out = tree[node_id].values[val_index];
	/*for(int i = 0; i < tree[node_id].values[val_index].size(); i++ ){
		std::cout << out[i] << " " ;
	}*/
	//std::cout << std::endl;
	return out;
}

std::vector<ArrayOct> generate_tree(size_t n_levels, std::vector<double> lower_left, std::vector<double> upper_rigth) {
	size_t parent_counter = 0;	
	size_t current_stack_head = 0;

	std::vector<ArrayOct> nodes;
	nodes.push_back({});
	nodes[0].x_min = lower_left[0];
	nodes[0].x_max = upper_rigth[0];
	nodes[0].y_min = lower_left[1];
	nodes[0].y_max = upper_rigth[1];
	nodes[0].z_min = lower_left[2];
	nodes[0].z_max = upper_rigth[2];

	current_stack_head++;

	// 8^0 + 8^1 + ... + 8^(n-1)
	size_t limit = 0;
	for (size_t i = 0; i < n_levels; i++)
		limit += static_cast<size_t>(std::pow(8, i));

	std::cout << "limit " << limit << std::endl;

	for (; parent_counter < limit; parent_counter++) {
		// hardcode
		double x_min = nodes[parent_counter].x_min;//nodes[current_stack_head].parent_id
		double x_max = nodes[parent_counter].x_max;
		double y_min = nodes[parent_counter].y_min;
		double y_max = nodes[parent_counter].y_max;
		double z_min = nodes[parent_counter].z_min;
		double z_max = nodes[parent_counter].z_max;

		double x = (x_max + x_min)/2.0;
		double y = (y_max + y_min)/2.0;
		double z = (z_max + z_min)/2.0;

		nodes[parent_counter].children_ids.push_back(current_stack_head );
			nodes.push_back({parent_counter});
		
		nodes[current_stack_head + 0].x_min = x_min;
		nodes[current_stack_head + 0].x_max = x;
		nodes[current_stack_head + 0].y_min = y_min;
		nodes[current_stack_head + 0].y_max = y;
		nodes[current_stack_head + 0].z_min = z_min;
		nodes[current_stack_head + 0].z_max = z;

		nodes[parent_counter].children_ids.push_back(current_stack_head + 1);
			nodes.push_back({parent_counter});

		nodes[current_stack_head + 1].x_min = x_min;
		nodes[current_stack_head + 1].x_max = x;
		nodes[current_stack_head + 1].y_min = y_min;
		nodes[current_stack_head + 1].y_max = y;
		nodes[current_stack_head + 1].z_min = z;
		nodes[current_stack_head + 1].z_max = z_max;

		nodes[parent_counter].children_ids.push_back(current_stack_head + 2);
			nodes.push_back({parent_counter});

		nodes[current_stack_head + 2].x_min = x_min;
		nodes[current_stack_head + 2].x_max = x;
		nodes[current_stack_head + 2].y_min = y;
		nodes[current_stack_head + 2].y_max = y_max;
		nodes[current_stack_head + 2].z_min = z_min;
		nodes[current_stack_head + 2].z_max = z;

		nodes[parent_counter].children_ids.push_back(current_stack_head + 3);
			nodes.push_back({parent_counter});

		nodes[current_stack_head + 3].x_min = x_min;
		nodes[current_stack_head + 3].x_max = x;
		nodes[current_stack_head + 3].y_min = y;
		nodes[current_stack_head + 3].y_max = y_max;
		nodes[current_stack_head + 3].z_min = z;
		nodes[current_stack_head + 3].z_max = z_max;

		nodes[parent_counter].children_ids.push_back(current_stack_head + 4);
			nodes.push_back({parent_counter});

		nodes[current_stack_head + 4].x_min = x;
		nodes[current_stack_head + 4].x_max = x_max;
		nodes[current_stack_head + 4].y_min = y_min;
		nodes[current_stack_head + 4].y_max = y;
		nodes[current_stack_head + 4].z_min = z_min;
		nodes[current_stack_head + 4].z_max = z;

		nodes[parent_counter].children_ids.push_back(current_stack_head + 5);
			nodes.push_back({parent_counter});

		nodes[current_stack_head + 5].x_min = x;
		nodes[current_stack_head + 5].x_max = x_max;
		nodes[current_stack_head + 5].y_min = y_min;
		nodes[current_stack_head + 5].y_max = y;
		nodes[current_stack_head + 5].z_min = z;
		nodes[current_stack_head + 5].z_max = z_max;

		nodes[parent_counter].children_ids.push_back(current_stack_head + 6);
			nodes.push_back({parent_counter});

		nodes[current_stack_head + 6].x_min = x;
		nodes[current_stack_head + 6].x_max = x_max;
		nodes[current_stack_head + 6].y_min = y;
		nodes[current_stack_head + 6].y_max = y_max;
		nodes[current_stack_head + 6].z_min = z_min;
		nodes[current_stack_head + 6].z_max = z;

		nodes[parent_counter].children_ids.push_back(current_stack_head + 7);
			nodes.push_back({parent_counter});

		nodes[current_stack_head + 7].x_min = x;
		nodes[current_stack_head + 7].x_max = x_max;
		nodes[current_stack_head + 7].y_min = y;
		nodes[current_stack_head + 7].y_max = y_max;
		nodes[current_stack_head + 7].z_min = z;
		nodes[current_stack_head + 7].z_max = z_max;
		/*for (size_t old_stack_head = current_stack_head; 
			 old_stack_head < current_stack_head + 8;
			 old_stack_head++) {
			nodes[parent_counter].children_ids.push_back(old_stack_head);
			nodes.push_back({parent_counter});
		}*/
		current_stack_head += 8;
	}

	return nodes;
}


