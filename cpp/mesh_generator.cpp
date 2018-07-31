#include "vtu_converter.cpp"
#include "exfile_reader.cpp"
//#include "vtu_write.cpp"
#include "prolate_spheroidal_to_euclidian.cpp"
#include <string>
class mesh_generator{
	const char* exfile = "heart.exfile";std::string out_file = "heart_mesh.vtu", out_file_fibers = "fibers.vtu";
	 int ref1 = 0, ref2 = 0, ref3 = 0, ref1_fibers = 0, ref2_fibers = 0, ref3_fibers = 0, 
	 number_of_buckets = 10000;
	bool surface = false; double fiber_length = 3;
public:
	void set_fiber_length(double);
	void set_fibers_surface(bool);
	void set_refinement(int,int,int);
	void set_refinement_fibers(int,int,int);
	void set_exfile(const char*);
	void set_out_file(std::string);
	void set_out_file_fibers(std::string);
	void set_number_of_buckets(int);
	void generate_mesh(bool);

};

void mesh_generator::set_fiber_length(double d){
	fiber_length = d;
}
void mesh_generator::set_fibers_surface(bool srfc){
	surface = srfc;
}
void mesh_generator::set_refinement(int x, int y, int z){
	ref1 = x;
	ref2 = y;
	ref3 = z;
}
void mesh_generator::set_refinement_fibers(int x, int y, int z){
	ref1_fibers = x;
	ref2_fibers = y;
	ref3_fibers = z;
}
void mesh_generator::set_exfile(const char* str){
	exfile = str;
}
void mesh_generator::set_out_file(std::string str){
	out_file = str;
}
void mesh_generator::set_number_of_buckets(int buckets){
	number_of_buckets = buckets;
}
void mesh_generator::set_out_file_fibers(std::string str){
	out_file_fibers = str;
}
void mesh_generator::generate_mesh(bool delete_data){

	exfile_reader ex_reader;
	ex_reader.set_file(exfile);
	ex_reader.set_number_of_buckets(number_of_buckets);
	vtu_converter vtu_conv;

	ex_reader.create_interpolation("interpolation.exfile", ref1, ref2, ref3);
	ex_reader.create_directional_field("fibers.exfile", ref1_fibers, ref2_fibers, ref3_fibers, surface);
	exfile_reader ex_reader_fibers;
	ex_reader_fibers.set_number_of_buckets(number_of_buckets);
	ex_reader_fibers.set_file("fibers.exfile");
	ex_reader.set_file("interpolation.exfile");

	std::cout << "converting exfile..." << std::endl;

	ex_reader_fibers.create_node_coordinates("coordinates_fibers_prolate_spheroidal");
	ex_reader.create_node_coordinates("coordinates_prolate_spheroidal");
	ex_reader.create_node_connectivity("connectivity");
	ex_reader.create_fiber_values("fiber_values");
	std::cout << "checking..." << std::endl;
	prolate_spheroidal_to_euclidian("coordinates_prolate_spheroidal", "coordinates", 35.25);
	
	prolate_spheroidal_to_euclidian("coordinates_fibers_prolate_spheroidal", "coordinates_fibers", 35.25);

	vtu_conv.fibers_to_vtu("coordinates_fibers", out_file_fibers, fiber_length);
	std::cout << "works" << std::endl;
	vtu_conv.heart_to_vtu("connectivity", "coordinates", out_file);



}
int main(){

	mesh_generator mesh_gen;
	mesh_gen.set_exfile("heart1.exfile");
	mesh_gen.set_out_file("heart222.vtu");
	mesh_gen.set_number_of_buckets(10000);
	mesh_gen.set_out_file_fibers("fibers220.vtu");
	mesh_gen.set_refinement_fibers(2,2,0);
	mesh_gen.set_refinement(2,2,2);
	mesh_gen.generate_mesh(true);

	return 0;

}