
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <typeinfo>
void split(std::string str, char delimiter, std::vector<std::string> &vec_out) {
	//seperate an input string at every delimiter char and write the reult in a given vector
  vec_out.clear();
  std::stringstream ss(str); // Turn the string into a stream.
  std::string token;
  //int i = 0;
  while(getline(ss, token, delimiter)) {
    vec_out.push_back(token);
    //std::cout << "token is of type "<< typeid(token).name() <<  " " << vec_out[i]<< std::endl;
    //i++;
  }
 
}
std::vector<int> heart_to_file(std::string connectivity, std::string coordinates){

	std::ofstream cell_type_file;
	std::ofstream region_file;
	std::ofstream offset_file;
	std::ofstream cell_file;
	std::ifstream connectivity_file;
	std::ifstream coordinates_file;

	cell_type_file.open("cell_type", std::ios::trunc);
	region_file.open("region", std::ios::trunc);
	offset_file.open("offset", std::ios::trunc);
	cell_file.open("cells", std::ios::trunc);
	connectivity_file.open(connectivity);
	coordinates_file.open(coordinates);

	int offset = 0;
	int region = 0;
	int num_points = 0;
	int num_cells = 0;

	std::string line;
	std::getline(coordinates_file, line);
	while(line.size() > 0){
		//count number of vertices to get offset for the vtu file format.
		offset += 1;
		offset_file << offset << '\n';
		cell_file << num_points << '\n';
		num_points += 1;
		cell_type_file << "1 \n"; //cell type 1 represents a single vertex in vtu file format
		region_file << region << '\n';
		line.clear();
		std::getline(coordinates_file, line); 
	}
	std::vector<std::string> s;
	std::getline(connectivity_file, line);
	while(line.size() > 0){
		//std::cout << line << std::endl;
		split(line, ' ', s);
		//std::cout << s.size() << std::endl;
		//std::cout << s[0] << s[1] << s[2] << s[3] << s[4] << s[5] << s[6] << s[7]<<  std::endl; 
		//if last two nodes in an element are undefined (element identifier == -1) we hvae an 
		//edge collapse and use wedge shape type.
		if(s[7] == "-1"){
			num_cells += 1;
			cell_file << (std::stoi(s[0]) - 1) << " " << (std::stoi(s[1]) - 1) << " " << (std::stoi(s[2]) - 1)
			 << " "<< (std::stoi(s[3]) - 1) << " " << (std::stoi(s[4]) - 1) << " " << (std::stoi(s[5]) - 1) << '\n';
			cell_type_file << "13 \n";
			offset += 6;
			offset_file << offset << '\n';
			region_file << region << '\n';
		}
		else{
			num_cells += 1;
			cell_file << (std::stoi(s[0]) - 1) << " " << (std::stoi(s[4]) - 1) << " " << (std::stoi(s[5]) - 1)
			<< " " << (std::stoi(s[1]) - 1) << " " << (std::stoi(s[2]) - 1) << " " <<(std::stoi(s[6]) - 1) << " "
			<< (std::stoi(s[7]) - 1) << " " << (std::stoi(s[3]) - 1) << '\n';
			cell_type_file << "12 \n";
			offset += 8;
			offset_file << offset << '\n';
			region_file << region << '\n';
		}
		line.clear();
		std::getline(connectivity_file, line);
	} 
	std::vector<int> result {num_points, num_cells, region};
	return result;
}

std::vector<int> fiber_directions(std::string coordinates, double d){
	//function takes euclidian coordinates and file from exfile_reader::create_directional_field
	//to create the edges to visualize the fiber directions.
	std::ofstream cell_type_file;
	std::ofstream region_file;
	std::ofstream offset_file;
	std::ofstream cell_file;
	std::ifstream coordinates_file;
	std::ofstream vertices;

	cell_type_file.open("cell_type_fibers", std::ios::trunc);
	region_file.open("region_fibers", std::ios::trunc);
	offset_file.open("offset_fibers", std::ios::trunc);
	cell_file.open("cells_fibers", std::ios::trunc);
	coordinates_file.open(coordinates);
	vertices.open("vertices_fibers", std::ios::trunc);

	int offset = 0;
	int region = 0;
	int num_cells = 0;
	int num_points = 0;
	bool toggle = false;

	std::string line;
	std::vector<std::string> coords;
	std::vector<std::string> old_coords(3,"");
	std::vector<std::vector<int>> connectivity;
	std::vector<int> edge(2,0);
	
	double x;
	double y;
	double z;
	double x2;
	double y2;
	double z2;
	double x3;
	double y3;
	double z3;
	double x1;
	double y1;
	double z1;
	double v1;
	double v2;
	double v3;
	double lam;
	double b;

	std::getline(coordinates_file, line);
	while(line.size() > 0){
		split(line, ' ', coords);

		cell_type_file << "1 \n";
		offset += 1;
		offset_file << offset << '\n';
		region_file << region << '\n';
		cell_file << num_points << '\n';

		if(toggle){
			edge[0] = num_points - 1;
			edge[1] = num_points;
			connectivity.push_back(edge);

			x = std::stof(old_coords[0]);
			y = std::stof(old_coords[1]);
			z = std::stof(old_coords[2]);

			x1 = std::stof(coords[0]);
			y1 = std::stof(coords[1]);
			z1 = std::stof(coords[2]);
			//calculate the direction vector of the straight line that approximates the 
			//fiber órientation
			v1 = x1 - x;
			v2 = y1 - y;
			v3 = z1 - z;
			//create two points on the straight line in the fiber direction. the vertex where the
			//fiber orientation is to be visualized must be in the middle of both points.
			//norm the distance between the two points to d from the input values
			b = v1*v1 + v2*v2 + v3*v3;
			lam = 0;
			if(!(b == 0)){
				lam = d/(std::sqrt(b));
			}
			x2 = x - 0.5*lam*v1;
			y2 = y - 0.5*lam*v2;
			z2 = y - 0.5*lam*v3;

			x3 = x - 0.5*lam*v1;
			y3 = y - 0.5*lam*v2;
			z3 = y - 0.5*lam*v3;

			vertices << x2 << " " << y2 << " " << z2 << '\n';
			vertices << x3 << " " << y3 << " " << z3 << '\n';

		}
		toggle = !toggle;
		num_points += 1;
		std::getline(coordinates_file, line);
		old_coords[0] = coords[0];
		old_coords[1] = coords[1];
		old_coords[2] = coords[2];

	}

	for(std::size_t i = 0; i < connectivity.size(); i++){
		num_cells += 1;
		cell_type_file << "3 \n"; //cell shape type of a straight line (edge) of vtu file format.
		offset += 2;
		offset_file << offset << '\n';
		cell_file << connectivity[i][0] << " " << connectivity[i][1] << '\n';
		region_file << "0 \n";
	}

	std::vector<int> result {num_points, num_cells, 0};
	return result;
}

void create_file(std::string coordinates, int num_points, int num_cells, int region,std::string out_file, bool fibers){
	//function that takes as arguments the file name of a coordinates file, the number of vertices,
	//cells and subsets(regions) as integers a boolean that describes if a fiber visualization is
	//to be created or not. if fibers == false, the heart mesh is to be created.
	//the result gets written to out_file as a vtu-file format. this can be read by promesh/paraview etc. 
	std::cout << "exporting " << out_file << "..." << std::endl;

	std::ofstream  vtu_out;

	vtu_out.open(out_file, std::ios::trunc);
	std::string line;
	//hardcode syntax for VTU-file format

	vtu_out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"> \n";
	vtu_out << "<UnstructuredGrid> \n";
	vtu_out << "<Piece NumberOfPoints=\"" << num_points<<"\" NumberOfCells=\"" << num_cells <<"\">\n";
	vtu_out << "<Points>\n";
	vtu_out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	//write vertex coordinates to VTU file:
	std::ifstream coordinates_file;
	coordinates_file.open(coordinates);
	std::getline(coordinates_file, line);
	while(!(line.size() == 0)){
		vtu_out << line << " \n";
		line.clear();
		std::getline(coordinates_file, line);
	}

	vtu_out << "</DataArray>\n";
	vtu_out << "</Points>\n";
	vtu_out << "<Cells>\n";
	vtu_out << " <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

	std::ifstream connectivity_file;
	if(fibers){connectivity_file.open("cells_fibers");}
	else{connectivity_file.open("cells");}

	line.clear();
	std::getline(connectivity_file, line);
	while(!(line.size() == 0)){
		vtu_out << line << " \n";
		line.clear();
		std::getline(connectivity_file, line);
	}
	vtu_out << "</DataArray>\n";
	vtu_out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

	std::ifstream offset_file;
	if(fibers){offset_file.open("offset_fibers");}
	else{offset_file.open("offset");}

	line.clear();
	std::getline(offset_file, line);
	while(!(line.size() == 0)){
		vtu_out << line << " \n";
		line.clear();
		std::getline(offset_file, line);
	}
	vtu_out << "</DataArray>\n";
	vtu_out << "<DataArray type=\"Int8\" Name=\"types\" format=\"ascii\">\n";

	std::ifstream cell_type_file;
	if(fibers){cell_type_file.open("cell_type_fibers");}
	else{cell_type_file.open("cell_type");}

	line.clear();
	std::getline(cell_type_file, line);
	while(!(line.size() == 0)){
		vtu_out << line << " \n";
		line.clear();
		std::getline(cell_type_file, line);
	}
	vtu_out << "</DataArray>\n";
	vtu_out << "</Cells>\n";
	vtu_out << "<PointData>\n";
	vtu_out << "</PointData>\n";
	vtu_out << "<CellData>\n";
	vtu_out << "<DataArray type=\"Int32\" Name=\"regions\" NumberOfComponents=\"1\" format=\"ascii\">\n";

	std::ifstream region_file;
	if(fibers){region_file.open("region_fibers");}
	else{region_file.open("region");}

	line.clear();
	std::getline(region_file, line);
	while(!(line.size() == 0)){
		vtu_out << line << " \n";
		line.clear();
		std::getline(region_file, line);
	}
	vtu_out << "</DataArray>\n";
	vtu_out << "</CellData>\n";
	vtu_out << "<RegionInfo Name=\"regions\">\n";
	vtu_out << "<Region Name=\"subset\"></Region>\n";
	for(int i = 1; i <= region; i++){
		vtu_out << "<Region Name=\"cell " << i << "\"></Region>\n";
	}
	vtu_out << "</RegionInfo>\n";
	vtu_out << "</Piece>";
	vtu_out << "</UnstructuredGrid> \n";
	vtu_out << "</VTKFile>";

}	
/*int main(){
	heart_to_file("test1.txt", "test.txt");
	return 0;
}*/