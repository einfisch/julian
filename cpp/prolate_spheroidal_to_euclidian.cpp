#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
void split1(std::string str, char delimiter, std::vector<std::string> &vec_out) {
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
void prolate_spheroidal_to_euclidian(std::string in_file_string, std::string out_file_string, double focus){
	//takes a coordinates file in prolate spheroidal coordinates and creates a coordinate file
	//in euclididan coordinates.
	std::ofstream out_file;
	std::ifstream in_file;

	out_file.open(out_file_string, std::ios::trunc);
	in_file.open(in_file_string);

	std::string line;
	std::getline(in_file, line);
	std::vector<std::string> coords;
	double lamda;
	double mu;
	double theta;
	double x;
	double y;
	double z;

	while(!(line.size() == 0)){
		split1(line, ' ', coords);
		lamda = std::stof(coords[0]);
		mu = std::stof(coords[1]);
		theta = std::stof(coords[2]);
		x = focus * std::cosh(lamda) * std::cos(mu);
		y = focus * std::sinh(lamda) * std::sin(mu) * std::cos(theta);
		z = focus * std::sinh(lamda) * std::sin(mu) * std::sin(theta);
		out_file << x << " " << y << " " << z << '\n';
		line.clear();
		std::getline(in_file, line);
	}
}
/*int main(){
	return 0;
}*/