import vtu_converter
import exfile_reader
from vtu_write import directional_field
import os
from prolate_spheroidal_to_euclidian import prolate_spheroid_to_euclidian

class mesh_generator:
	def __init__(self, exfile = "heart.exfile", out_file = "heart_mesh.vtu", refinement1 = 1,refinement2 = 1, refinement3 = 1):
		self.exfile = exfile
		self.refinement1 = refinement1
		self.refinement2 = refinement2
		self.refinement3 = refinement3	
		self.refinement1_fibers = refinement1
		self.refinement2_fibers = refinement2
		self.refinement3_fibers = refinement3
		self.out_file = out_file
	def set_refinement(self, x,y,z):
		self.refinement1 = x
		self.refinement2 = y
		self.refinement3 = z
	def set_refinement_fibers(self, x,y,z):
		self.refinement1_fibers = x
		self.refinement2_fibers = y
		self.refinement3_fibers = z
	def set_exfile(self, exfile):
		self.exfile = exfile
	def set_outfile(self, out_file):
		self.out_file = out_file
	def generate_mesh(self, delete_data = True):
		ex_reader = exfile_reader.exfile_reader(self.exfile)
		vtu_conv = vtu_converter.VTU_converter()
		
		if (self.refinement1 > 0 or self.refinement2 > 0 or self.refinement3 > 0):
			ex_reader.kd_interpolation("interpolation.exfile", self.refinement1, self.refinement2, self.refinement3)
			ex_reader.create_directional_field("fibres.exfile", self.refinement1_fibers, self.refinement2_fibers, self.refinement3_fibers)
			ex_reader_fibres = exfile_reader.exfile_reader("fibres.exfile")
			#ex_reader.set_file("interpolation.exfile")
			ex_reader = exfile_reader.exfile_reader("interpolation.exfile")
			print("converting exfile...")
			ex_reader_fibres.create_node_coordinates("fibres_nodes.txt")
		else:
			ex_reader.create_directional_field("fibres.exfile", self.refinement1_fibers, self.refinement2_fibers, self.refinement3_fibers)
			ex_reader_fibres = exfile_reader.exfile_reader("fibres.exfile")
			print("converting exfile...")
			ex_reader_fibres.create_node_coordinates("fibres_nodes.txt")

		ex_reader.create_node_coordinates()
		ex_reader.create_node_connectivity()
		ex_reader.create_fibre_values()
		prolate_spheroid_to_euclidian("coordinates_prolate_spheroidal.txt", "coordinates_euclidian.txt", 0.3525E+02)
		prolate_spheroid_to_euclidian("fibres_nodes.txt", "fibres_nodes_euclidian.txt", 0.3525E+02)
		coords = "coordinates_euclidian.txt"
		cons = "connectivity.txt"
		vtu_conv.write_directional_field( "fibres_nodes_euclidian.txt","fibresx.vtu", 3)

		vtu_conv.heart_to_vtu(cons, coords, self.out_file)
		if (delete_data == True):
			directory = os.getcwd()
			os.remove(os.path.join(directory, "connectivity.txt"))
			os.remove(os.path.join(directory, "coordinates_euclidian.txt"))
			os.remove(os.path.join(directory, "coordinates_prolate_spheroidal.txt"))
		print("Done.")
mesh_generator = mesh_generator("heart1.exfile",  "heart_000.vtu")
mesh_generator.set_refinement_fibers(3,2,0)
mesh_generator.set_refinement(3,3,3)
mesh_generator.generate_mesh(True)
