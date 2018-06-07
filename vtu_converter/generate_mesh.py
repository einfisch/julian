import vtu_converter
import exfile_reader
import os
from prolate_spheroidal_to_euclidian import prolate_spheroid_to_euclidian

class mesh_generator:
	def __init__(self, exfile = "heart.exfile", refinement1 = 1,refinement2 = 1, refinement3 = 1, out_file = "heart_mesh.vtu"):
		self.exfile = exfile
		self.refinement1 = refinement1
		self.refinement2 = refinement2
		self.refinement3 = refinement3	
		self.out_file = out_file
	def setRefinement(self, x):
		self.refinement = x
	def setExfile(self, exfile):
		self.exfile = exfile

	def generate_mesh(self, delete_data = True):
		ex_reader = exfile_reader.exfile_reader(self.exfile)
		vtu_conv = vtu_converter.VTU_converter()
		print(self.refinement1)
		if (self.refinement1 > 0 or self.refinement2 > 0 or self.refinement3 > 0):
			ex_reader.test_interpolation("interpolation.exfile", self.refinement1, self.refinement2, self.refinement3)
			ex_reader = exfile_reader.exfile_reader("interpolation.exfile")
		ex_reader.create_node_coordinates()
		ex_reader.create_node_connectivity()

		prolate_spheroid_to_euclidian("coordinates_prolate_spheroidal.txt", "coordinates_euclidian.txt", 0.3525E+02)
		coords = "coordinates_euclidian.txt"
		cons = "connectivity.txt"
		vtu_conv.heart_to_vtu(cons, coords, self.out_file)
		"""if (delete_data == True):
			directory = os.getcwd()
			os.remove(os.path.join(directory, "connectivity.txt"))
			os.remove(os.path.join(directory, "coordinates_euclidian.txt"))
			os.remove(os.path.join(directory, "coordinates_prolate_spheroidal.txt"))"""

mesh_generator = mesh_generator("heart1.exfile", 1,1,1, "heart_1.vtu")#
mesh_generator.generate_mesh()