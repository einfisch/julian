#object for parsing a VTU-file. writes the point coordinates, connectivity, offset, types, regions in seperate .txt-files
import vtu_write
class VTU_converter:
	
	def parse (vtu_file):
		pass
	def heart_to_vtu (self, con, coo, out = "heart_mesh.vtu"):

		conFile = open(con, "r")
		coordFile = open(coo, "r")

		num = vtu_write.wideMeshWedges(conFile, coordFile) #wideMeshWedges returns number of vertices, cells and subsets in a list
		numPoints = num[0]
		numCells = num[1] + num[0]
		subsets = num[2]
		vtu_write.create_file(coo, numPoints, numCells, subsets, out, False)


	def write_directional_field (self, coo, out = "fibres.vtu", r=1):
		coordFile = open(coo, "r")
		
		num = vtu_write.directional_field(coordFile, r)
		numPoints = num[0]
		numCells = num[1] + num[0]
		subsets = num[2]
		vtu_write.create_file("vertices_angles.txt", numPoints, numCells, subsets, out, True)
#converter = VTU_converter()

"""coords = "nodes_euclidian.txt"
cons = "connectivity.txt"
#converter.heart_to_vtu_fine(cons, coords)
converter.heart_to_vtu(cons, coords, "test_1connectivity.vtu")"""
