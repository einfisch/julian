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

	def heart_to_vtu_fine (self, con, coo, out = "heart_mesh_fine.vtu"):
		"""uses coordinate and connectivity data-sets to construct a VTU-file. Uses 8 Hexahedron per 27 vertices and
		places them in a subset"""
		conFile = open(con, "r")
		coordFile = open(coo, "r")
		num = vtu_write.fineMeshWedges(conFile, coordFile) #wideMeshWedges returns number of vertices, cells and subsets in a list
		numPoints = num[0]
		numCells = num[1] + num[0]
		subsets = num[2]
		vtu_write.create_file(coo, numPoints, numCells, subsets, out)
	def directional_field (self, coo, angles, out = "directional_field.vtu"):
		coordFile = open(coo, "r")
		angleFile = open(angles, "r")
		num = vtu_write.directional_field(angleFile, coordFile)
		numPoints = num[0]
		numCells = num[1] + num[0]
		subsets = num[2]
		vtu_write.create_file("vertices_angles.txt", numPoints, numCells, subsets, out, True)
#converter = VTU_converter()

"""coords = "nodes_euclidian.txt"
cons = "connectivity.txt"
#converter.heart_to_vtu_fine(cons, coords)
converter.heart_to_vtu(cons, coords, "test_1connectivity.vtu")"""
