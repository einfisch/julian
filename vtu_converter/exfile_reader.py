import os
import math
from opencmiss.zinc.context import Context
from opencmiss.zinc.status import OK as ZINC_OK
def distance(coord1, coord2):
	focus = 0.3525E+02
	x1 = focus * math.cosh(coord1[0]) * math.cos(coord1[1])
	y1 = focus * math.sinh(coord1[0]) * math.sin(coord1[1]) * math.cos(coord1[2])
	z1 = focus * math.sinh(coord1[0]) * math.sin(coord1[1]) * math.sin(coord1[2])
	x2 = focus * math.cosh(coord2[0]) * math.cos(coord2[1])
	y2 = focus * math.sinh(coord2[0]) * math.sin(coord2[1]) * math.cos(coord2[2])
	z2 = focus * math.sinh(coord2[0]) * math.sin(coord2[1]) * math.sin(coord2[2])

	d = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
	return(d)
def to_euclidian(coords):
	focus = 0.3525E+02
	x = focus * math.cosh(coords[0]) * math.cos(coords[1])
	y = focus * math.sinh(coords[0]) * math.sin(coords[1]) * math.cos(coords[2])
	z = focus * math.sinh(coords[0]) * math.sin(coords[1]) * math.sin(coords[2])
	return([x,y,z])
class exfile_reader:
	def __init__(self, in_file):
		self.in_file = in_file
	def set_file(in_file):
		self.in_fille = in_file
	def create_node_coordinates (self, out_file = "coordinates_prolate_spheroidal.txt"):
		"""takes an .exfile document, extracts the prolate spheroidal coordinates at each nodes and write them in
		order into a new file. The function works by creating a node iterator from the nodeset wich contains all
		nodes """
		if os.path.isfile("/tmp/" + out_file):
			os.remove("/tmp/" + out_file)
		out = open(out_file, "w+")
		context = Context("heart")
		region = context.getDefaultRegion()
		region.readFile(self.in_file)
		fieldmodule = region.getFieldmodule()
		field = fieldmodule.findFieldByName("coordinates")
		field = field.castFiniteElement()
		#test = fieldmodule.createFieldDerivative(field, 2)
		cache = fieldmodule.createFieldcache()
		node_set = fieldmodule.findNodesetByName("nodes")
		node_iter = node_set.createNodeiterator()	#create the node iterator
		counter = 0
		node = node_iter.next()
		print (field)
		while node.isValid():
			cache.setNode(node) 	#sets the fieldcache-position to the current node
			"""test = field.getNodeParameters(cache, 1, node.VALUE_LABEL_D_DS1, 1, 10)
			print (test)"""
			result, out_values = field.evaluateReal(cache, 3)#evaluates real coordinates at the current position of cache
			#test1 = test.evaluateReal(cache, 2)
			#print(test.)
			if result == ZINC_OK:
				out.write(str(out_values[0]) + " " + str(out_values[1]) + " " +str(out_values[2]) + "\r\n")
			else:
				break
			node = node_iter.next()
			counter += 1
		print(str(counter) + " nodes evaluated successfully")
	def create_node_connectivity(self, out_file = "connectivity.txt"):
		if os.path.isfile("/tmp/" + out_file):
			os.remove("/tmp/" + out_file)
		out = open(out_file, "w+")
		context = Context("heart")
		region = context.getDefaultRegion()
		region.readFile(self.in_file)
		fieldmodule = region.getFieldmodule()
		field = fieldmodule.findFieldByName("coordinates")
		cache = fieldmodule.createFieldcache()
		mesh = fieldmodule.findMeshByDimension(3)
		ele_iter = mesh.createElementiterator()
		#d_ds1 = mesh.getChartDifferentialoperator(1, 2) 
		#print(d_ds1.isValid())
		counter = 0
		element = ele_iter.next()
		ele_template = element.getElementfieldtemplate(field, 1)
		
		while element.isValid():
			for i in range(1,9):
				
				current_node = element.getNode(ele_template, i)
				cache.setElement(element)
				#test = field.evaluateDerivative(d_ds1, cache, 1)
				#print(test)
				local_index = current_node.getIdentifier()
				out.write(str(local_index) + " ")
			out.write("\r\n")
			element = ele_iter.next()
			ele_template = element.getElementfieldtemplate(field, 1)



	def create_fibre_values (self, out_file = "fibre_values_spheroidal.txt"):
		"""takes an .exfile document, extracts the prolate spheroidal coordinates at each nodes and write them in
		order into a new file. The function works by creating a node iterator from the nodeset wich contains all
		nodes """
		if os.path.isfile("/tmp/" + out_file):
			os.remove("/tmp/" + out_file)
		out = open(out_file, "w+")
		context = Context("heart")
		region = context.getDefaultRegion()
		region.readFile(self.in_file)
		fieldmodule = region.getFieldmodule()
		field = fieldmodule.findFieldByName("fibres")
		cache = fieldmodule.createFieldcache()
		node_set = fieldmodule.findNodesetByName("nodes")
		node_iter = node_set.createNodeiterator()	#create the node iterator
		counter = 0
		node = node_iter.next()

		while node.isValid():
			cache.setNode(node) 	#sets the fieldcache-position to the current node
			result, out_values = field.evaluateReal(cache, 3)#evaluates real coordinates at the current position of cache
			if result == ZINC_OK:
				out.write(str(out_values[0]) + " " + str(out_values[1]) + " " +str(out_values[2]) + "\r\n")
			else:
				break
			node = node_iter.next()
			counter += 1
		print(str(counter) + " nodes evaluated successfully")
	

	def test_interpolation (self, out_file = "interpolation.exfile", refinement1 = 1, refinement2 = 1, refinement3 = 1):
		
		if os.path.isfile("/tmp/" + out_file):
			os.remove("/tmp/" + out_file)
		out = open(out_file, "w+")
		context = Context("heart")
		region = context.getDefaultRegion()
		region.readFile(self.in_file)
			
		fieldmodule = region.getFieldmodule()
		field = fieldmodule.findFieldByName("coordinates")
		field_fibres = fieldmodule.findFieldByName("fibres")

		mesh = fieldmodule.findMeshByDimension(3)
		

		field = field.castFiniteElement()
		field_fibres = field_fibres.castFiniteElement()
		cache = fieldmodule.createFieldcache()
		node_set = fieldmodule.findNodesetByName("nodes")
		node_template = node_set.createNodetemplate()
		node_template.defineField(field)
		node_template.defineField(field_fibres)
		node_template.setValueNumberOfVersions(field,  1, 2, 1)#define, that d/ds1, d/ds2, d2/ds1ds2 should be saved in the node lamda
		node_template.setValueNumberOfVersions(field,  1, 3, 1)#mu
		node_template.setValueNumberOfVersions(field,  1, 4, 1)#theta
		node_template.setValueNumberOfVersions(field_fibres,  1, 2, 1)#fibre angle
		node_template.setValueNumberOfVersions(field_fibres,  3, 2, 1)#sheet angle
		node_template.setValueNumberOfVersions(field_fibres,  3, 3, 1)
		node_template.setValueNumberOfVersions(field_fibres,  3, 4, 1)
		node_iter = node_set.createNodeiterator()	#create the node iterator

		l2 = 1.0/2**refinement2
		l3 = 1.0/2**refinement3
		l1 = 1.0/2**refinement1	#length of sides of local cells after refinement
		
		counter = 0
		
		ele_counter = 0
		num_ele = mesh.getSize()
		mesh_size = mesh.getSize()

		buckets = {}
		bucket_values = []
		#node_coordinates = []#list that contains all node coordinates to check for doubles.
		node1 = node_iter.next()
		while node1.isValid(): #get number of total nodes befor refinement
			counter += 1
			node1 = node_iter.next()
			cache.setNode(node1)
			result, coords = field.evaluateReal(cache, 3)
			coords = to_euclidian(coords)
			coord_value = coords[0] + coords[1] + coords[2]
			buckets[coord_value] = [node1]
			bucket_values.append(coord_value)
			#node_coordinates.append(coords)
		#print(buckets)
		bucket_values = sorted(bucket_values)
		print(bucket_values)
		#print (bucket_values)
		c = 0
		fieldmodule.beginChange()
		while c < mesh_size:
			element = mesh.findElementByIdentifier(c +1)
			#get elementfieldtemplates to create new elements of the refined mesh
			ele_field_template1 = element.getElementfieldtemplate(field,1)
			
			ele_field_template_lamda = element.getElementfieldtemplate(field, 1)
			ele_field_template_mu = element.getElementfieldtemplate(field, 2)
			ele_field_template_theta = element.getElementfieldtemplate(field, 3)
			ele_field_template_fibre = element.getElementfieldtemplate(field_fibres, 1)
			ele_field_template_imbrication = element.getElementfieldtemplate(field_fibres, 2)
			ele_field_template_sheet = element.getElementfieldtemplate(field_fibres, 3)
			#iterte over x,y and z axis, depending on refinement for x,y and z axis
			node_set_test = []
			for i in range(0, ele_field_template1.getNumberOfLocalNodes() ):
				current_node = element.getNode(ele_field_template1, i + 1)
				node_set_test.append(current_node)
			#print (node_set_test)
			for i in range(0, 2**refinement1):
				for j in range(0, 2**refinement2 ):
					for k in range(0, 2**refinement3 ):

						node_identifiers = [] #list of node identifiers to be used for the new element
						#local coordinates of one new element, based on local coordinates of parent element
						local_node_coordinates = [[0 + i*l1,0 + j*l2,0 + k*l3], [l1 + i*l1 ,0 + j*l2,0 + k*l3],\
						 [0 + i*l1,l2 + j*l2,0 + k*l3],[l1 + i*l1,l2 + j*l2,0 + k*l3], [0 + i*l1,0 + j*l2,l3 + k*l3],\
						 [l1 + i*l1,0 + j*l2,l3 + k*l3],\
						  [0 + i*l1 ,l2 + j*l2,l3 + k*l3],[l1 + i*l1 ,l2 + j*l2,l3 + k*l3]]	
						"""node_coordinates = [[0 + i*le,0 + j*le,0 + k*le], [le + i*le ,0 + j*le,0 + k*le],\
						 [0 + i*le,le + j*le,0 + k*le],[le + i*le,le + j*le,0 + k*le], [0 + i*le,0 + j*le,le + k*le],\
						 [le + i*le,0 + j*le,le + k*le],\
						  [0 + i*le ,le + j*le,le + k*le],[le + i*le ,le + j*le,le + k*le]]	"""


						#print(local_node_coordinates)
						#create a new node for every local coordinate
						
						for local_coords in local_node_coordinates:
							#interpolate global coordinates of parent element to get global coordinates of new element by using the local coordinates

							cache.setMeshLocation(element, local_coords)
							result, global_coords = field.evaluateReal(cache, 3)
							euclidian_coords = to_euclidian(global_coords)
							coord_value = euclidian_coords[0] + euclidian_coords[1] + euclidian_coords[2]
							#do the same for for the fibre-angles
							result, global_angles = field_fibres.evaluateReal(cache, 3)
							#print(global_coords)
							bucket = 0
							#print(coord_value)
							while (bucket < len(bucket_values) -1):
								if coord_value > bucket_values[bucket]:
									bucket += 1
								else:
									break
							
							#print(bucket)
							#check if coordinates belong to an already existing node by comparing euclidian distance
							already_existing = False

							unified_buckets = []
							if bucket == 0:
								unified_buckets = buckets[bucket_values[bucket]] + buckets[bucket_values[bucket +1]]
							elif bucket == 98:
								unified_buckets = buckets[bucket_values[bucket -1]] + buckets[bucket_values[bucket]]
							else:
								unified_buckets = buckets[bucket_values[bucket -1]] + buckets[bucket_values[bucket]] + buckets[bucket_values[bucket +1]]
							for node_test in unified_buckets:
								
								#print(distance(global_coords,coordinates))
								cache.clearLocation()
								cache.setNode(node_test)
								coord_test = []
								result, coord_test = field.evaluateReal(cache, 3)
								if (distance(global_coords,coord_test) < 1E-5):
									already_existing = True
									#print("TRIGGERED")
									node_identifiers.append(node_test.getIdentifier())
									
									break
							


							#create new node with the interpolated global coordinates and interpolated fibre-angles
							if (already_existing == False):
								node = node_set.createNode(counter + 1, node_template)
								#print(node.isValid())
								#cache.clearLocation()
								cache.setNode(node)
								field.assignReal(cache, global_coords)
								field_fibres.assignReal(cache, global_angles)
								node_identifiers.append(counter + 1)
								buckets[bucket_values[bucket]].append(node)
								#print(len(buckets[bucket_values[bucket]]))

								counter += 1
								
								#node_set_test.append(node)
							#test = field.getNodeParameters(cache, 2, 2,1,1)
							
						#print(node_identifiers)

						#create element template and set basis-functions and shape
						Elementbasis = fieldmodule.createElementbasis(2, 1)
						element_template = mesh.createElementtemplate()
						#if first and second node identifier are the same -> wedge shape
						if(node_identifiers[0] == node_identifiers[1]):
							a = node_identifiers	
							#set the last two nodes to "undefined" by setting their node identifier as 0						
							node_identifiers = [a[0], a[2], a[3],a[4],a[6], a[7],0,0]
							

						
						node_indexes = [1,2,3,4,5,6,7,8]
						element_template.setElementShapeType(element.SHAPE_TYPE_CUBE)
						element_template.setNumberOfNodes(8)
							
						basislamda = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)#cubic hermite
						basismu = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)#lagrange
						basistheta = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)#decreasing in x1 fehlt
						
						
						

						#element_template.defineField(field,-1, ele_field_template_lamda)
						#element_template.defineField(field,2, ele_field_template_mu)
						#element_template.defineField(field,3, ele_field_template_theta)
						#element_template.defineField(field_fibres,1, ele_field_template_fibre)
						#element_template.defineField(field_fibres,2, ele_field_template_imbrication)
						#element_template.defineField(field_fibres,3, ele_field_template_sheet)

						element_template.defineFieldSimpleNodal(field, 1, basislamda, node_indexes)
						#
						element_template.defineFieldSimpleNodal(field, 2, basismu, node_indexes)
						element_template.defineFieldSimpleNodal(field, 3, basistheta, node_indexes)
						#create new element, based on element-template using the new nodes.
						z = 1
						#print(node_identifiers)
						for identifier in node_identifiers:
							node = node_set.findNodeByIdentifier(identifier)
							
							element_template.setNode(z, node)
							z += 1
						
						
						new_element = mesh.createElement(num_ele +1,element_template )
						
						num_ele += 1
			print(str(num_ele) + " of " + str(60 * 2**refinement1 * 2**refinement2 * 2**refinement3 -60) + " elements created.")
			mesh.destroyElement(element)
			c += 1

		region.writeFile(out_file)
		print("done")
		fieldmodule.endChange()		
	
"""if ([round(elem, 6) for elem in coordinates] in finished):
									x3 += 0.5
									
									continue"""
#reader = exfile_reader("plswor.exfile")
#reader.create_node_coordinates()
#reader = exfile_reader("heart.exfile")
#reader.test_interpolation()
#*
"""while x1 <= 1:
				while x2 <= 1:
					while x3 <= 1:
						
						xi = [x1, x2, x3]
						
						
						if (0.5 not in xi):
							x3 += 0.5
							continue

						else:
							print(xi)

							cache.setMeshLocation(element, xi)
							result, coordinates = field.evaluateReal(cache,3)
							node = node_set.createNode(counter + 1, node_template)
							if ([round(elem, 6) for elem in coordinates] in finished):
								x3 += 0.5
								print("lol")
								continue
							cache.setNode(node)
							finished.append([round(elem, 6) for elem in coordinates])
							field.assignReal(cache, coordinates)
							counter += 1
							out.write(str(coordinates[0]) + " " + str(coordinates[1] )+ " " + str(coordinates[2]) +'\r\n' )

						x3 += 0.5
					x2 += 0.5
					x3 = 0
				x1 += 0.5
				x2 = 0
				x3 = 0
			x1 = 0
			x2 = 0
			x3 = 0
			element = el_iter.next()
			ele_field_template = element.getElementfieldtemplate(field,3)
			basislamda = fieldmodule.createElementbasis(3, 7)#cubic hermite
			basismu = fieldmodule.createElementbasis(3, 2)#lagrange
			basistheta = fieldmodule.createElementbasis(3, 2)#decreasing in x1 fehlt
			
			for i in range(2,4):

				for j in range(2,4):
					for k in range(2,4):
						num_ele += 1
						ele_template = mesh.createElementtemplate()
						#new_element = mesh.createElement(num_ele, ele_template)
						
						ele_template.setElementShapeType(element.SHAPE_TYPE_CUBE)
						ele_template.setNumberOfNodes(8)
						#ele_template.defineField(field,1, ele_field_template)

						local_index = 1
						node_indexes = []
						
						for x1 in range(i-2,i):	#first new element
							for x2 in range(j-2,j):
								for x3 in range(k-2,k):
									xi = [x1/2.0, x2/2.0, x3/2.0]
									#print(xi)			
									node_indexes.append(local_index)		
									cache.setMeshLocation(element, xi)
									result, coordinates = field.evaluateReal(cache, 3)
									
									node = node_set.createNode(counter + 1, node_template)
									cache.setNode(node)
									
									#print(coordinates)
									field.assignReal(cache, coordinates)

									parameters = field.getNodeParameters(cache, 1,1,1,1 )
									#print(parameters)
									counter += 1
									local_index += 1
									ele_template.setNode(local_index, node)
										
						mesh.defineElement(num_ele, ele_template)
						ele_template.defineFieldSimpleNodal(field, 1, basislamda, range(1,9))
						ele_template.defineFieldSimpleNodal(field, 2, basismu, range(1,9))
						ele_template.defineFieldSimpleNodal(field, 3, basistheta, range(1,9))
			"""

"""def test_interpolation (self, out_file = "interpolation.exfile", refinement1 = 1, refinement2 = 1, refinement3 = 1):
		
		if os.path.isfile("/tmp/" + out_file):
			os.remove("/tmp/" + out_file)
		out = open(out_file, "w+")
		context = Context("heart")
		region = context.getDefaultRegion()
		region.readFile(self.in_file)
			
		fieldmodule = region.getFieldmodule()
		field = fieldmodule.findFieldByName("coordinates")
		field_fibres = fieldmodule.findFieldByName("fibres")

		mesh = fieldmodule.findMeshByDimension(3)
		

		field = field.castFiniteElement()
		field_fibres = field_fibres.castFiniteElement()
		cache = fieldmodule.createFieldcache()
		node_set = fieldmodule.findNodesetByName("nodes")
		node_template = node_set.createNodetemplate()
		node_template.defineField(field)
		node_template.defineField(field_fibres)
		node_template.setValueNumberOfVersions(field,  1, 2, 1)#define, that d/ds1, d/ds2, d2/ds1ds2 should be saved in the node lamda
		node_template.setValueNumberOfVersions(field,  1, 3, 1)#mu
		node_template.setValueNumberOfVersions(field,  1, 4, 1)#theta
		node_template.setValueNumberOfVersions(field_fibres,  1, 2, 1)#fibre angle
		node_template.setValueNumberOfVersions(field_fibres,  3, 2, 1)#sheet angle
		node_template.setValueNumberOfVersions(field_fibres,  3, 3, 1)
		node_template.setValueNumberOfVersions(field_fibres,  3, 4, 1)

		node_iter = node_set.createNodeiterator()	#create the node iterator
		l2 = 1.0/2**refinement2
		l3 = 1.0/2**refinement3
		l1 = 1.0/2**refinement1	#length of sides of local cells after refinement
		
		counter = 0
		node1 = node_iter.next()
		ele_counter = 0
		num_ele = mesh.getSize()
		mesh_size = mesh.getSize()


		while node1.isValid():
			counter += 1
			node1 = node_iter.next()
		
		c = 0
		fieldmodule.beginChange()
		while c < mesh_size:
			element = mesh.findElementByIdentifier(c +1)
			ele_field_template_lamda = element.getElementfieldtemplate(field, 1)
			ele_field_template_mu = element.getElementfieldtemplate(field, 2)
			ele_field_template_theta = element.getElementfieldtemplate(field, 3)
			ele_field_template_fibre = element.getElementfieldtemplate(field_fibres, 1)
			ele_field_template_imbrication = element.getElementfieldtemplate(field_fibres, 2)
			ele_field_template_sheet = element.getElementfieldtemplate(field_fibres, 3)

			for i in range(0, 2**refinement1):
				for j in range(0, 2**refinement2 ):
					for k in range(0, 2**refinement3 ):



						node_coordinates = [[0 + i*l1,0 + j*l2,0 + k*l3], [l1 + i*l1 ,0 + j*l2,0 + k*l3],\
						 [0 + i*l1,l2 + j*l2,0 + k*l3],[l1 + i*l1,l2 + j*l2,0 + k*l3], [0 + i*l1,0 + j*l2,l3 + k*l3],\
						 [l1 + i*l1,0 + j*l2,l3 + k*l3],\
						  [0 + i*l1 ,l2 + j*l2,l3 + k*l3],[l1 + i*l1 ,l2 + j*l2,l3 + k*l3]]	
						node_coordinates = [[0 + i*le,0 + j*le,0 + k*le], [le + i*le ,0 + j*le,0 + k*le],\
						 [0 + i*le,le + j*le,0 + k*le],[le + i*le,le + j*le,0 + k*le], [0 + i*le,0 + j*le,le + k*le],\
						 [le + i*le,0 + j*le,le + k*le],\
						  [0 + i*le ,le + j*le,le + k*le],[le + i*le ,le + j*le,le + k*le]]	
						print(node_coordinates)
						for local_coords in node_coordinates:
							cache.setMeshLocation(element, local_coords)
							result, global_coords = field.evaluateReal(cache, 3)
							result, global_angles = field_fibres.evaluateReal(cache, 3)
							#print(global_coords)
							node = node_set.createNode(counter + 1, node_template)
							cache.setNode(node)
							field.assignReal(cache, global_coords)
							field_fibres.assignReal(cache, global_angles)

							#test = field.getNodeParameters(cache, 2, 2,1,1)
							
							counter += 1
						Elementbasis = fieldmodule.createElementbasis(2, 1)
						element_template = mesh.createElementtemplate()
						element_template.setElementShapeType(element.SHAPE_TYPE_CUBE)
						element_template.setNumberOfNodes(8)
						basislamda = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)#cubic hermite
						basismu = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)#lagrange
						basistheta = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)#decreasing in x1 fehlt
						

						
						node_indexes = [1,2,3,4,5,6,7,8]
						#element_template.defineField(field,-1, ele_field_template_lamda)
						#element_template.defineField(field,2, ele_field_template_mu)
						#element_template.defineField(field,3, ele_field_template_theta)
						#element_template.defineField(field_fibres,1, ele_field_template_fibre)
						#element_template.defineField(field_fibres,2, ele_field_template_imbrication)
						#element_template.defineField(field_fibres,3, ele_field_template_sheet)

						element_template.defineFieldSimpleNodal(field, 1, basislamda, node_indexes)
						#
						element_template.defineFieldSimpleNodal(field, 2, basismu, node_indexes)
						element_template.defineFieldSimpleNodal(field, 3, basistheta, node_indexes)
						
						for z in range(0,8):
							node = node_set.findNodeByIdentifier(counter -7 +z)
							element_template.setNode(z + 1, node)
						
						
						mesh.defineElement(num_ele +1,element_template )
						num_ele += 1
			print(num_ele)
			mesh.destroyElement(element)
			c += 1
			
			
				
			
		region.writeFile(out_file)
		print("done")
		fieldmodule.endChange()		"""