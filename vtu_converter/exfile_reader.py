import os
import math
import numpy as np
import time
import sys
from opencmiss.zinc.context import Context
from opencmiss.zinc.status import OK as ZINC_OK
def update_progress(progress, program_status):
	#takes a progress as float and prints a progressbar to the console
	bar = 20
	
	output = ""
	if (progress >= 1):
		program_status = "Finished...                                              \r\n"
	block_count = int(round(bar*progress))
	text = "\rPercent: [{0}] {1}% {2}".format( "#"*block_count + "-"*(bar-block_count), round(progress*100), program_status)
	sys.stdout.write(text)
	sys.stdout.flush()
	
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

def to_prolate_spheroidal(coords):
	focus = 0.3525E+02
	x = coords[0]
	y = coords[1]
	z = coords[2]

	theta = math.atan(z/y)
	mu = math.acos(1/(2*focus)*((y*y + z*z + (x + focus)**2)**0.5 - (y*y + z*z + (x - focus)**2)**0.5))
	lam = math.acosh(1/(2*focus)*((y*y + z*z + (x + focus)**2)**0.5 + (y*y + z*z + (x - focus)**2)**0.5))

	return([lam, mu, theta])
class exfile_reader:
	def __init__(self, in_file):
		self.in_file = in_file
	def set_file(in_file):
		self.in_file = in_file


	def get_fiber_values(self, coords, precision):
		#returns the fiber vector at the specified euclidian coordinates (list like [x,y,z])
		ps_coords = to_prolate_spheroidal(coords)#convert to prolate spheroidal coordinates
		print(ps_coords)
		
		context = Context("heart")
		region = context.getDefaultRegion()
		region.readFile(self.in_file)
		fieldmodule = region.getFieldmodule()
		field = fieldmodule.findFieldByName("coordinates")
		field_fibres = fieldmodule.findFieldByName("fibres")
		mesh = fieldmodule.findMeshByDimension(3)
		#field = field.castFiniteElement()
		ele_iter = mesh.createElementiterator()
		field_fibres = field_fibres.castFiniteElement()
		cache = fieldmodule.createFieldcache()
		node_set = fieldmodule.findNodesetByName("nodes")
		minimum = 1E+5
		ele_index = -1
		x = 0.5
		y = 0.5
		z = 0.5
		x_min = x
		y_min = y
		z_min = z

		local_node_coordinates = [x, y, z]
		#search for element, in wich the input coordinates are
		element = ele_iter.next()
		while(element.isValid()):
			cache.setMeshLocation(element, local_node_coordinates)
			result, ele_coords = field.evaluateReal(cache, 3)
			
			d = distance(ele_coords, ps_coords)
			#print(d)
			if (d < minimum):
				minimum = d
				ele_index = element.getIdentifier()
			element = ele_iter.next()

		element = mesh.findElementByIdentifier(ele_index)
		#split the element in 8 smaller cubes, and search for the cube with the center closest
		#to the input coordinates.
		#continue with the closest cube until the given precision is reached.
		#Example analogous in 2D: x marks the center of each of the four smaller squares.
		#o marks the input coordinates.
		#continue with the new square until the given precision is reached.
		###########################################
		#      |  x  |  x  |		 |  x  |  x  |#
		#STEP1 |_____|_____| -STEP2->|_____|_____|#
		#      |  x  |  x  |		 |  x  |  x  |#
		#      |_____|o____| 		 |o____|_____|#
		###########################################
		refinement = 2
		while((d > precision) and (refinement < 100)):
			for i in [-1,1]:
				for j in [-1,1]:
					for k in [-1,1]:
						x1 = x + i*0.5**(refinement)
						y1 = y + j*0.5**(refinement)
						z1 = z + k*0.5**(refinement) 
						local_node_coordinates = [x1, y1, z1]
						cache.setMeshLocation(element, local_node_coordinates)
						result, ele_coords = field.evaluateReal(cache, 3)
						d = distance(ele_coords, ps_coords)
						if (d < minimum):
							minimum = d
							x_min = x1
							y_min = y1
							z_min = z1
			x = x_min
			y = y_min
			z = z_min
			print(x,y,z)
			refinement += 1
		
		#cache.setFieldReal(field,  ps_coords)
		
		cache.setMeshLocation(element, [x,y,z])
		result, fiber_values = field_fibres.evaluateReal(cache, 3)
		result, ddd = field.evaluateReal(cache, 3)
		#print result.isValid()
		print(distance(ddd, ps_coords))
		return(fiber_values)

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
		#print (field)
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
		#print(str(counter) + " nodes evaluated successfully")
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
		n = ele_template.getNumberOfLocalNodes()
		while element.isValid():
			for i in range(1,n+1):
				
				current_node = element.getNode(ele_template, i)
				cache.setElement(element)
				#test = field.evaluateDerivative(d_ds1, cache, 1)
				#print(test)
				local_index = current_node.getIdentifier()
				out.write(str(local_index) + " ")
			out.write("\r\n")
			element = ele_iter.next()
			ele_template = element.getElementfieldtemplate(field, 1)

	def create_directional_field (self, out_file = "fibres.exfile", refinement1 = 1, refinement2 = 1, refinement3 = 1 ,surface=False):
		fine_mesh_size = 60*2**refinement1 * 2**refinement2 * 2**refinement3 
		print("generating fibres...")
		program_status = ""
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
		#mesh_1d = fieldmodule.findMeshByDimension(1)
		field = field.castFiniteElement()
		field_fibres = field_fibres.castFiniteElement()
		cache = fieldmodule.createFieldcache()
		node_set = fieldmodule.findNodesetByName("nodes")
		node_template = node_set.createNodetemplate()
		node_template.defineField(field)
		node_template.defineField(field_fibres)

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
		node = node_iter.next()
		while node.isValid(): #get number of total nodes befor refinement
			counter += 1
			
			cache.setNode(node)
			result, coords = field.evaluateReal(cache, 3)
			coords = to_euclidian(coords)
			coord_value = coords[0] + coords[1] + coords[2]
			
			bucket_values.append(coord_value)
			node = node_iter.next()
			#node_coordinates.append(coords)
		#print(buckets)
		bucket_max = max(bucket_values)
		bucket_min = min(bucket_values)
		bucket_values = []
		bucket_interval = (bucket_max - bucket_min)/10000.0
		#bucket_values = sorted(bucket_values)
		buckets = {}
		for i in range(0,10000):
			#bucket_value = i*bucket_interval + bucket_min 
			buckets[i] = []
			#bucket_values.append(bucket_value)
		#print (len(bucket_values))
		
		node_iter = node_set.createNodeiterator()
		node = node_iter.next()


		#print(buckets)
		c = 0
		fieldmodule.beginChange()
		node_set.destroyAllNodes()
		#print(node_set.getSize())
		fibers_created = 0



		while c < mesh_size:
			
			element = mesh.findElementByIdentifier(c +1)
			#get elementfieldtemplates to create new elements of the refined mesh
			ele_field_template1 = element.getElementfieldtemplate(field,1)
			
			
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
						"""local_node_coordinates = [[0 + i*l1,0 + j*l2,0 + k*l3], [l1 + i*l1 ,0 + j*l2,0 + k*l3],\
						 [0 + i*l1,l2 + j*l2,0 + k*l3],[l1 + i*l1,l2 + j*l2,0 + k*l3], [0 + i*l1,0 + j*l2,l3 + k*l3],\
						 [l1 + i*l1,0 + j*l2,l3 + k*l3],\
						  [0 + i*l1 ,l2 + j*l2,l3 + k*l3],[l1 + i*l1 ,l2 + j*l2,l3 + k*l3]]"""
						if(surface):
							local_node_coordinates = [[0.5*l1+i*l1, 0.5*l2+j*l2, 1*l3+k*l3]]
						else:	
							local_node_coordinates = [[0.5*l1+i*l1, 0.5*l2+j*l2, 0.5*l3+k*l3]]
                                                #TODO: put this declaration out of loop, maybe do a function flag wether surface or center values are returned...

						#print(local_node_coordinates)
						#create a new node for every local coordinate
						
						for local_coords in local_node_coordinates:
							#get global angle between local coordinate axes
							"""cache.setMeshLocation(element, [0.1,0,0])
							result, xi1 = field.evaluateReal(cache, 3)
							cache.setMeshLocation(element, [0,0.1,0])
							result, xi2 = field.evaluateReal(cache, 3)

							xi1 =  to_euclidian(xi1)
							xi2 = to_euclidian(xi2)
							scalar = 0
							xi1_len = 0
							xi2_len = 0
							for m in range(0,3):
								scalar += xi1[m] * xi2[m]
								xi1_len += xi1[m] * xi1[m]
								xi2_len += xi2[m] * xi2[m]
							xi1_len = xi1_len**0.5
							xi2_len = xi2_len**0.5
							
							gamma = math.acos(scalar/(xi2_len*xi1_len))
							print(gamma)"""
							linear_basis = fieldmodule.createElementbasis(3, 2)
							node_indexes = [1,2]
							#fibre-angle is angle between x1 and x2 axis of local element. add new node to visualize this angle with distance
							#or length of the edge at vec_len
							edge = []
							cache.setMeshLocation(element, local_coords)
							result, global_coords = field.evaluateReal(cache, 3)
							result, global_angles = field_fibres.evaluateReal(cache, 3)
							euclidian_coords = to_euclidian(global_coords)
							coord_value = euclidian_coords[0] + euclidian_coords[1] + euclidian_coords[2]
							#do the same for for the fibre-angles
							result, global_angles = field_fibres.evaluateReal(cache, 3)
							fibre_angle = global_angles[0]
							
							#if (fibre_angle > 3.1415):
							#print(fibre_angle)
							x_angle = local_coords[2]#math.sin(fibre_angle) +  local_coords[2]
							#y_angle =  0.00001*math.sin(  gamma -1.57079633 + fibre_angle) + local_coords[1]
							
							slope = math.tan(fibre_angle)
							
							d = 0.0001*math.cos(fibre_angle)
							##y_angle =  0.00001*math.sin(fibre_angle) + local_coords[1]
							##z_angle =  0.00001*math.cos(fibre_angle) + local_coords[0]
							y_angle = local_coords[1]  + slope*(d)
							z_angle = local_coords[0] + d
							
							
							bucket = 0
							if(y_angle>10):
								print(slope)
							bucket_value = (coord_value - bucket_min)//bucket_interval
							if bucket_value < 0:
								bucket_value = 0
							if bucket_value > 9999:
								bucket_value = 9999
							#print(bucket)
							#check if coordinates belong to an already existing node by comparing euclidian distance
							already_existing = False
							#print(bucket_value)
							unified_buckets = []
							if bucket_value == 0:
								unified_buckets = buckets[bucket_value] + buckets[bucket_value+1]
							elif bucket_value == 9999:
								unified_buckets = buckets[bucket_value-1] + buckets[bucket_value]
							else:
								unified_buckets = buckets[bucket_value-1] + buckets[bucket_value] + buckets[bucket_value +1]
							for node_test in unified_buckets:
								
								#print(distance(global_coords,coordinates))
								cache.clearLocation()
								cache.setNode(node_test)
								coord_test = []
								result, coord_test = field.evaluateReal(cache, 3)
								if (distance(global_coords,coord_test) < 1E-4):
									already_existing = True
									
									break
							


							#create new node with the interpolated global coordinates and interpolated fibre-angles
							if (already_existing == False):
								#create interpolated node
								fibers_created+= 1
								node = node_set.createNode(counter + 1, node_template)

								cache.setNode(node)
								#print(global_coords)
								field.assignReal(cache, global_coords)
								field_fibres.assignReal(cache, global_angles)
								edge.append(counter + 1)
								buckets[bucket_value].append(node)
								counter += 1
								#create node to visualize fibre angle
								coords_angle = [z_angle, y_angle, x_angle]

								#print(coords_angle)
								cache.setMeshLocation(element, coords_angle)
								result, global_coords = field.evaluateReal(cache, 3)
								#print(global_coords)
								node_edge = node_set.createNode(counter + 1, node_template)
								cache.setNode(node_edge)
								field.assignReal(cache, global_coords)
								edge.append(counter + 1)
								counter += 1

									
			#print(str(num_ele) + " of " + str(60 * 2**refinement1 * 2**refinement2 * 2**refinement3 -60) + " elements created.")
			mesh.destroyElement(element)
			progress = (1.0)*(fibers_created)/ (fine_mesh_size)
			update_progress(progress, program_status)
			c += 1
		for i in range(1,100):
			node_dest = node_set.findNodeByIdentifier(i)
			node_set.destroyNode(node_dest)
		region.writeFile(out_file)
		#print("done")
		fieldmodule.endChange()		
		
		


		num_nodes = node_set.getSize()	#get number of vertices in nodeset


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

		mesh = fieldmodule.findMeshByDimension(3)
		ele_iter = mesh.createElementiterator()
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
		#print(str(counter) + " nodes evaluated successfully")
	def kd_interpolation (self, out_file = "interpolation.exfile", refinement1 = 1, refinement2 = 1, refinement3 = 1):
		fine_mesh_size = 60*2**refinement1 * 2**refinement2 * 2**refinement3 
		print("generating mesh...")
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
		node = node_iter.next()
		while node.isValid(): #get number of total nodes befor refinement
			counter += 1
			
			cache.setNode(node)
			result, coords = field.evaluateReal(cache, 3)
			coords = to_euclidian(coords)
			coord_value = coords[0] + coords[1] + coords[2]
			
			bucket_values.append(coord_value)
			node = node_iter.next()
			#node_coordinates.append(coords)
		#print(buckets)
		bucket_max = max(bucket_values)
		bucket_min = min(bucket_values)
		bucket_values = []
		bucket_interval = (bucket_max - bucket_min)/10000.0
		#bucket_values = sorted(bucket_values)
		buckets = {}
		for i in range(0,10000):
			#bucket_value = i*bucket_interval + bucket_min 
			buckets[i] = []
			#bucket_values.append(bucket_value)
		#print (len(bucket_values))
		#print(bucket_values[0])
		#print(bucket_min)
		node_iter = node_set.createNodeiterator()
		node = node_iter.next()
		while node.isValid(): #get number of total nodes befor refinement
			cache.setNode(node)
			result, coords = field.evaluateReal(cache, 3)
			coords = to_euclidian(coords)
			coord_value = coords[0] + coords[1] + coords[2]
			bucket_value = (coord_value - bucket_min)//bucket_interval
			#bucket_value = bucket_value * bucket_interval
			buckets[bucket_value].append(node)
			node = node_iter.next()
			#print(bucket_value)
			#print(coord_value)


		#print(buckets)
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
							#print(coord_value)
							bucket_value = (coord_value - bucket_min)//bucket_interval
							if bucket_value < 0:
								bucket_value = 0
							if bucket_value > 9999:
								bucket_value = 9999
							#print(bucket)
							#check if coordinates belong to an already existing node by comparing euclidian distance
							already_existing = False
							#print(bucket_value)
							unified_buckets = []
							if bucket_value == 0:
								unified_buckets = buckets[bucket_value] + buckets[bucket_value+1]
							elif bucket_value == 9999:
								unified_buckets = buckets[bucket_value-1] + buckets[bucket_value]
							else:
								unified_buckets = buckets[bucket_value-1] + buckets[bucket_value] + buckets[bucket_value +1]
							for node_test in unified_buckets:
								
								#print(distance(global_coords,coordinates))
								cache.clearLocation()
								cache.setNode(node_test)
								coord_test = []
								result, coord_test = field.evaluateReal(cache, 3)
								if (distance(global_coords,coord_test) < 1E-4):
									already_existing = True
									#print("TRIGGERED")
									node_identifiers.append(node_test.getIdentifier())
									#print(node_test.isValid())
									break
							


							#create new node with the interpolated global coordinates and interpolated fibre-angles
							if (already_existing == False):
								node = node_set.createNode(counter + 1, node_template)
								#print(node.isValid())
								#print(node.isValid())
								#cache.clearLocation()
								cache.setNode(node)
								field.assignReal(cache, global_coords)
								field_fibres.assignReal(cache, global_angles)
								node_identifiers.append(counter + 1)
								buckets[bucket_value].append(node)
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
			program_status = str(num_ele - 60) + " of " + str(fine_mesh_size) + " elements created."
			progress = (1.0)*(num_ele - 60)/ (fine_mesh_size)
			update_progress(progress, program_status)
			mesh.destroyElement(element)
			c += 1

		region.writeFile(out_file)
		#print("done")
		fieldmodule.endChange()		
		
#ex_reader = exfile_reader("heart1.exfile")
#values = ex_reader.get_fiber_values([-22.762729289,11.152001135,-39.54201334], 1E-5)
#print (values)
