import os
import math
import numpy as np
import time
import sys
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

def to_prolate_spheroidal(coords):
	focus = 0.3525E+02
	x = coords[0]
	y = coords[1]
	z = coords[2]

	theta = math.atan(z/y)
	mu = math.acos(1/(2*focus)*((y*y + z*z + (x + focus)**2)**0.5 - (y*y + z*z + (x - focus)**2)**0.5))
	lam = math.acosh(1/(2*focus)*((y*y + z*z + (x + focus)**2)**0.5 + (y*y + z*z + (x - focus)**2)**0.5))

	return([lam, mu, theta])
def get_fiber_values( x,y,z, precision):
		#returns the fiber vector at the specified euclidian coordinates (list like [x,y,z])
	coords = [x,y,z]
	ps_coords = to_prolate_spheroidal(coords)#convert to prolate spheroidal coordinates
	print(coords)
	print(ps_coords)
		
	context = Context("heart")
	region = context.getDefaultRegion()
	region.readFile("heart1.exfile")
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
	print(fiber_values)
	return(fiber_values)


#ex_reader = exfile_reader("heart1.exfile")
#values = ex_reader.get_fiber_values([-22.762729289,11.152001135,-39.54201334], 1E-5)
#print (values)
