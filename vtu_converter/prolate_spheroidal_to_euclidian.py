import os
import math
def prolate_spheroid_to_euclidian(inFile, outFile, focus):
	inCoords = open(inFile)
	if os.path.isfile("/tmp/" + outFile):
		os.remove("/tmp/" + outFile)
	outCoords = open(outFile , "w+")
	line = inCoords.readline()
	while (line):
		coords = line.split(" ")
		lamda = float(coords[0])
		mu = float(coords[1])
		theta = float(coords[2])
		x = focus * math.cosh(lamda) * math.cos(mu)
		y = focus * math.sinh(lamda) * math.sin(mu) * math.cos(theta)
		z = focus * math.sinh(lamda) * math.sin(mu) * math.sin(theta)
		outCoords.write(str(x) + " " + str(y) + " " + str(z) + " " + "\r\n")
		line = inCoords.readline()
#prolate_spheroid_to_euclidian("coordinates_prolate_spheroidal.txt", "nodes_euclidian.txt", 0.3525E+02)