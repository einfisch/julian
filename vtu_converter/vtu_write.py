import os
import exfile_reader
import math
def fineMeshWedges (con, coo):
    #creates Data for a mesh with 8 vertices per 27-element cell
    # uses wedges for the bottom cells
    if os.path.isfile("/tmp/region.txt"):
        os.remove("/tmp/region.txt")
    if os.path.isfile("/tmp/cellTypes.txt"):
        os.remove("/tmp/cellTypesVTU.txt")
    if os.path.isfile("/tmp/offset.txt"):
        os.remove("/tmp/coffset.txt")
    if os.path.isfile("/tmp/cells.txt"):
        os.remove("/tmp/cells.txt")
    cellTypeFile = open("cellTypes.txt", "w+")
    regionFile = open("region.txt", "w+")
    offsetFile = open("offset.txt", "w+")
    cellFile = open("cells.txt", "w+")
    offset = 0
    region = 0
    numPoints = 0
    numCells = 0
    line = coo.readline()
    #regionFile.write("0 ")#test
    #cellTypeFile.write("1 ")#testi
    #count number of vertices to get offset
    while (line):
        
        offset += 1
        offsetFile.write(str(offset) + " ")
        cellFile.write(str(numPoints) + " " )
        numPoints += 1
        cellTypeFile.write("1 ")
        regionFile.write("0 ")
        line = coo.readline()
        
    #iterate for every line in the connection-file:
    line = con.readline()
    
    while (line):
        region += 1 #new subset
        reS = str(region)+ " " #number of subset of the active cell
        s = line.split(" ")
        if s[0] == s[1]: #if first and second node are the same, we have wedge cell
            numCells += 8
            cell0 = str(int(s[0])-1)+" " +str(int(s[3])-1)+" " +str(int(s[4])-1) + " "+str(int(s[9])-1)+" "+str(int(s[12])-1)+" "+str(int(s[13])-1)
            cell1 = str(int(s[1])-1)+" " +str(int(s[4])-1)+" " +str(int(s[5])-1) + " "+str(int(s[9])-1)+" "+str(int(s[13])-1)+" "+str(int(s[14])-1)
     
            cell2 = str(int(s[9])-1)+" " +str(int(s[12])-1)+" " +str(int(s[13])-1) + " "+str(int(s[18])-1)+" "+str(int(s[21])-1)+" "+str(int(s[22])-1)
            cell3 = str(int(s[9])-1)+" " +str(int(s[13])-1)+" " +str(int(s[14])-1) + " "+str(int(s[19])-1)+" "+str(int(s[22])-1)+" "+str(int(s[23])-1)
    
            cell4 = str(int(s[3])-1)+" " +str(int(s[12])-1)+" " +str(int(s[13])-1) + " "+str(int(s[4])-1)+" "+str(int(s[6])-1)+" "+str(int(s[15])-1)+" "+str(int(s[16])-1)+" "+str(int(s[7])-1)
            cell5 = str(int(s[4])-1)+" " +str(int(s[13])-1)+" " +str(int(s[14])-1) + " "+str(int(s[5])-1)+" "+str(int(s[7])-1)+" "+str(int(s[16])-1)+" "+str(int(s[17])-1)+" "+str(int(s[8])-1)
            cell6 = str(int(s[12])-1)+" " +str(int(s[21])-1)+" " +str(int(s[22])-1) + " "+str(int(s[13])-1)+" "+str(int(s[15])-1)+" "+str(int(s[24])-1)+" "+str(int(s[25])-1)+" "+str(int(s[16])-1)
            cell7 = str(int(s[13])-1)+" " +str(int(s[22])-1)+" " +str(int(s[23])-1) + " "+str(int(s[14])-1)+" "+str(int(s[16])-1)+" "+str(int(s[25])-1)+" "+str(int(s[26])-1)+" "+str(int(s[17])-1)
            cellTypeFile.write("13 13 13 13 12 12 12 12 \r\n")
            offset += 6
            offsetFile.write(str(offset) + " ")
            offset += 6
            offsetFile.write(str(offset) + " ")
            offset += 6
            offsetFile.write(str(offset) + " ")
            offset += 6
            offsetFile.write(str(offset) + " ")
            offset += 8
            offsetFile.write(str(offset) + " ")
            offset += 8
            offsetFile.write(str(offset) + " ")
            offset += 8
            offsetFile.write(str(offset) + " ")
            offset += 8
            offsetFile.write(str(offset) + " \r\n")
            
            regionFile.write(reS + reS + reS + reS + reS + reS + reS +reS +"\r\n")
            
            
        else:
            numCells += 8
            cell0 = str(int(s[0])-1)+" " +str(int(s[9])-1)+" " +str(int(s[10])-1) + " "+str(int(s[1])-1)+" "+str(int(s[3])-1)+" "+str(int(s[12])-1)+" "+str(int(s[13])-1)+" "+str(int(s[4])-1)
            cell1 = str(int(s[1])-1)+" " +str(int(s[10])-1)+" " +str(int(s[11])-1) + " "+str(int(s[2])-1)+" "+str(int(s[4])-1)+" "+str(int(s[13])-1)+" "+str(int(s[14])-1)+" "+str(int(s[5])-1)
            cell2 = str(int(s[3])-1)+" " +str(int(s[12])-1)+" " +str(int(s[13])-1) + " "+str(int(s[4])-1)+" "+str(int(s[6])-1)+" "+str(int(s[15])-1)+" "+str(int(s[16])-1)+" "+str(int(s[7])-1)
            cell3 = str(int(s[4])-1)+" " +str(int(s[13])-1)+" " +str(int(s[14])-1) + " "+str(int(s[5])-1)+" "+str(int(s[7])-1)+" "+str(int(s[16])-1)+" "+str(int(s[17])-1)+" "+str(int(s[8])-1)
    
            cell4 = str(int(s[9])-1)+" " +str(int(s[18])-1)+" " +str(int(s[19])-1) + " "+str(int(s[10])-1)+" "+str(int(s[12])-1)+" "+str(int(s[21])-1)+" "+str(int(s[22])-1)+" "+str(int(s[13])-1)
            
            cell5 = str(int(s[10])-1)+" " +str(int(s[19])-1)+" " +str(int(s[20])-1) + " "+str(int(s[11])-1)+" "+str(int(s[13])-1)+" "+str(int(s[22])-1)+" "+str(int(s[23])-1)+" "+str(int(s[14])-1)
            cell6 = str(int(s[12])-1)+" " +str(int(s[21])-1)+" " +str(int(s[22])-1) + " "+str(int(s[13])-1)+" "+str(int(s[15])-1)+" "+str(int(s[24])-1)+" "+str(int(s[25])-1)+" "+str(int(s[16])-1)
            cell7 = str(int(s[13])-1)+" " +str(int(s[22])-1)+" " +str(int(s[23])-1) + " "+str(int(s[14])-1)+" "+str(int(s[16])-1)+" "+str(int(s[25])-1)+" "+str(int(s[26])-1)+" "+str(int(s[17])-1)
            cellTypeFile.write("12 12 12 12 12 12 12 12 \r\n")
            
            offset += 8
            offsetFile.write(str(offset) + " ")
            offset += 8
            offsetFile.write(str(offset) + " ")
            offset += 8
            offsetFile.write(str(offset) + " ")
            offset += 8
            offsetFile.write(str(offset) + " ")
            offset += 8
            offsetFile.write(str(offset) + " ")
            offset += 8
            offsetFile.write(str(offset) + " ")
            offset += 8
            offsetFile.write(str(offset) + " \r\n")
            offset += 8
            offsetFile.write(str(offset) + " \r\n")
            regionFile.write(reS + reS + reS + reS + reS + reS + reS + reS + "\r\n")
            
        cell = cell0 + "\r\n" + cell1 + "\r\n" + cell2 + "\r\n" + cell3 + "\r\n"+ cell4 + "\r\n" + cell5 + "\r\n" + cell6 + "\r\n" + cell7 + "\r\n"    
        cellFile.write(cell + '\n')
        
        line = con.readline()
    return [numPoints, numCells, region]



def wideMeshWedges (con, coo):
    #creates Data for a mesh with 8 vertices per 27-element cell
    # uses wedges for the bottom cells
    current_directory = os.getcwd()
    directory = os.path.join(current_directory, r"vtu_data")
    if not os.path.exists(directory):
        os.makedirs(directory)
    if os.path.isfile(os.path.join(directory, "region.txt")):
        os.remove(os.path.join(directory, "region.txt"))
    if os.path.isfile(os.path.join(directory, "cellTypesVTU.txt")):
        os.remove(os.path.join(directory, "cellTypesVTU.txt"))
    if os.path.isfile(os.path.join(directory, "offset.txt")):
        os.remove(os.path.join(directory, "offset.txt"))
    if os.path.isfile(os.path.join(directory, "cells.txt")):
        os.remove(os.path.join(directory, "cells.txt"))
    cellTypeFile = open(os.path.join(directory, "cellTypesVTU.txt"), "w+")
    regionFile = open(os.path.join(directory, "region.txt"), "w+")
    offsetFile = open(os.path.join(directory, "offset.txt"), "w+")
    cellFile = open(os.path.join(directory, "cells.txt"), "w+")
    offset = 0
    region = 0
    numPoints = 0
    numCells = 0
    line = coo.readline()
    #count number of vertices to get offset
    while (line):
        
        offset += 1
        offsetFile.write(str(offset) + " ")
        cellFile.write(str(numPoints) + " " )
        numPoints += 1
        cellTypeFile.write("1 ")
        regionFile.write("0 ")
        line = coo.readline()
        
    #iterate for every line in the connection-file:
    line = con.readline()
    while (line):
        
        reS = str(region) + " "
        s = line.split(" ")
        if s[7] == "-1": #if last two nodes are -1 we have a wedge-shaped cell
            numCells += 1
            cell = str(int(s[0])-1)+" " +str(int(s[1])-1)+" " +str(int(s[2])-1) + " "+str(int(s[3])-1)+" "+str(int(s[4])-1)+" "+str(int(s[5])-1) + " "
            #cell = str(int(s[0])-1)+" " +str(int(s[6])-1)+" " +str(int(s[8])-1) + " "+str(int(s[18])-1)+" "+str(int(s[24])-1)+" "+str(int(s[26])-1) + " "
            cellTypeFile.write("13 ")
            offset += 6
            offsetFile.write(str(offset) + " ")
            
            
        else:
            numCells += 1
            cell = str(int(s[0])-1)+" " +str(int(s[4])-1)+" " +str(int(s[5])-1) + " "+str(int(s[1])-1)+" "+str(int(s[2])-1)+" "+str(int(s[6])-1)+" "+str(int(s[7])-1)+" "+str(int(s[3])-1) + " "
            cellTypeFile.write("12 ")
            offset += 8
            offsetFile.write(str(offset) + " " )
            
            
        cellFile.write(cell + '\n')
        regionFile.write(reS + "\n")
        line = con.readline()
    return [numPoints, numCells, region]

def directional_field(coo, r):
    #function takes euclidian coordinates and fibre angles and returns data for vtu_coverter for the directional field
    current_directory = os.getcwd()
    directory = os.path.join(current_directory, r"vtu_data")
    if not os.path.exists(directory):
        os.makedirs(directory)
    if os.path.isfile(os.path.join(directory, "region_angles.txt")):
        os.remove(os.path.join(directory, "region_angles.txt"))
    if os.path.isfile(os.path.join(directory, "vertices_angles.txt")):
        os.remove(os.path.join(directory, "vertices_angles.txt"))
    if os.path.isfile(os.path.join(directory, "cellTypesVTU_angles.txt")):
        os.remove(os.path.join(directory, "cellTypesVTU_angles.txt"))
    if os.path.isfile(os.path.join(directory, "offset_angles.txt")):
        os.remove(os.path.join(directory, "offset_angles.txt"))
    if os.path.isfile(os.path.join(directory, "cells_angles.txt")):
        os.remove(os.path.join(directory, "cells_angles.txt"))
    cellTypeFile = open(os.path.join(directory, "cellTypesVTU_angles.txt"), "w+")
    regionFile = open(os.path.join(directory, "region_angles.txt"), "w+")
    offsetFile = open(os.path.join(directory, "offset_angles.txt"), "w+")
    cellFile = open(os.path.join(directory, "cells_angles.txt"), "w+")
    vert = open(os.path.join(directory, "vertices_angles.txt"), "w+")
    offset = 0
    region = 0
    numPoints = 0
    numCells = 0
    line_coords = coo.readline()
    
    connectivity = []
    toggle = False
    while (line_coords):
        coords = line_coords.split(" ")

        cellTypeFile.write("1 ")
        offset += 1
        offsetFile.write(str(offset ) + " ")
        
        regionFile.write(str(0)  + " ")
        
        cellFile.write(str(numPoints  ) + " " )
        
        
        if (toggle == True):
            connectivity.append([numPoints - 1 , numPoints -0 ])
            #use linear interpolation to create vertex with distance equal to r
            
            x = float(old_coords[0])
            y = float(old_coords[1])
            z = float(old_coords[2])

            #vert.write(str(x) + " " + str(y) + " " + str(z) +"\r\n")

            x1 = float(coords[0])
            y1 = float(coords[1])
            z1 = float(coords[2])
            v1 = x1 - x
            v2 = y1 - y
            v3 = z1 - z

            
            ###
            b = v1*v1 + v2*v2 + v3*v3
            ###
            lam = 0
            if (not(b == 0)):
                lam = r/(b**0.5)
            #print(lam) 
            x2 = x - 0.5*lam*v1
            y2 = y - 0.5*lam*v2
            z2 = z - 0.5*lam*v3
            x3 = x+ 0.5*lam*v1
            y3 = y+0.5*lam*v2
            z3 = z+ 0.5*lam*v3

            #print(((x-x3)**2+(y-y3)**2+(z-z3)**2)**0.5)
            vert.write(str(x2) + " " + str(y2) + " " + str(z2) +"\r\n")
            vert.write(str(x3) + " " + str(y3) + " " + str(z3) +"\r\n")
        toggle = not toggle
        numPoints += 1
        line_coords = coo.readline()
        old_coords = coords
    for edge in connectivity:
        numCells += 1
        cellTypeFile.write("3 ")
        offset += 2
        offsetFile.write(str(offset - 0) + " " )
        cellFile.write(str(edge[0]) +" " +str(edge[1]) + "\r\n")
        regionFile.write("0 ")
    return [numPoints, numCells, 0]



    
def create_file(coo, numPoints, numCells, subsets, out, angles):    
    
    print("exporting " + out +" ...")
    current_directory = os.getcwd()
    directory = os.path.join(current_directory, r"vtu_data")
    vtuOut = open(out, "w+")
   

   

    #hardcode syntax of VTU-file format

    vtuOut.write("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"> \n")
    vtuOut.write("<UnstructuredGrid> \n")

    vtuOut.write("<Piece NumberOfPoints=\""+str(numPoints)+"\" NumberOfCells=\""+str(numCells)+"\">\n")
    vtuOut.write("<Points>\n")
    vtuOut.write("<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n")
    #write point coordinates to VTU-File
    if angles:
        coordFile = open(os.path.join(directory, coo))
    else:
        coordFile = open(coo, "r")
    line = coordFile.readline()

    while line:
        vtuOut.write(line)
        line = coordFile.readline()
        

    vtuOut.write("</DataArray>\n")
    vtuOut.write("</Points>\n")
    vtuOut.write("<Cells>\n")
    vtuOut.write(" <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
    if angles:
        connectivity = open(os.path.join(directory, "cells_angles.txt"))
    else:
        connectivity = open(os.path.join(directory, "cells.txt"))
    line = connectivity.readline()
    while line:
        vtuOut.write(line)
        line = connectivity.readline()

    vtuOut.write("</DataArray>\n")
    vtuOut.write("<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
    if angles:
        offset = open(os.path.join(directory, "offset_angles.txt"))
    else:
        offset = open(os.path.join(directory, "offset.txt"))
    line = offset.readline()
    while line:
        vtuOut.write(line)
        line = offset.readline()
    vtuOut.write("</DataArray>\n")
    vtuOut.write("<DataArray type=\"Int8\" Name=\"types\" format=\"ascii\">\n")
    
    if angles:
        cellTypes = open(os.path.join(directory, "cellTypesVTU_angles.txt"))
    else:    
        cellTypes = open(os.path.join(directory, "cellTypesVTU.txt"))
    line = cellTypes.readline()
    while line:
        vtuOut.write(line)
        line = cellTypes.readline()
    vtuOut.write("</DataArray>\n")
    vtuOut.write("</Cells>\n")
    vtuOut.write("<PointData>\n")
    vtuOut.write("</PointData>\n")
    vtuOut.write("<CellData>\n")
    vtuOut.write("<DataArray type=\"Int32\" Name=\"regions\" NumberOfComponents=\"1\" format=\"ascii\">\n")
    if angles:
        region = open(os.path.join(directory, "region_angles.txt"))
    else:
        region = open(os.path.join(directory, "region.txt"))
    line = region.readline()
    while line:
        vtuOut.write(line)
        line = region.readline()
    vtuOut.write("</DataArray>\n")
    vtuOut.write("</CellData>\n")
    vtuOut.write("<RegionInfo Name=\"regions\">\n")
    vtuOut.write("<Region Name=\"subset\"></Region>\n")#possibly not needed
    i = 0
    while i < subsets:
        i += 1
        vtuOut.write("<Region Name=\"cell " +str(i) + "\"></Region>\n")
       
    vtuOut.write("</RegionInfo>\n")
    vtuOut.write("</Piece>")
    vtuOut.write("</UnstructuredGrid> \n")
    vtuOut.write("</VTKFile>")

