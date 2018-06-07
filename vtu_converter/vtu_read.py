

def vtu_read(vtu_file):
    file = open(vtu_file)
    
    line = file.readline()
    
    list = line.split(" ")
    i = 0
    while (list[0] == '' and i < len(list)): #ignore inentation
        i += 1
        list.pop(0)
    print(list)
    line = file.readline()
    
    list = line.split(" ")
    i = 0
    while (list[0] == '' and i < len(list)): #ignore inentation
        list.pop(0)
        i += 1
    print(list)
    line = file.readline()
    
    list = line.split(" ")
    i = 0
    while (list[0] == '' and i < len(list)): #ignore inentation
        list.pop(0)
        i += 1
    print(list)
    
vtu_read("newObject.vtu")