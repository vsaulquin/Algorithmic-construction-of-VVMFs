# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================



#==========================    
#A4xZ2 Clebsch-Gordan table
#========================== 

def SymbolToIntegerA4Z2(rep):
    if len(rep) < 4:
        return 0
    elif rep[-1] == "'":
        return 1
    else:
        return 2

def IntegerToSymbolA4Z2(integer):
    if integer == 0:
        return ""
    elif integer == 1:
        return "'"
    else:
        return '"'
    
def ProductA4Z2(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim1 == 1 and dim2 == 1:
        copy1 = int(rep1[2])
        copy2 = int(rep2[2])
        copy = str((copy1+copy2)%2)
        Symbol1 = SymbolToIntegerA4Z2(rep1)
        Symbol2 = SymbolToIntegerA4Z2(rep2)
        Symbol = (Symbol1+Symbol2)%3
        Symbol = IntegerToSymbolA4Z2(Symbol)
        Component = [Components1[0]*Components2[0]]
        rep = "1_"+copy+Symbol
        return [[rep, Component]]
    elif dim1 == 1 and dim2 == 3:
        copy1 = int(rep1[2])
        copy2 = int(rep2[2])
        copy = str((copy1+copy2)%2)
        Symbol = SymbolToIntegerA4Z2(rep1)
        rep = "3_"+copy
        Components = []
        for i in range(3):
            Components.append(Components1[0]*Components2[(i-Symbol)%3])
        return [[rep,Components]]
    elif dim1 == 3 and dim2 == 1:
        copy1 = int(rep1[2])
        copy2 = int(rep2[2])
        copy = str((copy1+copy2)%2)
        Symbol = SymbolToIntegerA4Z2(rep2)
        rep = "3_"+copy
        Components = []
        for i in range(3):
            Components.append(Components2[0]*Components1[(i-Symbol)%3])
        return [[rep,Components]]
    elif dim1 == 3 and dim2 == 3:
        copy1 = int(rep1[2])
        copy2 = int(rep2[2])
        copy = str((copy1+copy2)%2)
        Decomposition = []
        Decomposition.append(["1_"+copy,[Components1[0]*Components2[0]+Components1[1]*Components2[2]+Components1[2]*Components2[1]]])
        Decomposition.append(["1_"+copy+"'",[Components1[2]*Components2[2]+Components1[0]*Components2[1]+Components1[1]*Components2[0]]])
        Decomposition.append(["1_"+copy+'"',[Components1[1]*Components2[1]+Components1[0]*Components2[2]+Components1[2]*Components2[0]]])
        Components3S = [2*Components1[0]*Components2[0]-Components1[1]*Components2[2]-Components1[2]*Components2[1]]
        Components3S.append(2*Components1[2]*Components2[2]-Components1[0]*Components2[1]-Components1[1]*Components2[0])
        Components3S.append(2*Components1[1]*Components2[1]-Components1[0]*Components2[2]-Components1[2]*Components2[0])
        Decomposition.append(["3_"+copy,Components3S])
        Components3A = [Components1[2]*Components2[1]-Components1[1]*Components2[2]]
        Components3A.append(Components1[1]*Components2[0]-Components1[0]*Components2[1])
        Components3A.append(Components1[0]*Components2[2]-Components1[2]*Components2[0])
        Decomposition.append(["3_"+copy,Components3A])
        return Decomposition