# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================



#================================    
#Gamma3 = A4 Clebsch-Gordan table
#================================ 

def SymbolToIntegerGamma3(rep):
    if len(rep) < 2:
        return 0
    elif rep[-1] == "'":
        return 1
    else:
        return 2

def IntegerToSymbolGamma3(integer):
    if integer == 0:
        return ""
    elif integer == 1:
        return "'"
    else:
        return '"'
    
def ProductGamma3(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim1 == 1 and dim2 == 1:
        Symbol1 = SymbolToIntegerGamma3(rep1)
        Symbol2 = SymbolToIntegerGamma3(rep2)
        Symbol = (Symbol1+Symbol2)%3
        Symbol = IntegerToSymbolGamma3(Symbol)
        Component = [Components1[0]*Components2[0]]
        rep = "1"+Symbol
        return [[rep, Component]]
    elif dim1 == 1 and dim2 == 3:
        Symbol = SymbolToIntegerGamma3(rep1)
        rep = "3"
        Components = []
        for i in range(3):
            Components.append(Components1[0]*Components2[(i-Symbol)%3])
        return [[rep,Components]]
    elif dim1 == 3 and dim2 == 1:
        Symbol = SymbolToIntegerGamma3(rep2)
        rep = "3"
        Components = []
        for i in range(3):
            Components.append(Components2[0]*Components1[(i-Symbol)%3])
        return [[rep,Components]]
    elif dim1 == 3 and dim2 == 3:
        Decomposition = []
        Decomposition.append(["1",[Components1[0]*Components2[0]+Components1[1]*Components2[2]+Components1[2]*Components2[1]]])
        Decomposition.append(["1"+"'",[Components1[2]*Components2[2]+Components1[0]*Components2[1]+Components1[1]*Components2[0]]])
        Decomposition.append(["1"+'"',[Components1[1]*Components2[1]+Components1[0]*Components2[2]+Components1[2]*Components2[0]]])
        Components3S = [2*Components1[0]*Components2[0]-Components1[1]*Components2[2]-Components1[2]*Components2[1]]
        Components3S.append(2*Components1[2]*Components2[2]-Components1[0]*Components2[1]-Components1[1]*Components2[0])
        Components3S.append(2*Components1[1]*Components2[1]-Components1[0]*Components2[2]-Components1[2]*Components2[0])
        Decomposition.append(["3",Components3S])
        Components3A = [Components1[2]*Components2[1]-Components1[1]*Components2[2]]
        Components3A.append(Components1[1]*Components2[0]-Components1[0]*Components2[1])
        Components3A.append(Components1[0]*Components2[2]-Components1[2]*Components2[0])
        Decomposition.append(["3",Components3A])
        return Decomposition