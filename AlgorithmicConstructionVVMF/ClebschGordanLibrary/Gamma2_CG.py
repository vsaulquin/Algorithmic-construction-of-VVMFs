# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================



#================================    
#Gamma2 = S3 Clebsch-Gordan table
#================================
 
def SymbolToIntegerGamma2(rep):
    if len(rep) < 2:
        return 0
    else:
        return 1

def IntegerToSymbolGamma2(integer):
    if integer == 0:
        return ""
    elif integer == 1:
        return "'"

def ProductGamma2(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim1 == 1 and dim2 == 1:
        Symbol1 = SymbolToIntegerGamma2(rep1)
        Symbol2 = SymbolToIntegerGamma2(rep2)
        Symbol = (Symbol1+Symbol2)%2
        Symbol = IntegerToSymbolGamma2(Symbol)
        Component = [Components1[0]*Components2[0]]
        rep = "1"+Symbol
        return [[rep, Component]]
    elif dim1 == 1 and dim2 == 2:
        Symbol = SymbolToIntegerGamma2(rep1)
        Components = []
        if Symbol == 0:
            return [["2",[Components1[0]*Components2[0],Components1[0]*Components2[1]]]]
        else:
            return [["2",[-Components1[0]*Components2[1],Components1[0]*Components2[0]]]]
    elif dim1 == 2 and dim2 == 1:
        Symbol = SymbolToIntegerGamma2(rep1)
        if Symbol == 0:
            return [["2",[Components2[0]*Components1[0],Components2[0]*Components1[1]]]]
        else:
            return [["2",[-Components2[0]*Components1[1],Components2[0]*Components1[0]]]]
    elif dim1 == 2 and dim2 == 2:
        Decomposition = []
        Decomposition.append(["1",[Components1[0]*Components2[0]+Components1[1]*Components2[1]]])
        Decomposition.append(["1"+"'",[Components1[0]*Components2[1]-Components1[1]*Components2[0]]])
        Decomposition.append(["2",[Components1[0]*Components2[0]-Components1[1]*Components2[1],-Components1[0]*Components2[1]-Components1[1]*Components2[0]]])
        return Decomposition
