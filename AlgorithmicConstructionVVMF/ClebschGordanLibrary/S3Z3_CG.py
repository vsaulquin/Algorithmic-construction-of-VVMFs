# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================



#==========================    
#S3xZ3 Clebsch-Gordan table
#========================== 

def SymbolToIntegerS3Z3(rep):
    if len(rep) < 4:
        return 0
    else:
        return 1

def IntegerToSymbolS3Z3(integer):
    if integer == 0:
        return ""
    elif integer == 1:
        return "'"

def ProductS3Z3(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    copy1 = int(rep1[2])
    copy2 = int(rep2[2])
    copy = str((copy1+copy2)%3)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim1 == 1 and dim2 == 1:
        Symbol1 = SymbolToIntegerS3Z3(rep1)
        Symbol2 = SymbolToIntegerS3Z3(rep2)
        Symbol = (Symbol1+Symbol2)%2
        Symbol = IntegerToSymbolS3Z3(Symbol)
        Component = [Components1[0]*Components2[0]]
        rep = "1_"+copy+Symbol
        return [[rep, Component]]
    elif dim1 == 1 and dim2 == 2:
        Symbol = SymbolToIntegerS3Z3(rep1)
        Components = []
        if Symbol == 0:
            return [["2_"+copy,[Components1[0]*Components2[0],Components1[0]*Components2[1]]]]
        else:
            return [["2_"+copy,[-Components1[0]*Components2[1],Components1[0]*Components2[0]]]]
    elif dim1 == 2 and dim2 == 1:
        Symbol = SymbolToIntegerS3Z3(rep1)
        Components = []
        if Symbol == 0:
            return [["2_"+copy,[Components2[0]*Components1[0],Components2[0]*Components1[1]]]]
        else:
            return [["2_"+copy,[-Components2[0]*Components1[1],Components2[0]*Components1[0]]]]
    elif dim1 == 2 and dim2 == 2:
        Decomposition = []
        Decomposition.append(["1_"+copy,[Components1[0]*Components2[0]+Components1[1]*Components2[1]]])
        Decomposition.append(["1_"+copy+"'",[Components1[0]*Components2[1]-Components1[1]*Components2[0]]])
        Decomposition.append(["2_"+copy,[Components1[0]*Components2[0]-Components1[1]*Components2[1],-Components1[0]*Components2[1]-Components1[1]*Components2[0]]])
        return Decomposition
