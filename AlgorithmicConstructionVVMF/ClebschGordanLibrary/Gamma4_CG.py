# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================



#================================
#Gamma4 = S4 Clebsch-Gordan table
#================================

def SymbolToIntegerGamma4(rep):
    if rep[0] == "1" or rep[0]=="3":
        if len(rep) == 1:
            return 0
        else:
            return 1
    else:
        return 0

def IntegerToSymbolGamma4(Integer):
    if Integer == 0:
        return ""
    else:
        return "'"

def ProductGamma411(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerGamma4(rep1)
    Symbol2 = SymbolToIntegerGamma4(rep2)
    Symbol = IntegerToSymbolGamma4((Symbol1+Symbol2)%2)
    return [["1"+Symbol,[Components1[0]*Components2[0]]]]

def ProductGamma412(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerGamma4(rep1)
    if Symbol1 == 0:
        return [["2",[Components1[0]*Components2[0],Components1[0]*Components2[1]]]]
    else:
        return [["2",[Components1[0]*Components2[1],-Components1[0]*Components2[0]]]]

def ProductGamma413(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerGamma4(rep1)
    Symbol2 = SymbolToIntegerGamma4(rep2)
    Symbol = IntegerToSymbolGamma4((Symbol1+Symbol2)%2)
    return [["3"+Symbol,[Components1[0]*Components2[k] for k in range(3)]]]

def ProductGamma422(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Comp1 = [(Components1[0]*Components2[0]+Components1[1]*Components2[1])]
    Comp1Prime = [(Components1[0]*Components2[1]-Components1[1]*Components2[0])]
    Comp2 = [(Components1[1]*Components2[1]-Components1[0]*Components2[0]),(Components1[0]*Components2[1]+Components1[1]*Components2[0])]
    return [["1",Comp1],["1'",Comp1Prime],["2",Comp2]]
    
def ProductGamma423(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol2 = SymbolToIntegerGamma4(rep2)
    if Symbol2 == 0:
        Comp3 = [2*Components1[0]*Components2[0],3**0.5*Components1[1]*Components2[2]-Components1[0]*Components2[1],3**0.5*Components1[1]*Components2[1]-Components1[0]*Components2[2]]
        Comp3Prime = [-2*Components1[1]*Components2[0],3**0.5*Components1[0]*Components2[2]+Components1[1]*Components2[1],3**0.5*Components1[0]*Components2[1]+Components1[1]*Components2[2]]
        return [["3", Comp3],["3'",Comp3Prime]]
    if Symbol2 == 1:
        Comp3Prime = [2*Components1[0]*Components2[0],3**0.5*Components1[1]*Components2[2]-Components1[0]*Components2[1],3**0.5*Components1[1]*Components2[1]-Components1[0]*Components2[2]]
        Comp3 = [-2*Components1[1]*Components2[0],3**0.5*Components1[0]*Components2[2]+Components1[1]*Components2[1],3**0.5*Components1[0]*Components2[1]+Components1[1]*Components2[2]]
        return [["3", Comp3],["3'",Comp3Prime]]
    
def ProductGamma433(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerGamma4(rep1)
    Symbol2 = SymbolToIntegerGamma4(rep2)
    SymbolRep1 = IntegerToSymbolGamma4((Symbol1+Symbol2)%2)
    p = (Symbol1+Symbol2)%2
    if p == 0:
        Comp1 = [Components1[0]*Components2[0]+Components1[1]*Components2[2]+Components1[2]*Components2[1]]
        Comp2 = [2*Components1[0]*Components2[0]-Components1[1]*Components2[2]-Components1[2]*Components2[1],3**0.5*Components1[1]*Components2[1]+3**0.5*Components1[2]*Components2[2]]
        Comp3 = [Components1[2]*Components2[2]-Components1[1]*Components2[1],Components1[0]*Components2[2]+Components1[2]*Components2[0],-Components1[0]*Components2[1]-Components1[1]*Components2[0]]
        Comp3Prime = [Components1[2]*Components2[1]-Components1[1]*Components2[2],Components1[1]*Components2[0]-Components1[0]*Components2[1],Components1[0]*Components2[2]-Components1[2]*Components2[0]]
        return [["1"+SymbolRep1,Comp1],["2",Comp2],["3",Comp3],["3'", Comp3Prime]]
    if p == 1:
        Comp1 = [Components1[0]*Components2[0]+Components1[1]*Components2[2]+Components1[2]*Components2[1]]
        Comp2 = [3**0.5*Components1[1]*Components2[1]+3**0.5*Components1[2]*Components2[2],-2*Components1[0]*Components2[0]+Components1[1]*Components2[2]+Components1[2]*Components2[1]]
        Comp3Prime = [Components1[2]*Components2[2]-Components1[1]*Components2[1],Components1[0]*Components2[2]+Components1[2]*Components2[0],-Components1[0]*Components2[1]-Components1[1]*Components2[0]]
        Comp3 = [Components1[2]*Components2[1]-Components1[1]*Components2[2],Components1[1]*Components2[0]-Components1[0]*Components2[1],Components1[0]*Components2[2]-Components1[2]*Components2[0]]
        return [["1"+SymbolRep1,Comp1],["2",Comp2],["3",Comp3],["3'", Comp3Prime]]

def ProductGamma4(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet2 does not match the dimension of its representation")  
    if dim1 == 1 and dim2 == 1:
        return ProductGamma411(Multiplet1,Multiplet2)
    if dim1 == 1 and dim2 == 2:
        return ProductGamma412(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 1:
        return ProductGamma412(Multiplet2,Multiplet1)
    if dim1 == 1 and dim2 == 3:
        return ProductGamma413(Multiplet1, Multiplet2)
    if dim1 == 3 and dim2 == 1:
        return ProductGamma413(Multiplet2, Multiplet1)
    if dim1 == 2 and dim2 == 2:
        return ProductGamma422(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 3:
        return ProductGamma423(Multiplet1, Multiplet2)
    if dim1 == 3 and dim2 == 2:
        return ProductGamma423(Multiplet2,Multiplet1)
    if dim1 == 3 and dim2 == 3:
        return ProductGamma433(Multiplet1, Multiplet2)