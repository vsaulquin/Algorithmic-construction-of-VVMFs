# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================



#=======================================
#DoubleGamma4 = S4' Clebsch-Gordan table
#=======================================

def SymbolToIntegerDoubleGamma4(rep):
    if rep[0] == "1" or rep[0]=="3":
        if len(rep) == 1:
            return 0
        if len(rep) == 3:
            return 3
        if rep[1] == "'":
            return 2
        else:
            return 1
    if rep[0] == "2":
        return len(rep)-1

def IntegerToSymbolDoubleGamma4(Integer, dim):
    if dim == 1 or dim == 3:
        if Integer == 0:
            return ""
        if Integer == 1:
            return "^"
        if Integer == 2:
            return "'"
        else:
            return "^'"
    if dim == 2:
        if Integer == 0:
            return ""
        else:
            return "^"

def ProductDoubleGamma411(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerDoubleGamma4(rep1)
    Symbol2 = SymbolToIntegerDoubleGamma4(rep2)
    Symbol = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%4,1)
    return [["1"+Symbol,[Components1[0]*Components2[0]]]]

def ProductDoubleGamma412(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerDoubleGamma4(rep1)
    Symbol2 = SymbolToIntegerDoubleGamma4(rep2)
    Symbol = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%2,2)
    p = (Symbol1+Symbol2)%4
    if p == 0 or p == 1:
        return [["2"+Symbol,[Components1[0]*Components2[0],Components1[0]*Components2[1]]]]
    else:
        return [["2"+Symbol,[Components1[0]*Components2[1],-Components1[0]*Components2[0]]]]

def ProductDoubleGamma413(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerDoubleGamma4(rep1)
    Symbol2 = SymbolToIntegerDoubleGamma4(rep2)
    Symbol = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%4,3)
    return [["3"+Symbol,[Components1[0]*Components2[k] for k in range(3)]]]

def ProductDoubleGamma422(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerDoubleGamma4(rep1)
    Symbol2 = SymbolToIntegerDoubleGamma4(rep2)
    Symbol = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%2,2)
    SymbolPrime = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%2+2,1)
    p = Symbol1*Symbol2
    if p == 0:
        Comp1 = [(Components1[0]*Components2[0]+Components1[1]*Components2[1])]
        Comp1Prime = [(Components1[0]*Components2[1]-Components1[1]*Components2[0])]
        Comp2 = [(Components1[1]*Components2[1]-Components1[0]*Components2[0]),(Components1[0]*Components2[1]+Components1[1]*Components2[0])]
        return [["1"+Symbol,Comp1],["1"+SymbolPrime,Comp1Prime],["2"+Symbol,Comp2]]
    else:
        Comp1 = [(Components1[0]*Components2[1]-Components1[1]*Components2[0])]
        Comp1Prime = [(Components1[0]*Components2[0]+Components1[1]*Components2[1])]
        Comp2 =  [(Components1[0]*Components2[1]+Components1[1]*Components2[0]),(Components1[0]*Components2[0]-Components1[1]*Components2[1])]
        return [["1"+Symbol,Comp1],["1"+SymbolPrime,Comp1Prime],["2"+Symbol,Comp2]]
    
def ProductDoubleGamma423(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerDoubleGamma4(rep1)
    Symbol2 = SymbolToIntegerDoubleGamma4(rep2)
    Symbol = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%2,3)
    SymbolPrime = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%2+2,3)
    p = (Symbol1+Symbol2)%4
    if p == 0 or p == 1:
        Comp3 = [2*Components1[0]*Components2[0],3**0.5*Components1[1]*Components2[2]-Components1[0]*Components2[1],3**0.5*Components1[1]*Components2[1]-Components1[0]*Components2[2]]
        Comp3Prime = [-2*Components1[1]*Components2[0],3**0.5*Components1[0]*Components2[2]+Components1[1]*Components2[1],3**0.5*Components1[0]*Components2[1]+Components1[1]*Components2[2]]
        return [["3"+Symbol, Comp3],["3"+SymbolPrime,Comp3Prime]]
    if p == 2 or p == 3:
        Comp3Prime = [2*Components1[0]*Components2[0],3**0.5*Components1[1]*Components2[2]-Components1[0]*Components2[1],3**0.5*Components1[1]*Components2[1]-Components1[0]*Components2[2]]
        Comp3 = [-2*Components1[1]*Components2[0],3**0.5*Components1[0]*Components2[2]+Components1[1]*Components2[1],3**0.5*Components1[0]*Components2[1]+Components1[1]*Components2[2]]
        return [["3"+Symbol, Comp3],["3"+SymbolPrime,Comp3Prime]]
    
def ProductDoubleGamma433(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerDoubleGamma4(rep1)
    Symbol2 = SymbolToIntegerDoubleGamma4(rep2)
    SymbolRep1 = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%4, 1)
    SymbolRep2 = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%2, 2)
    SymbolRep3 = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%2, 3)
    SymbolRep3Prime = IntegerToSymbolDoubleGamma4((Symbol1+Symbol2)%2+2, 3)
    p = (Symbol1+Symbol2)%4
    if p == 0 or p == 1:
        Comp1 = [Components1[0]*Components2[0]+Components1[1]*Components2[2]+Components1[2]*Components2[1]]
        Comp2 = [2*Components1[0]*Components2[0]-Components1[1]*Components2[2]-Components1[2]*Components2[1],3**0.5*Components1[1]*Components2[1]+3**0.5*Components1[2]*Components2[2]]
        Comp3 = [Components1[2]*Components2[2]-Components1[1]*Components2[1],Components1[0]*Components2[2]+Components1[2]*Components2[0],-Components1[0]*Components2[1]-Components1[1]*Components2[0]]
        Comp3Prime = [Components1[2]*Components2[1]-Components1[1]*Components2[2],Components1[1]*Components2[0]-Components1[0]*Components2[1],Components1[0]*Components2[2]-Components1[2]*Components2[0]]
        return [["1"+SymbolRep1,Comp1],["2"+SymbolRep2,Comp2],["3"+SymbolRep3,Comp3],["3"+SymbolRep3Prime, Comp3Prime]]
    if p == 2 or p == 3:
        Comp1 = [Components1[0]*Components2[0]+Components1[1]*Components2[2]+Components1[2]*Components2[1]]
        Comp2 = [3**0.5*Components1[1]*Components2[1]+3**0.5*Components1[2]*Components2[2],-2*Components1[0]*Components2[0]+Components1[1]*Components2[2]+Components1[2]*Components2[1]]
        Comp3Prime = [Components1[2]*Components2[2]-Components1[1]*Components2[1],Components1[0]*Components2[2]+Components1[2]*Components2[0],-Components1[0]*Components2[1]-Components1[1]*Components2[0]]
        Comp3 = [Components1[2]*Components2[1]-Components1[1]*Components2[2],Components1[1]*Components2[0]-Components1[0]*Components2[1],Components1[0]*Components2[2]-Components1[2]*Components2[0]]
        return [["1"+SymbolRep1,Comp1],["2"+SymbolRep2,Comp2],["3"+SymbolRep3,Comp3],["3"+SymbolRep3Prime, Comp3Prime]]

def ProductDoubleGamma4(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet2 does not match the dimension of its representation")  
    if dim1 == 1 and dim2 == 1:
        return ProductDoubleGamma411(Multiplet1,Multiplet2)
    if dim1 == 1 and dim2 == 2:
        return ProductDoubleGamma412(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 1:
        return ProductDoubleGamma412(Multiplet2,Multiplet1)
    if dim1 == 1 and dim2 == 3:
        return ProductDoubleGamma413(Multiplet1, Multiplet2)
    if dim1 == 3 and dim2 == 1:
        return ProductDoubleGamma413(Multiplet2, Multiplet1)
    if dim1 == 2 and dim2 == 2:
        return ProductDoubleGamma422(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 3:
        return ProductDoubleGamma423(Multiplet1, Multiplet2)
    if dim1 == 3 and dim2 == 2:
        return ProductDoubleGamma423(Multiplet2,Multiplet1)
    if dim1 == 3 and dim2 == 3:
        return ProductDoubleGamma433(Multiplet1, Multiplet2)