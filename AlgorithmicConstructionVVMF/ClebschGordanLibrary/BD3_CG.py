# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================



#====================================
#Binary Dihedral Clebsch-Gordan table
#====================================

def SymbolToIntegerBD3(rep):
    if rep[0] == "1":
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

def IntegerToSymbolBD3(Integer, dim):
    if dim == 1:
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

def ProductBD311(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerBD3(rep1)
    Symbol2 = SymbolToIntegerBD3(rep2)
    Symbol = IntegerToSymbolBD3((Symbol1+Symbol2)%4,1)
    return [["1"+Symbol,[Components1[0]*Components2[0]]]]

def ProductBD312(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerBD3(rep1)
    Symbol2 = SymbolToIntegerBD3(rep2)
    Symbol = IntegerToSymbolBD3((Symbol1+Symbol2)%2,2)
    p = (Symbol1+Symbol2)%4
    if p == 0 or p == 1:
        return [["2"+Symbol,[Components1[0]*Components2[0],Components1[0]*Components2[1]]]]
    else:
        return [["2"+Symbol,[Components1[0]*Components2[1],-Components1[0]*Components2[0]]]]


def ProductBD322(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerBD3(rep1)
    Symbol2 = SymbolToIntegerBD3(rep2)
    Symbol = IntegerToSymbolBD3((Symbol1+Symbol2)%2,2)
    SymbolPrime = IntegerToSymbolBD3((Symbol1+Symbol2)%2+2,1)
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
    

def ProductBD3(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet2 does not match the dimension of its representation")  
    if dim1 == 1 and dim2 == 1:
        return ProductBD311(Multiplet1,Multiplet2)
    if dim1 == 1 and dim2 == 2:
        return ProductBD312(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 1:
        return ProductBD312(Multiplet2,Multiplet1)
    if dim1 == 2 and dim2 == 2:
        return ProductBD322(Multiplet1,Multiplet2)
