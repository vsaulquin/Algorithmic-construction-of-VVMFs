# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================



#==========================
#S4xZ3 Clebsch-Gordan table
#==========================

def SymbolToIntegerS4Z3(rep):
    if rep[0] == "1" or rep[0]=="3":
        a = int(rep[2])
        if len(rep) == 3:
            r = 0
        else:
            r = 1
    if rep[0] == "2":
        a = int(rep[2])
        r = None
    return (a,r)

def IntegerToSymbolS4Z3(a, r, dim):
    rep = str(dim)
    rep += "_"+str(a)
    if r == 1:
        rep +="'"
    return rep

def ProductS4Z311(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    a1,r1 = SymbolToIntegerS4Z3(rep1)
    a2,r2 = SymbolToIntegerS4Z3(rep2)
    a = (a1+a2)%3
    r = (r1+r2)%2
    rep = IntegerToSymbolS4Z3(a, r, 1)
    return [[rep,[Components1[0]*Components2[0]]]]

def ProductS4Z312(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    a1,r1 = SymbolToIntegerS4Z3(rep1)
    a2,r2 = SymbolToIntegerS4Z3(rep2)
    a = (a1+a2)%3
    rep = IntegerToSymbolS4Z3(a, None, 2)
    if r1 == 0:
        return [[rep,[Components1[0]*Components2[0],Components1[0]*Components2[1]]]]
    if r1 == 1:
        return [[rep,[Components1[0]*Components2[1],-Components1[0]*Components2[0]]]]
    
def ProductS4Z313(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    a1,r1 = SymbolToIntegerS4Z3(rep1)
    a2,r2 = SymbolToIntegerS4Z3(rep2)
    a = (a1+a2)%3
    r = (r1+r2)%2
    rep = IntegerToSymbolS4Z3(a, r, 3)
    return [[rep,[Components1[0]*Components2[0], Components1[0]*Components2[1],Components1[0]*Components2[2]]]]

def ProductS4Z322(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    a1,r1 = SymbolToIntegerS4Z3(rep1)
    a2,r2 = SymbolToIntegerS4Z3(rep2)
    a = (a1+a2)%3
    Comp1 = [(Components1[0]*Components2[0]+Components1[1]*Components2[1])]
    Comp1Prime = [(Components1[0]*Components2[1]-Components1[1]*Components2[0])]
    Comp2 = [(Components1[1]*Components2[1]-Components1[0]*Components2[0]),(Components1[0]*Components2[1]+Components1[1]*Components2[0])]
    return  [[IntegerToSymbolS4Z3(a,0,1),Comp1],[IntegerToSymbolS4Z3(a,1,1),Comp1Prime],[IntegerToSymbolS4Z3(a,0,2),Comp2]]

def ProductS4Z323(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    a1,r1 = SymbolToIntegerS4Z3(rep1)
    a2,r2 = SymbolToIntegerS4Z3(rep2)
    a = (a1+a2)%3
    Comp3 = [2*Components1[0]*Components2[0],3**0.5*Components1[1]*Components2[2]-Components1[0]*Components2[1],3**0.5*Components1[1]*Components2[1]-Components1[0]*Components2[2]]
    Comp3Prime = [-2*Components1[1]*Components2[0],3**0.5*Components1[0]*Components2[2]+Components1[1]*Components2[1],3**0.5*Components1[0]*Components2[1]+Components1[1]*Components2[2]]
    if r2 == 0:
        return  [[IntegerToSymbolS4Z3(a,0,3),Comp3],[IntegerToSymbolS4Z3(a,1,3),Comp3Prime]]
    if r2 == 1:
            return  [[IntegerToSymbolS4Z3(a,0,3),Comp3Prime],[IntegerToSymbolS4Z3(a,1,3),Comp3]]
    
def ProductS4Z333(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    a1,r1 = SymbolToIntegerS4Z3(rep1)
    a2,r2 = SymbolToIntegerS4Z3(rep2)
    a = (a1+a2)%3
    r = (r1+r2)%2
    if r == 0:
        Comp1 = [(Components1[0]*Components2[0]+Components1[1]*Components2[2]+Components1[2]*Components2[1])]
        Comp2 = [(2*Components1[0]*Components2[0]-Components1[1]*Components2[2]-Components1[2]*Components2[1]),3**0.5*(Components1[1]*Components2[1]+Components1[2]*Components2[2])]
        Comp3 = [(Components1[2]*Components2[2]-Components1[1]*Components2[1]),(Components1[0]*Components2[2]+Components1[2]*Components2[0]),(-Components1[0]*Components2[1]-Components1[1]*Components2[0])]
        Comp3Prime = [(Components1[2]*Components2[1]-Components1[1]*Components2[2]),(Components1[1]*Components2[0]-Components1[0]*Components2[1]),(Components1[0]*Components2[2]-Components1[2]*Components2[0])]
        return  [[IntegerToSymbolS4Z3(a,r,1),Comp1],[IntegerToSymbolS4Z3(a,0,2),Comp2],[IntegerToSymbolS4Z3(a,0,3),Comp3],[IntegerToSymbolS4Z3(a,1,3),Comp3Prime]]
    if r == 1:
        Comp1 = [(Components1[0]*Components2[0]+Components1[1]*Components2[2]+Components1[2]*Components2[1])]
        Comp2 = [3**0.5*(Components1[1]*Components2[1]+Components1[2]*Components2[2]),-(2*Components1[0]*Components2[0]-Components1[1]*Components2[2]-Components1[2]*Components2[1])]
        Comp3Prime = [(Components1[2]*Components2[2]-Components1[1]*Components2[1]),(Components1[0]*Components2[2]+Components1[2]*Components2[0]),(-Components1[0]*Components2[1]-Components1[1]*Components2[0])]
        Comp3 = [(Components1[2]*Components2[1]-Components1[1]*Components2[2]),(Components1[1]*Components2[0]-Components1[0]*Components2[1]),(Components1[0]*Components2[2]-Components1[2]*Components2[0])]
        return  [[IntegerToSymbolS4Z3(a,r,1),Comp1],[IntegerToSymbolS4Z3(a,0,2),Comp2],[IntegerToSymbolS4Z3(a,0,3),Comp3],[IntegerToSymbolS4Z3(a,1,3),Comp3Prime]]
    
def ProductS4Z3(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet2 does not match the dimension of its representation")  
    if dim1 == 1 and dim2 == 1:
        return ProductS4Z311(Multiplet1, Multiplet2)
    if dim1 == 1 and dim2 == 2:
        return ProductS4Z312(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 1:
        return ProductS4Z312(Multiplet2,Multiplet1)
    if dim1 == 1 and dim2 == 3:
        return ProductS4Z313(Multiplet1,Multiplet2)
    if dim1 == 3 and dim2 == 1:
        return ProductS4Z313(Multiplet2,Multiplet1)
    if dim1 == 2 and dim2 == 2:
        return ProductS4Z322(Multiplet1, Multiplet2)
    if dim1 == 2 and dim2 == 3:
        return ProductS4Z323(Multiplet1, Multiplet2)
    if dim1 == 3 and dim2 == 2:
        return ProductS4Z323(Multiplet2, Multiplet1)
    if dim1 == 3 and dim2 == 3:
        return ProductS4Z333(Multiplet1, Multiplet2)