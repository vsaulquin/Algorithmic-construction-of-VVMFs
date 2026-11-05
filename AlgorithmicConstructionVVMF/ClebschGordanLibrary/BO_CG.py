# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================



#====================================================
#Double Binary Octahedral Clebsch-Gordan Coefficients
#====================================================   

def SymbolToIntegerBO(rep):
    if rep[0] == "1":
        return len(rep)-1
    elif rep[0] == "2":
        return len(rep)-1
    elif rep[0] == "3":
        return len(rep)-1
    else:
        return 0

def IntegerToSymbolBO(Integer,dim):
    if dim == 1:
        if Integer == 0:
            return ""
        else:
            return "'"
    if dim == 2:
        if Integer == 0:
            return ""
        elif Integer == 1:
            return "^"
        else:
            return "^'"
    if dim == 3:
        if Integer == 0:
            return ""
        else:
            return "'"
    if dim == 4:
        return "^"
    
def Product1NBO(Multiplet1,MultipletN):
    rep1, Components1 = Multiplet1
    repN, ComponentsN = MultipletN
    if SymbolToIntegerBO(rep1) == 0:
        return [[repN,[Components1[0]*ComponentsN[k] for k in range(len(ComponentsN))]]]
    else:
        if repN[0] == "1":
            Integer = (SymbolToIntegerBO(repN)+1)%2
            Symbol = IntegerToSymbolBO(Integer,1)
            return [["1"+Symbol,[Components1[0]*ComponentsN[0]]]]
        elif repN[0] == "2":
            Integer = -SymbolToIntegerBO(repN)%3
            Symbol = IntegerToSymbolBO(Integer,2)
            return [["2"+Symbol, [Components1[0]*ComponentsN[1],-Components1[0]*ComponentsN[0]]]]
        elif repN[0] == "3":
            Integer = (SymbolToIntegerBO(repN)+1)%2
            Symbol = IntegerToSymbolBO(Integer, 3)
            return [["3"+Symbol,[Components1[0]*ComponentsN[k] for k in range(3)]]]
        elif repN[0] == "4":
            return [["4^",[Components1[0]*ComponentsN[3],-Components1[0]*ComponentsN[2],Components1[0]*ComponentsN[1],-Components1[0]*ComponentsN[0]]]]
        
def Product22BO(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerBO(rep1)
    Symbol2 = SymbolToIntegerBO(rep2)
    if Symbol1 == 0 and Symbol2 == 0:
        Comp1 = [Components1[0]*Components2[0]+Components1[1]*Components2[1]]
        Comp1Prime = [Components1[1]*Components2[0]-Components1[0]*Components2[1]]
        Comp2 = [Components1[1]*Components2[1]-Components1[0]*Components2[0],Components1[1]*Components2[0]+Components1[0]*Components2[1]]
        return [["1",Comp1],["1'",Comp1Prime],["2",Comp2]]
    if Symbol1 == 0 and Symbol2 == 1:
        return [["4^",[-Components1[0]*Components2[1],Components1[0]*Components2[0],Components1[1]*Components2[0],Components1[1]*Components2[1]]]]
    if Symbol1 == 1 and Symbol2 == 0:
        return Product22BO(Multiplet2,Multiplet1)
    if Symbol1 == 0 and Symbol2 == 2:
        return [["4^",[Components1[1]*Components2[0],Components1[1]*Components2[1],-Components1[0]*Components2[1],Components1[0]*Components2[0]]]]
    if Symbol1 == 2 and Symbol2 == 0:
        return Product22BO(Multiplet2, Multiplet1)
    if Symbol1 == 1 and Symbol2 == 1:
        Comp1 = [Components1[1]*Components2[0]-Components1[0]*Components2[1]]
        Comp3 = [Components1[1]*Components2[0]+Components1[0]*Components2[1],2**0.5*Components1[1]*Components2[1],-2**0.5*Components1[0]*Components2[0]]
        return [["1",Comp1],["3",Comp3]]
    if Symbol1 == 1 and Symbol2 == 2:
        Comp1Prime = [Components1[0]*Components2[0]+Components1[1]*Components2[1]]
        Comp3Prime = [Components1[0]*Components2[0]-Components1[1]*Components2[1],2**0.5*Components1[1]*Components2[0],2**0.5*Components1[0]*Components2[1]]
        return [["1'",Comp1Prime],["3'",Comp3Prime]]
    if Symbol1 == 2 and Symbol2 == 1:
        return Product22BO(Multiplet2,Multiplet1)
    if Symbol1 == 2 and Symbol2 == 2:
        Comp1 = [Components1[1]*Components2[0]-Components1[0]*Components2[1]]
        Comp3 = [Components1[1]*Components2[0]+Components1[0]*Components2[1],-2**0.5*Components1[0]*Components2[0],2**0.5*Components1[1]*Components2[1]]
        return [["1",Comp1],["3",Comp3]]

def Product23BO(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerBO(rep1)
    Symbol2 = SymbolToIntegerBO(rep2)
    if Symbol1 == 0:
        Triplet = [-2*Components1[0]*Components2[0],Components1[0]*Components2[1]+3**0.5*Components1[1]*Components2[2],3**0.5*Components1[1]*Components2[1]+Components1[0]*Components2[2]]
        Triplet2 = [-2*Components1[1]*Components2[0],Components1[1]*Components2[1]-3**0.5*Components1[0]*Components2[2],Components1[1]*Components2[2]-3**0.5*Components1[0]*Components2[1]]
        if Symbol2 == 0:
            return [["3",Triplet],["3'",Triplet2]]
        if Symbol2 == 1:
            return [["3",Triplet2],["3'",Triplet]]
    if Symbol1 == 1:
        if Symbol2 == 0:
            Comp2 = [Components1[0]*Components2[0]+2**0.5*Components1[1]*Components2[2],2**0.5*Components1[0]*Components2[1]-Components1[1]*Components2[0]]
            Comp4 = [-2**0.5*Components1[1]*Components2[0]-Components1[0]*Components2[1],Components1[1]*Components2[2]-2**0.5*Components1[0]*Components2[0],3**0.5*Components1[1]*Components2[1],3**0.5*Components1[0]*Components2[2]]
            return [["2^",Comp2],["4^",Comp4]]
        if Symbol2 == 1:
            Comp2 = [Components1[1]*Components2[0]-2**0.5*Components1[0]*Components2[1],Components1[0]*Components2[0]+2**0.5*Components1[1]*Components2[2]]
            Comp4 = [3**0.5*Components1[0]*Components2[2],-3**0.5*Components1[1]*Components2[1],Components1[1]*Components2[2]-2**0.5*Components1[0]*Components2[0],2**0.5*Components1[1]*Components2[0]+Components1[0]*Components2[1]]
            return [["2^'",Comp2],["4^",Comp4]]
    if Symbol1 == 2:
        if Symbol2 == 0: 
            Comp2 = [2**0.5*Components1[1]*Components2[1]+Components1[0]*Components2[0],-Components1[1]*Components2[0]+2**0.5*Components1[0]*Components2[2]]
            Comp4 = [3**0.5*Components1[1]*Components2[2],3**0.5*Components1[0]*Components2[1],-2**0.5*Components1[1]*Components2[0]-Components1[0]*Components2[2],-2**0.5*Components1[0]*Components2[0]+Components1[1]*Components2[1]]
            return [["2^'",Comp2],["4^",Comp4]]
        if Symbol2 == 1:
            Comp2 = [Components1[1]*Components2[0]-2**0.5*Components1[0]*Components2[2],Components1[0]*Components2[0]+2**0.5*Components1[1]*Components2[1]]
            Comp4 = [Components1[1]*Components2[1]-2**0.5*Components1[0]*Components2[0],2**0.5*Components1[1]*Components2[0]+Components1[0]*Components2[2],3**0.5*Components1[0]*Components2[1],-3**0.5*Components1[1]*Components2[2]]
            return [["2^",Comp2],["4^",Comp4]]
        

def Product24BO(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerBO(rep1)
    if Symbol1 == 0:
        Comp2 = [Components1[0]*Components2[1]+Components1[1]*Components2[2],Components1[1]*Components2[3]-Components1[0]*Components2[0]]
        Comp2Prime = [Components1[1]*Components2[0]+Components1[0]*Components2[3],Components1[1]*Components2[1]-Components1[0]*Components2[2]]
        Comp4 = [Components1[0]*Components2[0]+Components1[1]*Components2[3],Components1[0]*Components2[1]-Components1[1]*Components2[2],-Components1[1]*Components2[1]-Components1[0]*Components2[2],Components1[1]*Components2[0]-Components1[0]*Components2[3]]
        return [["2^",Comp2],["2^'",Comp2Prime],["4^",Comp4]]
    if Symbol1 == 1:
        Comp2 = [Components1[0]*Components2[0]+Components1[1]*Components2[1],Components1[1]*Components2[2]-Components1[0]*Components2[3]]
        Comp3 = [2**0.5*Components1[1]*Components2[1]-2**0.5*Components1[0]*Components2[0],Components1[1]*Components2[0]+3**0.5*Components1[0]*Components2[2],Components1[0]*Components2[1]-3**0.5*Components1[1]*Components2[3]]
        Comp3Prime = [-2**0.5*Components1[1]*Components2[2]-2**0.5*Components1[0]*Components2[3],3**0.5*Components1[0]*Components2[1]+Components1[1]*Components2[3],3**0.5*Components1[1]*Components2[0]-Components1[0]*Components2[2]]
        return [["2",Comp2],["3",Comp3],["3'",Comp3Prime]]
    if Symbol1 == 2:
        Comp2 = [Components1[0]*Components2[2]+Components1[1]*Components2[3],Components1[1]*Components2[0]-Components1[0]*Components2[1]]
        Comp3 = [2**0.5*Components1[1]*Components2[3]-2**0.5*Components1[0]*Components2[2],Components1[0]*Components2[3]-3**0.5*Components1[1]*Components2[1],3**0.5*Components1[0]*Components2[0]+Components1[1]*Components2[2]]
        Comp3Prime = [2**0.5*Components1[1]*Components2[0]+2**0.5*Components1[0]*Components2[1],Components1[0]*Components2[0]-3**0.5*Components1[1]*Components2[2],-Components1[1]*Components2[1]-3**0.5*Components1[0]*Components2[3]]
        return [["2",Comp2],["3",Comp3],["3'",Comp3Prime]]
    
def Product33BO(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerBO(rep1)
    Symbol2 = SymbolToIntegerBO(rep2)
    if Symbol1 == Symbol2:
        Comp1 = [Components1[0]*Components2[0]+Components1[2]*Components2[1]+Components1[1]*Components2[2]]
        Comp2 = [Components1[2]*Components2[1]+Components1[1]*Components2[2]-2*Components1[0]*Components2[0],3**0.5*Components1[1]*Components2[1]+3**0.5*Components1[2]*Components2[2]]
        Comp3 = [Components1[2]*Components2[1]-Components1[1]*Components2[2],Components1[1]*Components2[0]-Components1[0]*Components2[1],Components1[0]*Components2[2]-Components1[2]*Components2[0]]
        Comp3Prime = [Components1[2]*Components2[2]-Components1[1]*Components2[1],Components1[2]*Components2[0]+Components1[0]*Components2[2],-Components1[1]*Components2[0]-Components1[0]*Components2[1]]
        return [["1",Comp1],["2",Comp2],["3",Comp3],["3'",Comp3Prime]]
    if Symbol1 == 0 and Symbol2 == 1:
        Comp1 = [Components1[0]*Components2[0]+Components1[2]*Components2[1]+Components1[1]*Components2[2]]
        Comp2 = [3**0.5*Components1[1]*Components2[1]+3**0.5*Components1[2]*Components2[2],-Components1[2]*Components2[1]-Components1[1]*Components2[2]+2*Components1[0]*Components2[0]]
        Comp3 = [Components1[2]*Components2[2]-Components1[1]*Components2[1],Components1[2]*Components2[0]+Components1[0]*Components2[2],-Components1[1]*Components2[0]-Components1[0]*Components2[1]]
        Comp3Prime = [Components1[2]*Components2[1]-Components1[1]*Components2[2],Components1[1]*Components2[0]-Components1[0]*Components2[1],Components1[0]*Components2[2]-Components1[2]*Components2[0]]
        return [["1'",Comp1],["2",Comp2],["3",Comp3],["3'",Comp3Prime]]
    if Symbol1 == 1 and Symbol2 == 0:
        return Product33BO(Multiplet2,Multiplet1)

def Product34BO(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerBO(rep1)
    if Symbol1 == 0:
        Comp2 = [3**0.5*Components1[1]*Components2[3]-Components1[2]*Components2[0]-2**0.5*Components1[0]*Components2[1],Components1[1]*Components2[1]+3**0.5*Components1[2]*Components2[2]-2**0.5*Components1[0]*Components2[0]]
        Comp2Prime = [3**0.5*Components1[2]*Components2[1]-Components1[1]*Components2[2]-2**0.5*Components1[0]*Components2[3],3**0.5*Components1[1]*Components2[0]+Components1[2]*Components2[3]-2**0.5*Components1[0]*Components2[2]]
        Comp4 = [Components1[0]*Components2[0]+6**0.5*Components1[2]*Components2[2]-2*2**0.5*Components1[1]*Components2[1],-2*2**0.5*Components1[2]*Components2[0]-Components1[0]*Components2[1]-6**0.5*Components1[1]*Components2[3],6**0.5*Components1[1]*Components2[0]+3*Components1[0]*Components2[2],-6**0.5*Components1[2]*Components2[1]-3*Components1[0]*Components2[3]]
        Comp4_2 = [-Components1[0]*Components2[0]-2**0.5*Components1[1]*Components2[1],Components1[0]*Components2[1]-2**0.5*Components1[2]*Components2[0],Components1[0]*Components2[2]+2**0.5*Components1[2]*Components2[3],2**0.5*Components1[1]*Components2[2]-Components1[0]*Components2[3]]
        return [["2^",Comp2],["2^'",Comp2Prime],["4^",Comp4],["4^",Comp4_2]]
    if Symbol1 == 1:
        Comp2 = [3**0.5*Components1[1]*Components2[0]+Components1[2]*Components2[3]-2**0.5*Components1[0]*Components2[2],2**0.5*Components1[0]*Components2[3]+Components1[1]*Components2[2]-3**0.5*Components1[2]*Components2[1]]
        Comp2Prime = [Components1[1]*Components2[1]+3**0.5*Components1[2]*Components2[2]-2**0.5*Components1[0]*Components2[0],Components1[2]*Components2[0]+2**0.5*Components1[0]*Components2[1]-3**0.5*Components1[1]*Components2[3]]
        Comp4 = [Components1[0]*Components2[3]-2**0.5*Components1[1]*Components2[2],Components1[0]*Components2[2]+2**0.5*Components1[2]*Components2[3],2**0.5*Components1[2]*Components2[0]-Components1[0]*Components2[1],-Components1[0]*Components2[0]-2**0.5*Components1[1]*Components2[1]]
        Comp4_2 = [6**0.5*Components1[2]*Components2[1]+3*Components1[0]*Components2[3],6**0.5*Components1[1]*Components2[0]+3*Components1[0]*Components2[2],2*2**0.5*Components1[2]*Components2[0]+Components1[0]*Components2[1]+6**0.5*Components1[1]*Components2[3],Components1[0]*Components2[0]+6**0.5*Components1[2]*Components2[2]-2*2**0.5*Components1[1]*Components2[1]]
        return [["2^",Comp2],["2^'",Comp2Prime],["4^",Comp4],["4^",Comp4_2]]
    
def Product44BO(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Comp1 = [Components1[1]*Components2[0]-Components1[0]*Components2[1]+Components1[3]*Components2[2]-Components1[2]*Components2[3]]
    Comp1Prime = [Components1[2]*Components2[0]+Components1[0]*Components2[2]+Components1[3]*Components2[1]+Components1[1]*Components2[3]]
    Comp2 = [Components1[0]*Components2[1]-Components1[1]*Components2[0]+Components1[3]*Components2[2]-Components1[2]*Components2[3],Components1[2]*Components2[0]-Components1[0]*Components2[2]+Components1[3]*Components2[1]-Components1[1]*Components2[3]]
    Comp3 = [Components1[3]*Components2[2]+Components1[2]*Components2[3]-Components1[1]*Components2[0]-Components1[0]*Components2[1],2**0.5*Components1[0]*Components2[0]+2**0.5*Components1[3]*Components2[3],-2**0.5*Components1[1]*Components2[1]-2**0.5*Components1[2]*Components2[2]]
    Comp3_2 = [2*2**0.5*Components1[1]*Components2[0]+2*2**0.5*Components1[0]*Components2[1],3**0.5*Components1[2]*Components2[1]+3**0.5*Components1[1]*Components2[2]-Components1[0]*Components2[0]-3*Components1[3]*Components2[3],Components1[1]*Components2[1]+3*Components1[2]*Components2[2]+3**0.5*Components1[3]*Components2[0]+3**0.5*Components1[0]*Components2[3]]
    Comp3Prime = [Components1[3]*Components2[1]-Components1[1]*Components2[3]+Components1[0]*Components2[2]-Components1[2]*Components2[0],2**0.5*Components1[0]*Components2[3]-2**0.5*Components1[3]*Components2[0],2**0.5*Components1[1]*Components2[2]-2**0.5*Components1[2]*Components2[1]]
    Comp3Prime_2 = [2**0.5*Components1[3]*Components2[1]+2**0.5*Components1[1]*Components2[3]-2**0.5*Components1[2]*Components2[0]-2**0.5*Components1[0]*Components2[2],3**0.5*Components1[1]*Components2[1]+Components1[3]*Components2[0]+Components1[0]*Components2[3]-3**0.5*Components1[2]*Components2[2],3**0.5*Components1[3]*Components2[3]+Components1[2]*Components2[1]+Components1[1]*Components2[2]-3**0.5*Components1[0]*Components2[0]]
    return [["1",Comp1],["1'",Comp1Prime],["2",Comp2],["3",Comp3],["3",Comp3_2],["3'",Comp3Prime],["3'",Comp3Prime_2]]

def ProductBO(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet2 does not match the dimension of its representation")
    if dim1 == 1:
        return Product1NBO(Multiplet1, Multiplet2)
    if dim2 == 1:
        return Product1NBO(Multiplet2, Multiplet1)
    if dim1 == 2 and dim2 == 2:
        return Product22BO(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 3:
        return Product23BO(Multiplet1,Multiplet2)
    if dim1 == 3 and dim2 == 2:
        return Product23BO(Multiplet2, Multiplet1)
    if dim1 == 2 and dim2 == 4:
        return Product24BO(Multiplet1,Multiplet2)
    if dim1 == 4 and dim2 == 2:
        return Product24BO(Multiplet2,Multiplet1)
    if dim1 == 3 and dim2 == 3:
        return Product33BO(Multiplet1, Multiplet2)
    if dim1 == 3 and dim2 == 4:
        return Product34BO(Multiplet1, Multiplet2)
    if dim1 == 4 and dim2 == 3:
        return Product34BO(Multiplet2, Multiplet1)
    if dim1 == 4 and dim2 == 4:
        return Product44BO(Multiplet1,Multiplet2)