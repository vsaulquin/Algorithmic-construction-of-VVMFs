# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================
import numpy as np


#==========================    
#T'xZ4 Clebsch-Gordan table
#==========================
 
def SymbolToIntegerTZ4(rep):
    if len(rep) < 4:
        return 0
    elif rep[-1] == "'":
        return 1
    else:
        return 2

def IntegerToSymbolTZ4(integer):
    if integer == 0:
        return ""
    elif integer == 1:
        return "'"
    else:
        return '"'

def ProductTZ411(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    copy1 = int(rep1[2])
    copy2 = int(rep2[2])
    copy = str((copy1+copy2)%4)
    Symbol1 = SymbolToIntegerTZ4(rep1)
    Symbol2 = SymbolToIntegerTZ4(rep2)
    Symbol = (Symbol1+Symbol2)%3
    Symbol = IntegerToSymbolTZ4(Symbol)
    Component = [Components1[0]*Components2[0]]
    rep = "1_"+copy+Symbol
    return [[rep, Component]]

def ProductTZ412(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerTZ4(rep1)
    Symbol2 = SymbolToIntegerTZ4(rep2)
    Symbol = (Symbol1+Symbol2)%3
    Symbol = IntegerToSymbolTZ4(Symbol)
    copy1 = int(rep1[2])
    copy2 = int(rep2[2])
    copy = str((copy1+copy2)%4)
    return [["2_"+copy+Symbol,[Components1[0]*Components2[0],Components1[0]*Components2[1]]]]

def ProductTZ413(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol = SymbolToIntegerTZ4(rep1)
    copy1 = int(rep1[2])
    copy2 = int(rep2[2])
    copy = str((copy1+copy2)%4)
    rep = "3_"+copy
    Components = []
    for i in range(3):
        Components.append(Components1[0]*Components2[(i-Symbol)%3])
    return [[rep,Components]]

def ProductTZ422(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerTZ4(rep1)
    Symbol2 = SymbolToIntegerTZ4(rep2)
    Symbol = (Symbol1+Symbol2)%3
    copy1 = int(rep1[2])
    copy2 = int(rep2[2])
    copy = str((copy1+copy2)%4)
    ComponentSinglet = [Components1[0]*Components2[1]-Components1[1]*Components2[0]]
    ComponentsTriplet1 = 2**0.5*np.exp(1j*np.pi/6)*Components1[1]*Components2[1]
    ComponentsTriplet2 = np.exp(7j*np.pi/12)*(Components1[0]*Components2[1]+Components1[1]*Components2[0])
    ComponentsTriplet3 = 2**0.5*Components1[0]*Components2[0]
    if Symbol == 0:
        return [["1_"+copy+"'",ComponentSinglet],["3_"+copy,[ComponentsTriplet1, ComponentsTriplet2, ComponentsTriplet3]]]
    elif Symbol == 1:
        return [['1_'+copy+'"',ComponentSinglet],["3_"+copy,[ComponentsTriplet3, ComponentsTriplet1, ComponentsTriplet2]]]
    else:
        return [["1_"+copy,ComponentSinglet],["3_"+copy,[ComponentsTriplet2, ComponentsTriplet3, ComponentsTriplet1]]]
    
def ProductTZ423(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerTZ4(rep1)
    copy1 = int(rep1[2])
    copy2 = int(rep2[2])
    copy = str((copy1+copy2)%4)
    ComponentsFirstDoublet = [Components1[0]*Components2[0]-2**0.5*np.exp(7j*np.pi/12)*Components1[1]*Components2[1],-Components1[1]*Components2[0]+2**0.5*np.exp(5j*np.pi/12)*Components1[0]*Components2[2]]
    ComponentsSecondDoublet = [Components1[0]*Components2[1]-2**0.5*np.exp(7j*np.pi/12)*Components1[1]*Components2[2],-Components1[1]*Components2[1]+2**0.5*np.exp(5j*np.pi/12)*Components1[0]*Components2[0]]
    ComponentsThirdDoublet = [Components1[0]*Components2[2]-2**0.5*np.exp(7j*np.pi/12)*Components1[1]*Components2[0],-Components1[1]*Components2[2]+2**0.5*np.exp(5j*np.pi/12)*Components1[0]*Components2[1]]
    if Symbol1 == 0:
        return [["2_"+copy, ComponentsFirstDoublet],["2_"+copy+"'",ComponentsSecondDoublet],['2_'+copy+'"',ComponentsThirdDoublet]]
    elif Symbol1 == 1:
        return [["2_"+copy, ComponentsThirdDoublet],["2_"+copy+"'",ComponentsFirstDoublet],['2_'+copy+'"',ComponentsSecondDoublet]]
    else:
        return [["2_"+copy, ComponentsSecondDoublet],["2_"+copy+"'",ComponentsThirdDoublet],['2_'+copy+'"',ComponentsFirstDoublet]]
    
def ProductTZ433(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    copy1 = int(rep1[2])
    copy2 = int(rep2[2])
    copy = str((copy1+copy2)%4)
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
    
def ProductTZ4(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet2 does not match the dimension of its representation")  
    if dim1 == 1 and dim2 == 1:
        return ProductTZ411(Multiplet1,Multiplet2)
    if dim1 == 1 and dim2 == 2:
        return ProductTZ412(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 1:
        return ProductTZ412(Multiplet2,Multiplet1)
    if dim1 == 1 and dim2 == 3:
        return ProductTZ413(Multiplet1, Multiplet2)
    if dim1 == 3 and dim2 == 1:
        return ProductTZ413(Multiplet2, Multiplet1)
    if dim1 == 2 and dim2 == 2:
        return ProductTZ422(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 3:
        return ProductTZ423(Multiplet1, Multiplet2)
    if dim1 == 3 and dim2 == 2:
        return ProductTZ423(Multiplet2,Multiplet1)
    if dim1 == 3 and dim2 == 3:
        return ProductTZ433(Multiplet1, Multiplet2)
