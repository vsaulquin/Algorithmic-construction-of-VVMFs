# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================
import numpy as np


#======================================    
#DoubleGamma3 = T' Clebsch-Gordan table
#======================================

def SymbolToIntegerDoubleGamma3(rep):
    if len(rep) < 2:
        return 0
    elif rep[-1] == "'":
        return 1
    else:
        return 2

def IntegerToSymbolDoubleGamma3(integer):
    if integer == 0:
        return ""
    elif integer == 1:
        return "'"
    else:
        return '"'

def ProductDoubleGamma311(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerDoubleGamma3(rep1)
    Symbol2 = SymbolToIntegerDoubleGamma3(rep2)
    Symbol = (Symbol1+Symbol2)%3
    Symbol = IntegerToSymbolDoubleGamma3(Symbol)
    Component = [Components1[0]*Components2[0]]
    rep = "1"+Symbol
    return [[rep, Component]]

def ProductDoubleGamma312(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerDoubleGamma3(rep1)
    Symbol2 = SymbolToIntegerDoubleGamma3(rep2)
    Symbol = (Symbol1+Symbol2)%3
    Symbol = IntegerToSymbolDoubleGamma3(Symbol)
    return [["2"+Symbol,[Components1[0]*Components2[0],Components1[0]*Components2[1]]]]

def ProductDoubleGamma313(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol = SymbolToIntegerDoubleGamma3(rep1)
    rep = "3"
    Components = []
    for i in range(3):
        Components.append(Components1[0]*Components2[(i-Symbol)%3])
    return [[rep,Components]]

def ProductDoubleGamma322(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerDoubleGamma3(rep1)
    Symbol2 = SymbolToIntegerDoubleGamma3(rep2)
    Symbol = (Symbol1+Symbol2)%3
    ComponentSinglet = [Components1[0]*Components2[1]-Components1[1]*Components2[0]]
    ComponentsTriplet1 = 2**0.5*np.exp(1j*np.pi/6)*Components1[1]*Components2[1]
    ComponentsTriplet2 = np.exp(7j*np.pi/12)*(Components1[0]*Components2[1]+Components1[1]*Components2[0])
    ComponentsTriplet3 = 2**0.5*Components1[0]*Components2[0]
    if Symbol == 0:
        return [["1'",ComponentSinglet],["3",[ComponentsTriplet1, ComponentsTriplet2, ComponentsTriplet3]]]
    elif Symbol == 1:
        return [['1"',ComponentSinglet],["3",[ComponentsTriplet3, ComponentsTriplet1, ComponentsTriplet2]]]
    else:
        return [["1",ComponentSinglet],["3",[ComponentsTriplet2, ComponentsTriplet3, ComponentsTriplet1]]]
    
def ProductDoubleGamma323(Multiplet1, Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    Symbol1 = SymbolToIntegerDoubleGamma3(rep1)
    ComponentsFirstDoublet = [Components1[0]*Components2[0]-2**0.5*np.exp(7j*np.pi/12)*Components1[1]*Components2[1],-Components1[1]*Components2[0]+2**0.5*np.exp(5j*np.pi/12)*Components1[0]*Components2[2]]
    ComponentsSecondDoublet = [Components1[0]*Components2[1]-2**0.5*np.exp(7j*np.pi/12)*Components1[1]*Components2[2],-Components1[1]*Components2[1]+2**0.5*np.exp(5j*np.pi/12)*Components1[0]*Components2[0]]
    ComponentsThirdDoublet = [Components1[0]*Components2[2]-2**0.5*np.exp(7j*np.pi/12)*Components1[1]*Components2[0],-Components1[1]*Components2[2]+2**0.5*np.exp(5j*np.pi/12)*Components1[0]*Components2[1]]
    if Symbol1 == 0:
        return [["2", ComponentsFirstDoublet],["2'",ComponentsSecondDoublet],['2"',ComponentsThirdDoublet]]
    elif Symbol1 == 1:
        return [["2", ComponentsThirdDoublet],["2'",ComponentsFirstDoublet],['2"',ComponentsSecondDoublet]]
    else:
        return [["2", ComponentsSecondDoublet],["2'",ComponentsThirdDoublet],['2"',ComponentsFirstDoublet]]
    
def ProductDoubleGamma333(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
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
    
def ProductDoubleGamma3(Multiplet1,Multiplet2):
    rep1, Components1 = Multiplet1
    rep2, Components2 = Multiplet2
    dim1 = len(Components1)
    dim2 = len(Components2)
    if dim1 != int(rep1[0]):
        print("The number of components of Multiplet1 does not match the dimension of its representation")
    if dim2 != int(rep2[0]):
        print("The number of components of Multiplet2 does not match the dimension of its representation")  
    if dim1 == 1 and dim2 == 1:
        return ProductDoubleGamma311(Multiplet1,Multiplet2)
    if dim1 == 1 and dim2 == 2:
        return ProductDoubleGamma312(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 1:
        return ProductDoubleGamma312(Multiplet2,Multiplet1)
    if dim1 == 1 and dim2 == 3:
        return ProductDoubleGamma313(Multiplet1, Multiplet2)
    if dim1 == 3 and dim2 == 1:
        return ProductDoubleGamma313(Multiplet2, Multiplet1)
    if dim1 == 2 and dim2 == 2:
        return ProductDoubleGamma322(Multiplet1,Multiplet2)
    if dim1 == 2 and dim2 == 3:
        return ProductDoubleGamma323(Multiplet1, Multiplet2)
    if dim1 == 3 and dim2 == 2:
        return ProductDoubleGamma323(Multiplet2,Multiplet1)
    if dim1 == 3 and dim2 == 3:
        return ProductDoubleGamma333(Multiplet1, Multiplet2)