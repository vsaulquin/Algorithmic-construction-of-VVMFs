# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================

import numpy as np
import cmath
from fractions import Fraction

"This file aims to determine the minimal set of VVMFs by appling the procedure based on faithfulness of reps in the set."

#=====================
#Some useful functions
#=====================

def MatricesEquality(A,B):
    return np.max(np.abs(A-B)) <= 0.1

def makeL(IrrepT):
    #Given a matrix T, it constructs the L matrix associated i.e. the list of real numbers ri between 0 and 1 such that
    #the i-th eigenvalue of T can be written exp(2i*pi*ri).
    dim = len(IrrepT)
    T = np.matrix(IrrepT)
    EigvalT = np.linalg.eigvals(T)
    L = []
    for k in range(dim):
        r = cmath.phase((EigvalT[k]))%(2*np.pi)/(2*np.pi)
        if np.abs(r-1)<0.000001:
            r = 0
        L.append(r)
    return L


def FindRepInGroupTheory(rep, group):
    #Find the position of a rep in the list of irreps.
    Irreps, SIrreps, TIrreps, N, ConjugacyClassesSizes, CharacterTable, Id, Name = group
    for i in range(len(Irreps)):
        for k in range(len(Irreps[i])):
            if Irreps[i][k] == rep:
                return (i,k)
            

def AreOddWeightPossible(group):
    #Check if the group has irrep such that rep(S^2) = -Id. If so odd weight VVMFs exist.
    Irreps = group[0]
    SIrreps = group[1]
    for k in range(len(Irreps)):
        for j in range(len(Irreps[k])):
            S = np.matrix(SIrreps[k][j])
            I = np.eye(k+1)
            if MatricesEquality(S**2, (-1)*I):
                return True
    return False


def TypeDim4Rep(rep,group):
    #A dim 4 rep can be of two types: cyclic or non-cyclic.
    i1,i2 = FindRepInGroupTheory(rep,group)
    S = np.matrix(group[1][i1][i2])
    T = np.matrix(group[2][i1][i2])
    L = makeL(T)
    Trace = 0
    for k in range(len(L)):
        Trace += L[k]
    exponent = 3*Trace%2+1
    Sign = (-1)**exponent
    S = np.matrix(group[1][i1][i2])
    I = np.eye(4)
    CyclicType = MatricesEquality(np.dot(S,S), np.dot(Sign,I)) 
    return CyclicType

#===============
#Number of VVMFs
#===============

def MinimalWeight(rep, group):
    #Given a representation, it computes the weight of the first existing VVMF transforming under this representation.
    (i1,i2) = FindRepInGroupTheory(rep, group)
    IrrepT = group[2][i1][i2]
    dim = len(IrrepT)
    L_T = makeL(IrrepT)
    Tr_L = np.sum(L_T)
    if dim == 1:
        return round(12*Tr_L)
    elif dim == 2:
        return round(6*Tr_L-1)
    elif dim == 3:
        return round(4*Tr_L-2)
    elif dim == 4:
        Type = TypeDim4Rep(rep, group)
        if Type:
            return round(3*Tr_L-3)
        else:
            return round(3*Tr_L-2)
    elif dim == 5:    
        i = 12*Tr_L%5
        if i == 0:
            return round((12*Tr_L-20)/5)
        if i == 1:
            return round((12*Tr_L-16)/5)
        if i == 2:
            return round((12*Tr_L-12)/5)
        if i == 3:
            return round((12*Tr_L-8)/5)
        if i == 4:
            return round((12*Tr_L-14)/5)

def DimMRep(k, rep, group):
    #For the representation rep, it gives the dimensions of the vector space Mk(rep) containing every VVMF of weight k 
    # transforming under k.
    (dim_index, index) = FindRepInGroupTheory(rep, group)
    dim = dim_index+1
    S = np.matrix(group[1][dim_index][index])
    T = np.matrix(group[2][dim_index][index])
    ST = np.dot(S,T)
    if not(MatricesEquality(S**2,(-1)**k*np.eye(dim))):
           return 0
    elif k > 1:
        L = makeL(group[2][dim_index][index])
        omega = np.exp(2*np.pi*1j/3)
        a1 = np.trace(S**3)
        a2 = np.trace(np.dot(S**3,T))
        a3 = np.trace(ST**2)
        a4 = np.sum(L)
        return (round(np.abs((5+k)*dim/12+1j**k*a1/4+(-omega**2)**k*a2/(3*(1-omega))+omega**k*a3/(3*(1-omega**2))-a4)))
    else:
        k0 = MinimalWeight(rep, group)
        if k0 == 1:
            return 1
        else:
            return 0
    
def DimMGroup(k, group):
    #For a group, it gives the number of independent VVMF of weight k.
    Irreps = group[0]
    dim = 0
    for j in range(len(Irreps)):
        for i in range(len(Irreps[j])):
            dim += (j+1)*DimMRep(k,Irreps[j][i],group)
    return dim

def MinimalWeightGroup(group):
    #Gives the first weight such that there is a non-zero VVMF under this group.
    k = 1
    while DimMGroup(k,group) == 0:
        k += 1
    return k

#==================== 
#Minimal set of VVMFs
#====================

def FaithfullnessRepresentation(rep,group):
    #Test if a rep is faithful.
    (i1,i2)= FindRepInGroupTheory(rep, group)
    character = group[5][i1][i2]
    d = character[0]
    for k in range(1,len(character)):
        if character[k] == d:
            return False
    return True

def FaithfullnessSetOfRepresentations(reps, group):
    #Test if a set of reps is faithful in the sense that the sum is faithful.
    characters = []
    for rep in reps:
        (i1,i2) = FindRepInGroupTheory(rep, group)
        characters.append(group[5][i1][i2])
    character = []
    for k in range(len(characters[0])):
        Sum = 0
        for i in range(len(characters)):
            Sum += characters[i][k]
        character.append(Sum)
    d = character[0]
    for j in range(1,len(character)):
        if character[j] == d:
            return False
    return True


def getMinimalSet(group):
    #Determine the minimal set of VVMFs.
    OddWeightBool = AreOddWeightPossible(group)
    PossibleHighDimProblemBool = False
    if OddWeightBool:
        if len(group[0]) > 5:
            PossibleHighDimProblemBool = True
    k0 = MinimalWeightGroup(group)
    MinimalReps = []
    MinimalSet = []
    Irreps = []
    for i in range(len(group[0])):
        Irreps += group[0][i]
    for rep in Irreps:
        if DimMRep(k0, rep, group) > 0:
            MinimalReps.append(rep)
            i1, i2 = FindRepInGroupTheory(rep, group)
            L = makeL(group[2][i1][i2])
            MinimalSet.append((rep,k0,L))
    k = k0
    while not(FaithfullnessSetOfRepresentations(MinimalReps, group)):
        k += 1
        for rep in Irreps:
            if DimMRep(k, rep, group) > 0:
                MinimalReps.append(rep)
                i1, i2 = FindRepInGroupTheory(rep, group)
                L = makeL(group[2][i1][i2])
                MinimalSet.append((rep,k,L))
    return (OddWeightBool, PossibleHighDimProblemBool, MinimalSet)


def ToStringL(L):
    #Transforms L in a practical format.
    LString = "["
    for elt in L:
        LString += str(Fraction(elt).limit_denominator())
        LString += ","
    LString = LString[:-1]
    LString += "]"
    return LString
 

def DisplayMinimalVVMFProfile(group):
    #Display the result of the previous function.
    OddWeightBool, PossibleHighDimProblemBool, MinimalSet = getMinimalSet(group)
    MinimalWeights = []
    MinimalReps = []
    for k in range(len(MinimalSet)):
        MinimalWeights.append(MinimalSet[k][1])
        MinimalReps.append(MinimalSet[k][0]) 
    if OddWeightBool:
        print("This group has odd-weight VVMFs.")
        print()
    else:
        print("This group only has even-weight VVMFs.")
        print()
    if PossibleHighDimProblemBool:
            print("NOTE: This group has irreps of dimension greater than 6. If there exist weight 1 VVMFs transforming in one of this irrep, the code cannot know it since no formula exists.")
            print("In this case there will probably be missing VVMFs.")
            print()
    DifferentWeights = [MinimalWeights[0]]
    for k in range(1,len(MinimalWeights)):
        if MinimalWeights[k] != DifferentWeights[-1]:
            DifferentWeights.append(MinimalWeights[k])
    if len(DifferentWeights) > 1:
        print("Every VVMF of the group "+group[-1]+" is generated by "+str(len(DifferentWeights))+" sets of minimal VVMfs.")
    else:
        print("Every VVMF of the group "+group[-1]+" is generated by "+str(len(DifferentWeights))+" set of minimal VVMFs.")
    print()
    for k in range(len(DifferentWeights)):
        if k == 0:
            print("-->VVMFs of the first set have weight "+str(DifferentWeights[k])+".")
        elif k == 1:
            print("-->VVMFs of the second set have weight "+str(DifferentWeights[k])+".")
        elif k== 2:
            print("-->VVMFs of the third set have weight "+str(DifferentWeights[k])+".")
        else:
            print("VVMFs of the "+str(k+1)+"-th set have weight "+str(DifferentWeights[k])+".")
        print("In this set there are "+str(DimMGroup(DifferentWeights[k],group))+" independent minimal VVMFs.")
        print("They transform in:")
        for j in range(len(MinimalWeights)):
            if  MinimalWeights[j] == DifferentWeights[k]:
                rep = MinimalReps[j]
                (dim_index,index) = FindRepInGroupTheory(rep, group)
                L = makeL(group[2][dim_index][index])
                expr = "• One independent "+rep+" representation. This rep has L="+ToStringL(L)+"."
                if dim_index == 3:
                    CyclicType = TypeDim4Rep(rep, group)
                    if CyclicType:
                        expr += " It is a 4-dim rep of cyclic type."
                    else:
                        expr += " It is a 4-dim rep of non-cylic type."
                print(expr)
        print()
    return MinimalSet

