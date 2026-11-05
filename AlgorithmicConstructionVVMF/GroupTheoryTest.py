# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================

import numpy as np

"""This file checks that the objects Group_Theory gathered in FullGroupTheory are well constructed. In particular it checks that:
1) The matrices S and T are squared matrices of size the dimension of the representation.
2) The number of irreducible representations is compatible with the order of the group.
3) S and T fulfill S^4=(ST)^3=T^N=Id.
4) S and T are expressed in the T-diagonal basis.
5) The character table is a square table with size the number of conjugacy classes.
6) The number of conjugacy classes is equal to the number of irreducible representations.
7) The character table is orthogonal for the canonical scalar product on characters."""

def ShapeTest(Irreps, SIrreps, TIrreps):
    if len(Irreps) < 3 or len(SIrreps) < 3 or len(TIrreps) < 3:
        return (False, "Missing representations")
    else:
        for i in range(len(Irreps)):
            if len(SIrreps[i]) != len(Irreps[i]) or len(TIrreps[i]) != len(Irreps[i]):
                return (False, "Numbers of dim-"+str(i+1)+" irreps do not match.")
        for k in range(len(Irreps)):
            for j in range(len(Irreps[k])):
                if len(SIrreps[k][j]) != k+1 or len(TIrreps[k][j]) != k+1:
                    return (False, "The "+str(j)+"-th dim-"+str(k+1)+" irrep is not a "+str(k+1)+"x"+str(k+1)+" matrix.")
                elif k > 0:
                    for l in range(k+1):
                        if len(SIrreps[k][j][l]) != k+1 or len(TIrreps[k][j][l]) != k+1:
                            return (False, "The "+str(j)+"-th dim-"+str(k+1)+" irrep is not a "+str(k+1)+"x"+str(k+1)+" matrix.")
        return (True, None)

def MatricesEquality(A, B):
    return np.max(np.abs(A-B)) <= 0.001

def GroupCardinal(GroupId):
    Cardinal = ""
    i = 1
    while GroupId[i] != ",":
        Cardinal += GroupId[i]
        i += 1
    Cardinal = int(Cardinal)
    return Cardinal

def NumberOfIrreps(GroupTheory):
    Irreps = GroupTheory[0]
    Sum = 0
    for k in range(len(Irreps)):
        Sum += (k+1)**2*len(Irreps[k])
    Cardinal = GroupCardinal(GroupTheory[-2])
    if Sum < Cardinal:
        return(False,"There are not enough irreps.")
    elif Sum > Cardinal:
        return(False, "There are too much irreps")
    else:
        return(True, None)

def ConstraintsTest(SIrreps, TIrreps, N):
    for k in range(len(SIrreps)):
        for i in range(len(SIrreps[k])):
            I = np.eye(k+1)
            S = np.matrix(SIrreps[k][i])
            T = np.matrix(TIrreps[k][i])
            ST = np.dot(S, T)
            if not(MatricesEquality(S**4, I)):
                return (False, "The "+str(i)+"-th dim "+str(k+1)+" irrep does not satisfy S^4=Id")
            if not(MatricesEquality(T**N, I)):
                return (False, "The "+str(i)+"-th dim "+str(k+1)+" irrep does not satisfy T^"+str(N)+"=Id")
            if not(MatricesEquality(ST**3, I)):
                return (False, "The "+str(i)+"-th dim "+str(k+1)+" irrep does not satisfy ST^3=Id")
    return (True, None)

def IsTDiagonal(TIrreps):
    for k in range(len(TIrreps)):
        for l in range(len(TIrreps[k])):
            for i in range(len(TIrreps[k][l])):
                if len(TIrreps[k][l]) > 1:
                    for j in range(len(TIrreps[k][l][i])):
                        if i != j and TIrreps[k][l][i][j] != 0:
                            return (False, "The "+str(l)+"-th dim "+str(k+1)+" is not in the T-diagonal basis.")
    return (True,None)

def CharacterTableShapeTest(Irreps, CharacterTable, ConjugacyClassesSizes):
    NumberOfConjugacyClasses = len(ConjugacyClassesSizes)
    if len(CharacterTable) < 3 :
        return (False, "Missing representations in the character table.")
    else:
        for i in range(len(Irreps)):
            if len(CharacterTable[i]) != len(Irreps[i]):
                return (False, "Numbers of dim-"+str(i+1)+" characters do not match.")
        for k in range(len(Irreps)):
            for j in range(len(Irreps[k])):
                if len(CharacterTable[k][j]) != NumberOfConjugacyClasses:
                    return (False, "The "+str(j)+"-th dim-"+str(k+1)+" character does not have the good number of entries.")
        return (True, None)
    
def NumberOfConjugacyClasses(ConjugacyClassesSizes, GroupId):
    Cardinal = GroupCardinal(GroupId)
    Sum = 0
    for elt in ConjugacyClassesSizes:
        Sum += elt
    if Sum != Cardinal:
        return(False, "The sum of the sizes of the conjugacy classes does not match the cardinal of the group.")
    return (True, None)

def ScalarProduct(GroupId,Character1, Character2, ConjugacyClassesSizes):
    Cardinal = GroupCardinal(GroupId)
    Sum = 0
    for k in range(len(Character1)):
        Sum += ConjugacyClassesSizes[k]*np.conj(Character1[k])*Character2[k]
    return Sum/Cardinal

def OrhtogonalityCharacterTableTest(GroupId, CharacterTable, ConjugacyClassesSizes):
    Characters = []
    for k in range(len(CharacterTable)):
        for l in range(len(CharacterTable[k])):
            Characters.append(CharacterTable[k][l])
    for i in range(len(Characters)):
        for j in range(i+1):
            scalar_product = ScalarProduct(GroupId, Characters[i], Characters[j], ConjugacyClassesSizes)
            if i == j:
                if np.abs(scalar_product-1) > 0.0001:
                    return(False, "The "+str(i+1)+"-th character is not normalized to 1.")
            else:
                if np.abs(scalar_product) > 0.0001:
                    return(False, "The "+str(j+1)+"-th and the "+str(i+1)+"-th characters are not orthogonal.")
    return (True,None)
    
def CharacterTheoryTest(Group_Theory):
    Irreps, SIrreps, TIrreps, N, ConjugacyClassesSizes, CharacterTable, GroupId, name = Group_Theory
    ShapeBool, Message = CharacterTableShapeTest(Irreps, CharacterTable, ConjugacyClassesSizes)
    if not(ShapeBool):
        return ShapeBool, Message
    else:
        NumberBool, Message = NumberOfConjugacyClasses(ConjugacyClassesSizes, GroupId)
        if not(NumberBool):
            return(NumberBool, Message)
        else:
            OrthogonalityBool, Message = OrhtogonalityCharacterTableTest(GroupId, CharacterTable, ConjugacyClassesSizes)
            return(OrthogonalityBool, Message)

def RepresentationTheoryTest(Group_Theory):
    Irreps, SIrreps, TIrreps, N, ConjugacyClassesSizes, CharacterTable, Id, name = Group_Theory
    ShapeBool, Message = ShapeTest(Irreps, SIrreps, TIrreps)
    if not(ShapeBool):
        return ShapeBool, Message
    else:
        ConstraintsBool, Message = ConstraintsTest(SIrreps, TIrreps, N)
        if not(ConstraintsBool):
            return(ConstraintsBool, Message)
        else:
            TDiagonalBool, Message = IsTDiagonal(TIrreps)
            return (TDiagonalBool,Message)

def FullGroupTheoryTest(FullGroupTheory):
    Error = False
    if len(FullGroupTheory) != 26:
        print("There are "+str(len(FullGroupTheory))+" groups when 20 are expected.")
        Error = True
    for group in FullGroupTheory:
        if len(group) != 8:
            print("Not the good number of objects in group_theory of "+group[-1])
            Error = True
        else:
            Bool, Message = RepresentationTheoryTest(group)
            if not(Bool):
                print(group[-1] + " is not well constructed. Error: " + Message)
                Error = True
            Bool, Message = NumberOfIrreps(group)
            if not(Bool):
                print(group[-1] + " is not well constructed. Error: " + Message)
                Error = True
            Bool, Message = CharacterTheoryTest(group)
            if not(Bool):
                 print(group[-1] + " is not well constructed. Error: " + Message)
                 Error = True
    if not(Error):
        print("Everything is okay")
