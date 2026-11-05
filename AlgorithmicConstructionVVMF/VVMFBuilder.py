# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================

import numpy as np
import sympy
import cmath
from fractions import Fraction

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from GroupTheoryDictionnary import FullGroupTheory
from MinimalVVMF import *

q = sympy.symbols('q')

"This file gives the q-expansions of any finite modular group."


group_theory = FullGroupTheory[1]      #Change this to specify which group you want to study.
from ClebschGordanLibrary.Gamma3_CG import ProductGamma3 as Product      #Change this to inform the code which Clesbch-Gordan coefficients must be used.
weight = 10       #The weight of the VVMFs you want the q-expansions.
FirstUse = True     #Ensure that it is True when you open the file or study a new group for the first time. Once the minimal set has been built, you can turn it to False.

#=====================
#Product of expansions
#=====================

def RoundExponents(expression, q, tolerance=0.001):
    #After a product, an expansion can contain integer power or float power. If a float power is an integer, it 
    #transforms it in an integer type so that the simplifcation is made. AI-generated.
    def Replace(expr):
        if expr.is_Pow and expr.base == q and expr.exp.is_Float:
            exp = float(expr.exp)
            n = round(exp)
            if  abs(exp - n) < tolerance:
                return (expr, q**n)
        return (expr, expr)
    replacements = expression.atoms(sympy.Pow)
    for p in replacements:
        motif, remplacement = Replace(p)
        expression = expression.replace(motif, remplacement)
    return expression

def SimplifyMultiplet(Multiplet,group_theory,n_max):
    #Simplify the product of q-expansion of every component of the multiplet.
    dim = len(Multiplet[1])
    rep = Multiplet[0]
    (i1,i2) = FindRepInGroupTheory(rep, group_theory)
    T = group_theory[2][i1][i2]
    L = makeL(T)
    for k in range(dim):
        Multiplet[1][k] = sympy.simplify(Multiplet[1][k]/q**(L[k])).series(q,n=n_max)
        for j in range(n_max): 
            Multiplet[1][k] = RoundExponents(Multiplet[1][k],q)
        Coeffs = getCoefficient(Multiplet[1][k])
        SimplifiedMultiplet = 0
        for l in range(len(Coeffs)):
            SimplifiedMultiplet += Coeffs[l]*q**l
        Multiplet[1][k] = q**L[k]*(SimplifiedMultiplet)
    return Multiplet

def SimplifyProduct(Product,group_theory,n_max):
    #The result of a product of multiplet can contain several multiplet. It simplifies each of them.
    n = len(Product)
    for k in range(n):
        Product[k] = SimplifyMultiplet(Product[k], group_theory,n_max)
    return Product

#==============
#ArrayStructure
#==============

def EmptyArrayVVMF(group_theory):
    #Create an empty array of the same size as the one collecting the irreps of the group.
    #In each entry of the array, we will collect every independent multiplet of VVMF of a given weight
    #transforming under their corresponding irreductible representation.
    Irreps = group_theory[0]
    ArrayVVMF = [[] for _ in range(len(Irreps))]
    for k in range(len(Irreps)):
        ArrayVVMF[k] = [[] for _ in range(len(Irreps[k]))]
    return ArrayVVMF


def SortInArray(Multiplet, ArrayVVMF, group_theory):
    #Put a multiplet of VVMF in the good entry of the array.
    rep = Multiplet[0]
    (i1,i2) = FindRepInGroupTheory(rep, group_theory)
    ArrayVVMF[i1][i2].append(Multiplet)
        
#==========================
#Annulation and redudancies
#==========================

def getCoefficient(expression):
    #Gives the coefficient of a q-expansion.
    poly = sympy.Poly(RoundExponents(expression.expand().removeO(),q),q)
    coeff_expansion = poly.all_coeffs()[::-1]
    coeff = []
    for k in range(len(coeff_expansion)):
        if type(coeff_expansion[k]) == sympy.core.numbers.ComplexInfinity:
            coeff.append(0)
        else:
            if np.isnan(complex(coeff_expansion[k])):
                coeff.append(0)
            elif np.abs(coeff_expansion[k]) < 0.01:
                coeff.append(0)
            else:
                coeff.append(complex(coeff_expansion[k]))
    return list(coeff)

def NullExpansion(expression):
    #Check if the expression can be considered to be zero.
    if expression == sympy.nan:
        return True
    if expression == 0:
        return True
    if expression.removeO() == 0:
        return True
    coeff = getCoefficient(expression)
    if np.max(np.abs(coeff)) <0.0001:
        return True
    return False

def NullMultiplet(Multiplet, group_theory):
    #Check if a multiplet can be considered to be zero. Can be if every component of the multiplet is zero.
    dim = len(Multiplet[1])
    rep = Multiplet[0]
    (i1,i2) = FindRepInGroupTheory(rep, group_theory)
    T = group_theory[2][i1][i2]
    L = makeL(T)
    for k in range(dim):
        if not(NullExpansion((Multiplet[1][k]/q**(L[k])).expand().removeO())):
            return False
    return True

def ExtractLeadingTerm(expression,group_theory): 
    #Gives the first non-vanishing term of a q-expansion.
    Coeffs = getCoefficient(expression.expand().removeO())
    i = 0
    while i < len(Coeffs) and Coeffs[i] == 0 :
        i += 1
    if i == len(Coeffs):
        return 1
    GlobalPhase = cmath.phase(Coeffs[i]) 
    for k in range(len(Coeffs)):
        Coeffs[k] = (Coeffs[k]/np.exp(1j*GlobalPhase)).real
    if len(Coeffs) <= 1:
        return 1
    RoundCoeffs = []
    for j in range(len(Coeffs)):
        RoundCoeffs.append(np.abs(round(Coeffs[j])))
    max_sqrt = 1
    for k in sympy.divisors(round(Coeffs[i]**2)):
         Divide = True
         for j in range(len(Coeffs)):
             if np.abs(Coeffs[j]) != 0 and np.abs(Coeffs[j]) < 10**5 and (np.abs((round(np.abs(Coeffs[j])/k**0.5))/(np.abs(Coeffs[j])/k**0.5)-1)>0.001):
                 Divide = False
         if Divide:
             max_sqrt = k
    return np.exp(1j*GlobalPhase)*max_sqrt**0.5 
   
def MultipletToCoefficient(Multiplet, group_theory):
    #For a multiplet, gives the list of the coefficients of the q-expansion of each of its component once the leading term
    #has been factorised out.
    dim = len(Multiplet[1])
    rep = Multiplet[0]
    (i1,i2) = FindRepInGroupTheory(rep, group_theory)
    T = group_theory[2][i1][i2]
    L = makeL(T)
    MultipletCoeffs = []
    for k in range(dim):
        Expression = (Multiplet[1][k]/q**(L[k])).expand().removeO()
        LeadingTerm = ExtractLeadingTerm(Expression, group_theory)
        Coeffs = getCoefficient(Expression/LeadingTerm)
        for i in range(len(Coeffs)):
            Coeffs[i] = round(Coeffs[i].real)
        MultipletCoeffs.append(Coeffs)
    return MultipletCoeffs

def LinearlyIndependenceTest(ListOfMultiplet,group_theory):
    #Given a list of multiplets, this tests if they are linearly dependent or not.
    "We use the fact that a multiplet of modular forms are linearly independent iff their first entry are linearly independent"
    if len(ListOfMultiplet) <= 1: 
        return True
    MatrixOfCoefficients = []
    for k in range(len(ListOfMultiplet)):
        MatrixOfCoefficients.append(MultipletToCoefficient(ListOfMultiplet[k], group_theory)[0])
    MinSize = np.min([len(MatrixOfCoefficients[k]) for k in range(len(MatrixOfCoefficients))])
    for k in range(len(MatrixOfCoefficients)):
        MatrixOfCoefficients[k] = MatrixOfCoefficients[k][:MinSize]
    MatrixOfCoefficients = np.matrix(MatrixOfCoefficients)
    m = np.linalg.matrix_rank(MatrixOfCoefficients,tol=1e-6)
    n = len(ListOfMultiplet)
    return n == m

def Kernel(Matrix, tol=1e-6):
    #Gives a vector in the kernel of the matrix.
    "AI-generated"
    U, S, Vt = np.linalg.svd(Matrix, full_matrices=True)
    N = Vt.T[:, -1]
    return N

def EraseRedundanciesOfSameRep(ListOfMultiplet,group_theory):
    #Given a list of multiplet of the same rep, if they are linearly dependent, it erases some multiplets until the list becomes independent.
    if LinearlyIndependenceTest(ListOfMultiplet, group_theory):
        return ListOfMultiplet
    MatrixOfCoefficients = []
    for k in range(len(ListOfMultiplet)):
        MatrixOfCoefficients.append(MultipletToCoefficient(ListOfMultiplet[k],group_theory)[0])
    MinSize = np.min([len(MatrixOfCoefficients[k]) for k in range(len(MatrixOfCoefficients))])
    for k in range(len(MatrixOfCoefficients)):
        MatrixOfCoefficients[k] = MatrixOfCoefficients[k][:MinSize]
    MatrixOfCoefficients = np.matrix(MatrixOfCoefficients)
    MatrixOfCoefficients = np.transpose(MatrixOfCoefficients)
    EigVector = Kernel(MatrixOfCoefficients)
    IndexRedundantVVMF = 0
    while np.abs(EigVector[IndexRedundantVVMF]) < 0.001:
        IndexRedundantVVMF += 1
    ListOfMultiplet.pop(IndexRedundantVVMF)
    return EraseRedundanciesOfSameRep(ListOfMultiplet, group_theory)
    
def EraseRedundanciesAmongVVMFs(ListOfMultiplet, group_theory):
    #Given a list of multiplet, they are sorted by multiplet of the same rep before being erased
    NextWeightArray = EmptyArrayVVMF(group_theory)
    for i in range(len(ListOfMultiplet)):
        if not(NullMultiplet(ListOfMultiplet[i], group_theory)):
            SortInArray(ListOfMultiplet[i],NextWeightArray, group_theory)
    ListOfVVMF = []
    for k in range(len(NextWeightArray)):
        for j in range(len(NextWeightArray[k])):
            ListOfVVMF += EraseRedundanciesOfSameRep(NextWeightArray[k][j], group_theory)
    return ListOfVVMF
        
#==========
#Formatting
#==========  

def CanonicalForm(value,LastTour): 
    #Check if value is close to an integer or a square root of an integer. If complex try to put the phase in the form pi*fraction.
    value = complex(value)
    if np.abs(value)<0.1:
        return 0
    GlobalPhase = cmath.phase(value)/np.pi
    if np.abs(GlobalPhase) < 0.01:
        GlobalPhase = 0
    elif np.abs(GlobalPhase*12-round(GlobalPhase*12))<0.01:
        GlobalPhase = round(GlobalPhase*12)/12
    AbsoluteValue = np.abs(value)
    if np.abs(round(AbsoluteValue)-AbsoluteValue)<0.1: 
        if LastTour:
            GlobalPhase = sympy.nsimplify(GlobalPhase)
            return sympy.exp(GlobalPhase*sympy.I*sympy.pi)*(round(AbsoluteValue))
        return np.exp(GlobalPhase*1j*np.pi)*round(AbsoluteValue)
    if np.abs(round(AbsoluteValue**2)-AbsoluteValue**2) < 0.05 and LastTour:
        GlobalPhase = sympy.nsimplify(GlobalPhase)
        return sympy.exp(GlobalPhase*sympy.I*sympy.pi)*sympy.sqrt(round(AbsoluteValue**2))
    else:
        return value


def NormalizeMultiplet(Multiplet,LastTour,group_theory):
    #A multiplet contains different q-expansions. We want them to be all made of integer coefficients and that at least
    #one starts with a 1.
    rep = Multiplet[0]
    (i1,i2) = FindRepInGroupTheory(rep, group_theory)
    T = group_theory[2][i1][i2]
    L = makeL(T)
    FirstNonZeroTerm = []
    for k in range(len(Multiplet[1])):
        FirstNonZeroTerm.append(ExtractLeadingTerm((Multiplet[1][k]/q**L[k]).expand(), group_theory))
    OverallConstantIndex = np.argmin([np.abs(FirstNonZeroTerm[i]) for i in range(len(FirstNonZeroTerm))])
    OverallConstant = FirstNonZeroTerm[OverallConstantIndex]
    Coeffs = MultipletToCoefficient(Multiplet,group_theory)
    for k in range(len(Multiplet[1])): 
        Multiplet[1][k] = q**L[k]*sympy.simplify(Multiplet[1][k]/q**L[k]/FirstNonZeroTerm[k])
        ConstantExpression = CanonicalForm(FirstNonZeroTerm[k]/OverallConstant,LastTour)
        CanonicalMultiplet = 0
        for i in range(len(Coeffs[k])):
            CanonicalMultiplet = Coeffs[k][i]*q**i + CanonicalMultiplet
        Multiplet[1][k] = q**Fraction(L[k]).limit_denominator()*CanonicalMultiplet
        if ConstantExpression != 1:
            Multiplet[1][k] = ConstantExpression*Multiplet[1][k]
    return Multiplet    

#===============================
#Construction by tensor products
#===============================

def Combinatory(weight, MinimalWeights):
    #This function gives every possible tensor product powers to reach a given weight starting from VVMFs of weights collected 
    #in MinimalWeights.
    if len(MinimalWeights) == 1:
        if weight % MinimalWeights[-1] == 0:
            return [[weight//MinimalWeights[-1]]]
        else:
            return []
    PossiblePowers = []
    for a in range(weight//MinimalWeights[-1]+1):
        res = Combinatory(weight-MinimalWeights[-1]*a, MinimalWeights[:-1])
        for elt in res:
            PossiblePowers.append(elt+[a])
    return PossiblePowers

def TensorProduct(ListOfPowers, ListOfMinimalVVMF, group_theory, n_max):
    #Every possible powers are known, this function does all these possible products.
    ListOfVVMF = [[group_theory[0][0][0],[1]]]
    for k in range(len(ListOfPowers)):
        for j in range(ListOfPowers[k]):
            NewListOfVVMF = []
            for i in range(len(ListOfVVMF)):
                NewVVMF = Product(ListOfVVMF[i],ListOfMinimalVVMF[k])
                NewVVMF = SimplifyProduct(NewVVMF, group_theory, n_max)
                NewListOfVVMF += NewVVMF
                NewListOfVVMF = EraseRedundanciesAmongVVMFs(NewListOfVVMF, group_theory)
            for l in range(len(NewListOfVVMF)):
                NewListOfVVMF[l] = NormalizeMultiplet(NewListOfVVMF[l], False, group_theory)
            ListOfVVMF = NewListOfVVMF
    return ListOfVVMF
        
def CheckNumberOfVVMF(ArrayOfVVMF, weight, group_theory):
    #The result of the previous step must agree the number of independent VVMF of a given weight and rep that we can
    #compute explicitly (see VVMF.py).
    Problem = False
    for k in range(len(ArrayOfVVMF)):
        for j in range(len(ArrayOfVVMF[k])):
            if len(ArrayOfVVMF[k][j]) != DimMRep(weight,group_theory[0][k][j],group_theory):
                Problem = True
                print("There are "+str(len(ArrayOfVVMF[k][j]))+" VVMF of rep "+group_theory[0][k][j]+" of weight "+str(weight)+" when "+str(DimMRep(weight,group_theory[0][k][j],group_theory))+" are expected")
                if len(ArrayOfVVMF[k][j]) < DimMRep(weight,group_theory[0][k][j],group_theory):
                    print("This is probably due to a too small set of minimal VVMFs. Maybe some minimal VVMFs are in high dimensional irrep. An other possibility is that you are in a case where the lowest weight VVMFs are not enough to build every VVMFs. ")
                else:
                    print("Something went wrong when erasing redundant VVMFs. Please verify the linear independence of the result.")
                print()
    return Problem

def DisplayVVMF(ArrayOfVVMF, weight, group_theory,Problem):
    #Write an array of VVMF in a readable manner.
    number = 0
    for k in range(len(ArrayOfVVMF)):
        for j in range(len(ArrayOfVVMF[k])):
            for i in range(len(ArrayOfVVMF[k][j])):
                ArrayOfVVMF[k][j][i] = NormalizeMultiplet(ArrayOfVVMF[k][j][i], True, group_theory)
                number += len(ArrayOfVVMF[k][j][i][1])
    if Problem:
        print(str(number)+" VVMFs of the group "+group_theory[-1]+" of weight "+str(weight)+" have been found.")
    else:
        print("The group " + group_theory[-1]+" has "+str(number)+" VVMFs of weight "+str(weight))
    print("They are gathered in the following multiplets:\n")
    for k in range(len(ArrayOfVVMF)):
        for j in range(len(ArrayOfVVMF[k])):
            n = len(ArrayOfVVMF[k][j])
            if n > 0:
                if n==1:
                    print("--> "+str(n)+" multiplet transforming in the representation "+group_theory[0][k][j]+".")
                    print("Its q-expansion is:")
                    print()
                else:
                    print("--> "+str(n)+" multiplets transforming in the representation "+group_theory[0][k][j]+".")
                    print("Their q-expansion are:")
                    print()
            for i in range(n):
                if i > 0:
                    print()
                    print("And")
                    print()
                for l in range(len(ArrayOfVVMF[k][j][i][1])):
                    print(ArrayOfVVMF[k][j][i][1][l])
                    print()

def BuildVVMF(weight,group_theory):
    #Put everything together: first determine every tensor product possible, apply them, check if we get the good number of VVMFs
    # and display the result.
    global ListOfMinimalVVMF
    global ListOfMinimalWeights
    global FirstUse
    global order
    
    if FirstUse: 
        ListOfMinimalVVMF, ListOfMinimalWeights, order = ConstructMinimalVVMF(group_theory)
    if len(ListOfMinimalVVMF) == 0:
        print("No minimal VVMF has been built.")
        return
    else:
        ListOfPossiblePowers = Combinatory(weight, ListOfMinimalWeights)
        if len(ListOfPossiblePowers) == 0:
                print("There is no VVMF of this weight.")
        VVMFs = []
        for k in range(len(ListOfPossiblePowers)):
            VVMFs += TensorProduct(ListOfPossiblePowers[k],ListOfMinimalVVMF, group_theory, order)
            VVMFs = EraseRedundanciesAmongVVMFs(VVMFs, group_theory)
        ArrayVVMFs = EmptyArrayVVMF(group_theory)
        for i in range(len(VVMFs)):
            SortInArray(VVMFs[i], ArrayVVMFs, group_theory)
        Problem = CheckNumberOfVVMF(ArrayVVMFs, weight, group_theory)
        DisplayVVMF(ArrayVVMFs,weight, group_theory,Problem)
    return VVMFs


VVMF = BuildVVMF(weight,group_theory)

