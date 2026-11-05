# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================

import numpy as np
import cmath
import mpmath
import sympy
from fractions import Fraction

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from MinimalSet import *

q = sympy.symbols('q')

"This files builds the VVMFs of the minimal set. It says the user which Mathematica file to use and asks the coefficients of the q-expansions."
"The formulae used by Mathematica does not fix the overall constants of the entries of a VVMF. They are computed here."

#=========================
#From Coefficients to VVMF
#=========================

def Kernel(Matrix, tol=1e-10):
    #Gives a vector in the kernel of the matrix.
    "AI-generated"
    U, S, Vt = np.linalg.svd(Matrix, full_matrices=True)
    N = Vt.T[:, -1]
    return N


def ComputeConstants(Coeffs,L, IrrepS,weight):
    #At that point we know the q-expansions of the VVMF thanks to the hypergeometric series but not the overall constants.
    #To do so we will call the functions of the above section but we need a valuable function, not a formal expression.
    #This function does the conversion and then compute the constants.
    dim =  len(Coeffs)
    order = np.min([len(Coeffs[i]) for i in range(dim)])
    def VVMF(tau):
        x = np.exp(2*1j*np.pi*tau)
        f = [0 for _ in range(dim)]
        for j in range(order):
            for k in range(dim):
                f[k] += Coeffs[k][j]*x**j
        for j in range(dim):
            f[j] *= np.exp(2*1j*np.pi*tau*L[j])
        return[[f[i]] for i in range(dim)]
    S = np.matrix(IrrepS)
    real_tau = np.random.random()-0.5
    imaginary_tau = 3*np.random.random()+(1-real_tau**2)**0.5
    tau = real_tau + imaginary_tau*1j
    Y = np.diag([VVMF(tau)[k][0] for k in range(dim)])
    YS = np.diag([VVMF(-1/tau)[k][0] for k in range(dim)])
    A = YS-np.dot((-tau)**weight,np.dot(S,Y))
    NullEigenVector = Kernel(A)
    C = [complex(np.conj(NullEigenVector[k][0])) for k in range(len(NullEigenVector))]
    i = np.argmin([np.abs(C[k]) for k in range(len(C))])
    Constants  = [C[k]/C[i] for k in range(len(C))] 
    return NormalizeConstants(Constants)

def NormalizeConstants(Constants):
    #Put the constants calculated in a readable form if possible.
    for k in range(len(Constants)):
        if np.abs(Constants[k]-Constants[k].real) < 0.001:
            Constants[k] = float(Constants[k].real)
            if np.abs(round(Constants[k])-Constants[k]) < 0.001:
                Constants[k] = round(Constants[k])
            elif np.abs(round(Constants[k]**2)-Constants[k]**2) < 0.01:
                Signe = np.abs(Constants[k])/Constants[k]
                C = np.sqrt(round(Constants[k]**2))
                Constants[k] = Signe*C
            elif np.abs(round(1/Constants[k]**2)-1/Constants[k]**2) < 0.001:
                Signe = np.abs(Constants[k])/Constants[k]
                C = np.sqrt(round(1/Constants[k]**2))
                Constants[k] = Signe*C
    return Constants

def makeNormalizedf(Coeffs,L,Constants):
    #Construct the VVMF as a valuable function once the remaining constants are known.
    dim =  len(Coeffs)
    order = np.min([len(Coeffs[i]) for i in range(dim)])
    def Normalizedf(tau):
        x = np.exp(2*1j*np.pi*tau)
        f = [0 for _ in range(dim)]
        for j in range(order):
            for k in range(dim):
                f[k] += Coeffs[k][j]*x**j
        for j in range(dim):
            f[j] *= np.exp(2*1j*np.pi*tau*L[j])
        for j in range(dim):
            f[j] *= Constants[j]
        return[[f[i]] for i in range(dim)]
    return Normalizedf

def CheckVVMF(VVMF, IrrepT, IrrepS, k):
    #Check if the function VVMF transforms correcly under the given representation and weight.
    #To do this, it takes random values of the module in the fundamental domain and verifies the T and S transformations.
    T = np.matrix(IrrepT)
    S = np.matrix(IrrepS)
    for i in range(10000):
        real_tau = np.random.random()-0.5
        imaginary_tau = 3*np.random.random()+(1-real_tau**2)**0.5
        tau = real_tau + imaginary_tau*1j
        Y_tau = np.matrix(VVMF(tau))
        Y_Ttau = np.matrix(VVMF(tau+1))
        Y_Stau = np.matrix(VVMF(-1/tau))
        if not(MatricesEquality(Y_Ttau, np.dot(T,Y_tau))) or not(MatricesEquality(Y_Stau, (-tau)**k*np.dot(S,Y_tau))):
            return False
    return True    

#====================
#Minimal VVMF builder
#====================
def makeMinimalVVMF(Group,Rep,Coeffs,weight,ZeroVVMF=False):
    #Given a list of coefficients and the characteristics of the VVMF we want to build, it applies the functions above to
    #determine if the construction is possible and if yes builds the VVMF and shows the result.
    (i1,i2) = FindRepInGroupTheory(Rep, Group)
    T = Group[2][i1][i2]
    S = Group[1][i1][i2]
    L = makeL(T)
    dim = len(Coeffs)
    order = np.min([len(Coeffs[i]) for i in range(dim)])
    if ZeroVVMF:
        Constants = [1 for _ in range(len(Coeffs))]
        Bool = True
    else:
        Constants = ComputeConstants(Coeffs,L,S,weight)
        Normalizedf = makeNormalizedf(Coeffs,L, Constants,)
        Bool = CheckVVMF(Normalizedf,T,S,weight)
        if not(Bool):
            print("No overall constants have been found that permit to create a VVMF with the coefficients given.") 
        else:
            print("The overall constants are:")
            for k in range(dim):
                if np.abs(Constants[k]-Constants[k].real) < 0.001:
                    Constants[k] = float(Constants[k].real)
                    if np.abs(round(Constants[k])-Constants[k]) < 0.001:
                        print("C"+str(k+1)+" = "+str(round(Constants[k])))
                    elif np.abs(round(Constants[k]**2)-Constants[k]**2) < 0.01:
                        Signe = np.abs(Constants[k])/Constants[k]
                        C = sympy.sqrt(round(Constants[k]**2))
                        if Signe == 1.0:
                            expr = "C"+str(k+1)+" = "+str(C)
                        else:
                            expr = "C"+str(k+1)+" = -"+str(C)
                        print(expr)
                    elif np.abs(round(1/Constants[k]**2)-1/Constants[k]**2) < 0.001:
                        Signe = np.abs(Constants[k])/Constants[k]
                        C = sympy.sqrt(round(1/Constants[k]**2))
                        if Signe == 1.0:
                            expr = "C"+str(k+1)+" = 1/("+str(C)+")"
                        else:
                            expr = "C"+str(k+1)+" = -1/("+str(C)+")"
                        print(expr)
                else:
                    GlobalPhase = cmath.phase(Constants[k])/np.pi
                    if np.abs(GlobalPhase*12-round(GlobalPhase*12))<0.01:
                        GlobalPhase = round(GlobalPhase*12)/12
                        GlobalPhase = sympy.nsimplify(GlobalPhase)
                    modulus = np.abs(Constants[k])
                    if np.abs(round(modulus)/modulus-1)<0.01: 
                        modulus = round(modulus)
                    if np.abs(round(modulus**2)/modulus**2-1) < 0.01:
                        signe = modulus/np.abs(modulus)
                        value = signe*sympy.sqrt(round(modulus**2))
                    print("C"+str(k+1)+" = "+str(value*sympy.exp(sympy.I*GlobalPhase*sympy.pi)))
            print("It verifies the transforming properties of a VVMF for the representation "+Rep+" of the group "+Group[-1]+"\n")        
    somme = [0 for _ in range(dim)]
    for j in range(order):
        for k in range(dim):
            somme[k] += Coeffs[k][j]*q**j
    
    for j in range(dim):
        somme[j] *= q**L[j]
    for j in range(dim):
        somme[j] *= Constants[j]
    return (Bool, [Rep,[somme[i] for i in range(dim)]])


def ReshapeCoeffList(Coeffs):
    #Transforms the list of coeffs given by Mathematica in a usable format.
    if len(Coeffs) == 0 or Coeffs[0] != '{':
        return [0]
    ListCharacter = list(Coeffs)
    Dim = -1
    for elt in ListCharacter:
        if elt == '{':
            Dim += 1
    CoeffVVMF = []
    Coeff = ''
    for k in range(1,len(ListCharacter)-1):
        if ListCharacter[k] == '{':
            ListCoeff = []
        elif ListCharacter[k] == '}':
            CoeffVVMF.append(ListCoeff)
        elif ListCharacter[k] != ',':
            Coeff += ListCharacter[k]
        else:    
            ListCoeff.append(float(Coeff))
            Coeff = ''
    return CoeffVVMF

fileNames = ["Dim1VVMF_Expansion.nb","Dim2VVMF_Expansion.nb","Dim3VVMF_Expansion.nb","Dim4CyclicTypeVVMF_Expansion.nb","Dim4NonCyclicTypeVVMF_Expansion.nb"]



def ConstructMinimalVVMF(group_theory):
    #Asks the user to give the coefficients of the VVMFs in the minimal set copied from Mathematica and builds them.
    MinimalVVMFsProfile = DisplayMinimalVVMFProfile(group_theory)
    ListOfMinimalVVMF = []
    ListOfMinimalWeights = []
    order = 10000
    Irreps = group_theory[0]
    RepFirstApparition = [[] for _ in range(len(Irreps))]
    for k in range(len(Irreps)):
        RepFirstApparition[k] = [[True] for _ in range(len(Irreps[k]))]
    for k in range(len(MinimalVVMFsProfile)):
        rep = MinimalVVMFsProfile[k][0]
        (i1,i2) = FindRepInGroupTheory(rep, group_theory)
        dim = int(rep[0])
        Type = 0
        if dim == 4:
            Cyclicity = TypeDim4Rep(rep, group_theory)
        weight = MinimalVVMFsProfile[k][1]
        L = MinimalVVMFsProfile[k][2]
        if dim < 5:
            file = fileNames[int(rep[0])-1+Type]
            if RepFirstApparition[i1][i2]:
                print("Please to construct the minimal VVMF of weight "+str(weight)+" in rep "+rep+",")
                RepFirstApparition[i1][i2] = False
                Coeffs = input("enter L="+ToStringL(L)+" in the file "+file+" and copy the coefficients here:")
            else:
                if dim < 4 or (dim == 4 and Cyclicity):
                    print("This rep has already appeared in the minimal set. The associated VVMF is not the minimal one but a modular derivative.")
                    print("Please to construct the minimal VVMF of weight "+str(weight)+" in rep "+rep+",")
                    Coeffs = input("Copy the coefficients here:")
                else:
                    print("This rep has already appeared in the minimal set. The associated VVMF is not the minimal.")
                    print("Please to construct the minimal VVMF of weight "+str(weight)+" in rep "+rep+",")
                    Coeffs = input("Copy the coefficients here:")
        else:
            print("Please to construct the minimal VVMF of weight "+str(weight)+" in rep "+rep+",")
            Coeffs = input("copy the coefficients here:")
        print()
        CoeffVVMF = ReshapeCoeffList(Coeffs)
        if CoeffVVMF == [0]:
            CoeffVVMF = [[0 for _ in range(10)] for _ in range(dim)]
            (Bool, MinimalVVMF) = makeMinimalVVMF(group_theory, rep, CoeffVVMF, weight,ZeroVVMF = True)
        else:
            (Bool, MinimalVVMF) = makeMinimalVVMF(group_theory, rep, CoeffVVMF, weight)
            for k in range(len(CoeffVVMF)):
                if len(CoeffVVMF[k]) < order:
                    order = len(CoeffVVMF[k])
        if Bool:
            ListOfMinimalVVMF.append(MinimalVVMF)
            ListOfMinimalWeights.append(weight)
    print()    
    return ListOfMinimalVVMF, ListOfMinimalWeights, order
    



