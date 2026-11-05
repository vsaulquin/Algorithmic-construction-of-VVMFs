# -*- coding: utf-8 -*-
#==================
#@author: vsaulquin
#==================

import numpy as np

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from GroupTheoryTest import FullGroupTheoryTest

"""This file is a dictionnary of representations and characters for different quotient groups of SL(2,Z). 
For a given group, what we collect is:
--> The list of all irreducible representations  of the group.
--> For each of this representation the matrix representation of S.
--> For each of this representation the matrix representation of T. Express S and T in the T-diagonal basis.
--> The level of the group i.e. the smallest non-zero integer N such that T^N=1.
--> A list containing the sizes of the conjugacy classes of the group.
--> The character tables of the group.
--> The GAP Id of the group.
-->The name of the group.
All these data are collected in an object called Group_Theory. There is one of this object for each group of interest.
The list FullGroupTheory collects every Group_Theory in a single object. At the end, the file GroupTheoryTest is called to check that
the objects Group_Theory are well constructed."""

omega = np.exp(2*1j*np.pi/3)
zeta = np.exp(2*1j*np.pi/5)
phi = (1+5**0.5)/2
xi = np.exp(2*1j*np.pi/8)
eta = np.exp(2*1j*np.pi/12)
theta = np.exp(2*1j*np.pi/7)
delta = np.exp(2*1j*np.pi/13)

def MultScal(scalaire, matrice):
    result = [[element * scalaire for element in ligne] for ligne in matrice]
    return result


#=================================================
#Class of Gamma_N groups and their double covering
#=================================================

"Gamma_2 = S3. This group has no double covering"
IrrepsGamma2 = [[str(1), str(1)+"'"],[str(2)],[]]

TIrreps1Gamma2 = [[1],[-1]]
SIrreps1Gamma2 = [[1],[-1]]
TIrreps2Gamma2 = [[[1,0],[0,-1]]]
SIrreps2Gamma2 = [[[-1/2,-3**0.5/2],[-3**0.5/2,1/2]]] 
TIrreps3Gamma2 = []
SIrreps3Gamma2 = []

TIrrepsGamma2 = [TIrreps1Gamma2, TIrreps2Gamma2, TIrreps3Gamma2]
SIrrepsGamma2 = [SIrreps1Gamma2, SIrreps2Gamma2, SIrreps3Gamma2]

ConjugacyClassesSizesGamma2 = [1,3,2]

CharacterTable1Gamma2 = [[1,1,1],[1,-1,1]]
CharacterTable2Gamma2 = [[2,0,-1]]
CharacterTable3Gamma2 = []

CharacterTableGamma2 = [CharacterTable1Gamma2, CharacterTable2Gamma2, CharacterTable3Gamma2]

Gamma2_Theory = [IrrepsGamma2, SIrrepsGamma2, TIrrepsGamma2, 2, ConjugacyClassesSizesGamma2, CharacterTableGamma2, "[6,1]", "Gamma_2"]


"Gamma_3 = A4"
IrrepsGamma3 = [[str(1), str(1)+"'", str(1)+'"'],[],["3"]]

TIrreps1Gamma3 = [[1],[omega],[omega**2]]
SIrreps1Gamma3 = [[1],[1],[1]]
TIrreps2Gamma3 = []
SIrreps2Gamma3 = [] 
TIrreps3Gamma3 = [[[1,0,0],[0,omega,0],[0,0,omega**2]]]
SIrreps3Gamma3 = [MultScal(1/3,[[-1,2,2],[2,-1,2],[2,2,-1]])]

TIrrepsGamma3 = [TIrreps1Gamma3, TIrreps2Gamma3, TIrreps3Gamma3]
SIrrepsGamma3 = [SIrreps1Gamma3, SIrreps2Gamma3, SIrreps3Gamma3]

ConjugacyClassesSizesGamma3 = [1,4,4,3]

CharacterTable1Gamma3 = [[1,1,1,1],[1,omega,omega**2,1],[1,omega**2,omega,1]]
CharacterTable2Gamma3 = []
CharacterTable3Gamma3 = [[3,0,0,-1]]

CharacterTableGamma3 = [CharacterTable1Gamma3, CharacterTable2Gamma3, CharacterTable3Gamma3]

Gamma3_Theory = [IrrepsGamma3, SIrrepsGamma3, TIrrepsGamma3, 3, ConjugacyClassesSizesGamma3, CharacterTableGamma3, "[12,3]", "Gamma_3"]

"Double covering of Gamma_3: T' called DoubleGamma3 here"
IrrepsDoubleGamma3 = [[str(1), str(1)+"'", str(1)+'"'],[str(2),str(2)+"'", str(2)+'"'],["3"]]

A1 = [[omega,0],[0,1]]
A2 = [[-1j/3**0.5, -(2/3)**0.5*np.exp(1j*np.pi/12)],[(2/3)**0.5*np.exp(-1j*np.pi/12),1j/3**0.5]]

TIrreps1DoubleGamma3 = [[1],[omega],[omega**2]]
SIrreps1DoubleGamma3 = [[1],[1],[1]]
TIrreps2DoubleGamma3 = [A1, MultScal(omega, A1), MultScal(omega**2, A1)]
SIrreps2DoubleGamma3 = [A2,A2,A2]
TIrreps3DoubleGamma3 = [[[1,0,0],[0,omega,0],[0,0,omega**2]]]
SIrreps3DoubleGamma3 = [MultScal(1/3,[[-1,2,2],[2,-1,2],[2,2,-1]])]

TIrrepsDoubleGamma3 = [TIrreps1DoubleGamma3, TIrreps2DoubleGamma3, TIrreps3DoubleGamma3]
SIrrepsDoubleGamma3 = [SIrreps1DoubleGamma3, SIrreps2DoubleGamma3, SIrreps3DoubleGamma3]

ConjugacyClassesSizesDoubleGamma3 = [1,1,6,4,4,4,4]

CharacterTable1DoubleGamma3 = [[1,1,1,1,1,1,1],[1,1,1,omega,omega**2,omega,omega**2],[1,1,1,omega**2,omega,omega**2,omega]]
CharacterTable2DoubleGamma3 = [[2,-2,0,omega**2,-omega,-omega**2,omega],[2,-2,0,1,-1,-1,1],[2,-2,0,omega,-omega**2,-omega,omega**2]]
CharacterTable3DoubleGamma3 = [[3,3,-1,0,0,0,0]]

CharacterTableDoubleGamma3 = [CharacterTable1DoubleGamma3, CharacterTable2DoubleGamma3, CharacterTable3DoubleGamma3]

DoubleGamma3_Theory = [IrrepsDoubleGamma3, SIrrepsDoubleGamma3, TIrrepsDoubleGamma3, 3, ConjugacyClassesSizesDoubleGamma3, CharacterTableDoubleGamma3, "[24,3]", "Gamma_3'"]

"Gamma_4 = S4"
IrrepsGamma4 = [[str(1),str(1)+"'"],[str(2)],[str(3),str(3)+"'"]]

A1 = [[1,0],[0,-1]]
A2 = [[-1/2, 3**0.5/2],[3**0.5/2, 1/2]]
B1 = [[-1,0,0],[0,-1j,0],[0,0,1j]]
B2 = [[0,-1/2**0.5,-1/2**0.5],[-1/2**0.5,1/2,-1/2],[-1/2**0.5,-1/2,1/2]]

TIrreps1Gamma4 = [[1],[-1]]
SIrreps1Gamma4 = [[1],[-1]]
TIrreps2Gamma4 = [A1]
SIrreps2Gamma4 = [A2] 
TIrreps3Gamma4 = [B1,MultScal(-1,B1)]
SIrreps3Gamma4 = [B2,MultScal(-1,B2)]

TIrrepsGamma4 = [TIrreps1Gamma4, TIrreps2Gamma4, TIrreps3Gamma4]
SIrrepsGamma4 = [SIrreps1Gamma4, SIrreps2Gamma4, SIrreps3Gamma4]

ConjugacyClassesSizesGamma4 = [1,6,8,3,6]

CharacterTable1Gamma4 = [[1,1,1,1,1],[1,-1,1,1,-1]]
CharacterTable2Gamma4 = [[2,0,-1,2,0]]
CharacterTable3Gamma4 = [[3,1,0,-1,-1],[3,-1,0,-1,1]]

CharacterTableGamma4 = [CharacterTable1Gamma4, CharacterTable2Gamma4, CharacterTable3Gamma4]

Gamma4_Theory = [IrrepsGamma4, SIrrepsGamma4, TIrrepsGamma4, 4, ConjugacyClassesSizesGamma4, CharacterTableGamma4, "[24,12]", "Gamma_4"]

"Double covering of Gamma_4: S4'=SL(2,Z4) called DoubleGamma4 here."
IrrepsDoubleGamma4 = [[str(1),str(1)+"'",str(1)+"^", str(1)+"^'"],[str(2),str(2)+"^"],[str(3),str(3)+"'",str(3)+"^", str(3)+"^'"]]

TIrreps1DoubleGamma4 = [[1],[-1],[-1j],[1j]]
SIrreps1DoubleGamma4 = [[1],[-1],[1j],[-1j]]
TIrreps2DoubleGamma4 = [A1,MultScal(-1j,A1)]
SIrreps2DoubleGamma4 = [A2,MultScal(1j,A2)] 
TIrreps3DoubleGamma4 = [B1,MultScal(-1,B1),MultScal(-1j,B1),MultScal(1j,B1)]
SIrreps3DoubleGamma4 = [B2,MultScal(-1,B2),MultScal(1j,B2),MultScal(-1j,B2)]

TIrrepsDoubleGamma4 = [TIrreps1DoubleGamma4, TIrreps2DoubleGamma4, TIrreps3DoubleGamma4]
SIrrepsDoubleGamma4 = [SIrreps1DoubleGamma4, SIrreps2DoubleGamma4, SIrreps3DoubleGamma4]

ConjugacyClassesSizesDoubleGamma4 = [1,6,1,8,3,6,6,8,3,6]

CharacterTable1DoubleGamma4 = [[1,1,1,1,1,1,1,1,1,1],[1,-1,1,1,1,-1,-1,1,1,-1],[1,-1j,-1,1,1,1j,-1j,-1,-1,1j],[1,1j,-1,1,1,-1j,1j,-1,-1,-1j]]
CharacterTable2DoubleGamma4 = [[2,0,2,-1,2,0,0,-1,2,0],[2,0,-2,-1,2,0,0,1,-2,0]]
CharacterTable3DoubleGamma4 = [[3,1,3,0,-1,1,-1,0,-1,-1],[3,-1,3,0,-1,-1,1,0,-1,1],[3,-1j,-3,0,-1,1j,1j,0,1,-1j],[3,1j,-3,0,-1,-1j,-1j,0,1,1j]]

CharacterTableDoubleGamma4 = [CharacterTable1DoubleGamma4, CharacterTable2DoubleGamma4, CharacterTable3DoubleGamma4]

DoubleGamma4_Theory = [IrrepsDoubleGamma4, SIrrepsDoubleGamma4, TIrrepsDoubleGamma4, 4, ConjugacyClassesSizesDoubleGamma4, CharacterTableDoubleGamma4, "[48,30]", "Gamma_4'"]


"Gamma_5 = A5"
IrrepsGamma5 = [[str(1)],[],[str(3),str(3)+"'"],["4"],["5"]]

A1 = MultScal(1/5**0.5,[[1,-2**0.5,-2**0.5],[-2**0.5,-phi,phi-1],[-2**0.5,phi-1,-phi]])
A2 = MultScal(1/5**0.5,[[-1,2**0.5,2**0.5],[2**0.5,1-phi,phi],[2**0.5,phi,1-phi]])
B1 = MultScal(1/5**0.5,[[1,phi-1,phi,-1],[phi-1,-1,1,phi],[phi,1,-1,phi-1],[-1,phi,phi-1,1]])
C1 = [[-1,6**0.5,6**0.5,6**0.5,6**0.5],[6**0.5,(phi-1)**2,-2*phi,2*(phi-1),phi**2],[6**0.5,-2*phi,phi**2,(phi-1)**2,2*(phi-1)],[6**0.5,2*(phi-1),(phi-1)**2,phi**2,-2*phi],[6**0.5,phi**2,2*(phi-1),-2*phi,(phi-1)**2]]

TIrreps1Gamma5 = [[1]]
SIrreps1Gamma5 = [[1]]
TIrreps2Gamma5 = []
SIrreps2Gamma5 = [] 
TIrreps3Gamma5 = [[[1,0,0],[0,zeta,0],[0,0,zeta**4]],[[1,0,0],[0,zeta**2,0],[0,0,zeta**3]]]
SIrreps3Gamma5 = [A1,A2]
TIrreps4Gamma5 = [[[zeta,0,0,0],[0,zeta**2,0,0],[0,0,zeta**3,0],[0,0,0,zeta**4]]]
SIrreps4Gamma5 = [B1]
TIrreps5Gamma5 = [[[1,0,0,0,0],[0,zeta,0,0,0],[0,0,zeta**2,0,0],[0,0,0,zeta**3,0],[0,0,0,0,zeta**4]]]
SIrreps5Gamma5 = [MultScal(1/5, C1)]

TIrrepsGamma5 = [TIrreps1Gamma5, TIrreps2Gamma5, TIrreps3Gamma5, TIrreps4Gamma5, TIrreps5Gamma5]
SIrrepsGamma5 = [SIrreps1Gamma5, SIrreps2Gamma5, SIrreps3Gamma5, SIrreps4Gamma5, SIrreps5Gamma5]

ConjugacyClassesSizesGamma5 = [1,15,12,12,20]

CharacterTable1Gamma5 = [[1,1,1,1,1]]
CharacterTable2Gamma5 = []
CharacterTable3Gamma5 = [[3,-1,1-phi,phi,0],[3,-1,phi,1-phi,0]]
CharacterTable4Gamma5 = [[4,0,-1,-1,1]]
CharacterTable5Gamma5 = [[5,1,0,0,-1]]

CharacterTableGamma5 = [CharacterTable1Gamma5, CharacterTable2Gamma5, CharacterTable3Gamma5, CharacterTable4Gamma5, CharacterTable5Gamma5]

Gamma5_Theory = [IrrepsGamma5, SIrrepsGamma5, TIrrepsGamma5, 5, ConjugacyClassesSizesGamma5, CharacterTableGamma5, "[60,5]", "Gamma_5"]

"Double covering of Gamma5: A5'=SL(2,Z_5) called DoubleGamma5 here."
IrrepsDoubleGamma5 = [[str(1)],["2^","2^'"],[str(3),str(3)+"'"],["4","4^"],["5"],["6^"]]

A1 = MultScal(1j/5**(1/4), [[phi**0.5,(phi-1)**0.5],[(phi-1)**0.5,-phi**0.5]])
A2 = MultScal(1j/5**(1/4), [[(phi-1)**0.5,phi**0.5],[phi**0.5,-(phi-1)**0.5]])
B1 = MultScal(1/5**0.5,[[1,-2**0.5,-2**0.5],[-2**0.5,-phi,phi-1],[-2**0.5,phi-1,-phi]])
B2 = MultScal(1/5**0.5,[[-1,2**0.5,2**0.5],[2**0.5,1-phi,phi],[2**0.5,phi,1-phi]])
C1 = MultScal(1/5**0.5,[[1,phi-1,phi,-1],[phi-1,-1,1,phi],[phi,1,-1,phi-1],[-1,phi,phi-1,1]])
C2 = MultScal(1j/5**(3/4),[[-(2*phi+1)**0.5,(3*phi)**0.5,(3*(phi-1))**0.5,(2*phi-3)**0.5],[(3*phi)**0.5,(2*phi-3)**0.5,(2*phi+1)**0.5,(3*(phi-1))**0.5],[(3*(phi-1))**0.5,(2*phi+1)**0.5,-(2*phi-3)**0.5,-(3*phi)**0.5],[(2*phi-3)**0.5,(3*(phi-1))**0.5,-(3*phi)**0.5,(2*phi+1)**0.5]])
D1 = MultScal(1/5,[[-1,6**0.5,6**0.5,6**0.5,6**0.5],[6**0.5,(phi-1)**2,-2*phi,2*(phi-1),phi**2],[6**0.5,-2*phi,phi**2,(phi-1)**2,2*(phi-1)],[6**0.5,2*(phi-1),(phi-1)**2,phi**2,-2*phi],[6**0.5,phi**2,2*(phi-1),-2*phi,(phi-1)**2]])
E1 = MultScal(-1j/5**(3/4),[[phi**0.5,-(2*(phi-1))**0.5,(2*phi-3)**0.5,-(2*phi+1)**0.5,(2*phi)**0.5,(phi-1)**0.5],[-(2*(phi-1))**0.5,-(phi-1)**0.5,(2*(phi-1))**0.5,(2*phi)**0.5,phi**0.5,(2*phi)**0.5],[(2*phi-3)**0.5,(2*(phi-1))**0.5,phi**0.5,-(phi-1)**0.5,-(2*phi)**0.5,(2*phi+1)**0.5],[-(2*phi+1)**0.5,(2*phi)**0.5,-(phi-1)**0.5,-phi**0.5,(2*(phi-1))**0.5,(2*phi-3)**0.5],[(2*phi)**0.5,phi**0.5,-(2*phi)**0.5,(2*(phi-1))**0.5,(phi-1)**0.5,(2*(phi-1))**0.5],[(phi-1)**0.5,(2*phi)**0.5,(2*phi+1)**0.5,(2*phi-3)**0.5,(2*(phi-1))**0.5,-phi**0.5]])

TIrreps1DoubleGamma5 = [[1]]
SIrreps1DoubleGamma5 = [[1]]
TIrreps2DoubleGamma5 = [[[zeta**2,0],[0,zeta**3]],[[zeta,0],[0,zeta**4]]]
SIrreps2DoubleGamma5 = [A1,A2] 
TIrreps3DoubleGamma5 = [[[1,0,0],[0,zeta,0],[0,0,zeta**4]],[[1,0,0],[0,zeta**2,0],[0,0,zeta**3]]]
SIrreps3DoubleGamma5 = [B1,B2]
TIrreps4DoubleGamma5 = [[[zeta,0,0,0],[0,zeta**2,0,0],[0,0,zeta**3,0],[0,0,0,zeta**4]],[[zeta,0,0,0],[0,zeta**2,0,0],[0,0,zeta**3,0],[0,0,0,zeta**4]]]
SIrreps4DoubleGamma5 = [C1,C2]
TIrreps5DoubleGamma5 = [[[1,0,0,0,0],[0,zeta,0,0,0],[0,0,zeta**2,0,0],[0,0,0,zeta**3,0],[0,0,0,0,zeta**4]]]
SIrreps5DoubleGamma5 = [D1]
TIrreps6DoubleGamma5 = [[[1,0,0,0,0,0],[0,zeta,0,0,0,0],[0,0,zeta**2,0,0,0],[0,0,0,zeta**3,0,0],[0,0,0,0,zeta**4,0],[0,0,0,0,0,1]]]
SIrreps6DoubleGamma5 = [E1]

TIrrepsDoubleGamma5 = [TIrreps1DoubleGamma5, TIrreps2DoubleGamma5, TIrreps3DoubleGamma5, TIrreps4DoubleGamma5, TIrreps5DoubleGamma5, TIrreps6DoubleGamma5]
SIrrepsDoubleGamma5 = [SIrreps1DoubleGamma5, SIrreps2DoubleGamma5, SIrreps3DoubleGamma5, SIrreps4DoubleGamma5, SIrreps5DoubleGamma5, SIrreps6DoubleGamma5]

ConjugacyClassesSizesDoubleGamma5 = [1,12,12,30,20,20,12,1,12]

CharacterTable1DoubleGamma5 = [[1,1,1,1,1,1,1,1,1]]
CharacterTable2DoubleGamma5 = [[2,-phi,phi-1,0,1,-1,phi,-2,1-phi],[2,phi-1,-phi,0,1,-1,1-phi,-2,phi]]
CharacterTable3DoubleGamma5 = [[3,phi,1-phi,0-1,0,0,phi,3,1-phi],[3,1-phi,phi,-1,0,0,1-phi,3,phi]]
CharacterTable4DoubleGamma5 = [[4,-1,-1,0,1,1,-1,4,-1],[4,-1,-1,0,-1,1,1,-4,1]]
CharacterTable5DoubleGamma5 = [[5,0,0,1,-1,-1,0,5,0]]
CharacterTable6DoubleGamma5 = [[6,1,1,0,0,0,-1,-6,-1]]

CharacterTableDoubleGamma5 = [CharacterTable1DoubleGamma5, CharacterTable2DoubleGamma5, CharacterTable3DoubleGamma5, CharacterTable4DoubleGamma5, CharacterTable5DoubleGamma5, CharacterTable6DoubleGamma5]

DoubleGamma5_Theory = [IrrepsDoubleGamma5, SIrrepsDoubleGamma5, TIrrepsDoubleGamma5, 5, ConjugacyClassesSizesDoubleGamma5, CharacterTableDoubleGamma5, "[120,5]", "Gamma_5'"]

"Gamma_6 = Gamma_2xGamma_3=S3xA4. Representations of Gamma_6 are made from those of S3 and A4."
IrrepsGamma6 = [["1_0","1_1","1_2","1_0'","1_1'","1_2'"],[str(2),str(2)+"'",str(2)+'"'],[str(3),str(3)+"'"],[],[],["6"]]

A1 = [[1,0],[0,-1]]
A2 = [[-1/2, 3**0.5/2],[3**0.5/2,1/2]]
B1 = [[1,0,0],[0,omega,0],[0,0,omega**2]]
B2 = MultScal(1/3,[[-1,2,2],[2,-1,2],[2,2,-1]])
C1 = [[1,-2,-2,-3**0.5,2*3**0.5,2*3**0.5],[-2,1,-2,2*3**0.5,-3**0.5,2*3**0.5],[-2,-2,1,2*3**0.5,2*3**0.5,-3**0.5],[-3**0.5,2*3**0.5,2*3**0.5,-1,2,2],[2*3**0.5,-3**0.5,2*3**0.5,2,-1,2],[2*3**0.5,2*3**0.5,-3**0.5,2,2,-1]]

TIrreps1Gamma6 = [[1],[omega],[omega**2],[-1],[-omega],[-omega**2]]
SIrreps1Gamma6 = [[1],[1],[1],[-1],[-1],[-1]]
TIrreps2Gamma6 = [A1,MultScal(omega,A1),MultScal(omega**2,A1)] 
SIrreps2Gamma6 = [A2,A2,A2]
TIrreps3Gamma6 = [B1, MultScal(-1,B1)]
SIrreps3Gamma6 = [B2, MultScal(-1,B2)]
TIrreps4Gamma6 = []
SIrreps4Gamma6 = []
TIrreps5Gamma6 = []
SIrreps5Gamma6 = []
TIrreps6Gamma6 = [[[1,0,0,0,0,0],[0,omega,0,0,0,0],[0,0,omega**2,0,0,0],[0,0,0,-1,0,0],[0,0,0,0,-omega,0],[0,0,0,0,0,-omega**2]]]
SIrreps6Gamma6 = [MultScal(1/6,C1)]

TIrrepsGamma6 = [TIrreps1Gamma6, TIrreps2Gamma6, TIrreps3Gamma6, TIrreps4Gamma6, TIrreps5Gamma6, TIrreps6Gamma6]
SIrrepsGamma6 = [SIrreps1Gamma6, SIrreps2Gamma6, SIrreps3Gamma6, SIrreps4Gamma6, SIrreps5Gamma6, SIrreps6Gamma6]

ConjugacyClassesSizesGamma6 = [1,3,4,3,2,12,9,4,8,6,12,8]

CharacterTable1Gamma6 = [[1,1,1,1,1,1,1,1,1,1,1,1],[1,1,omega**2,1,1,omega**2,1,omega,omega**2,1,omega,omega],[1,1,omega,1,1,omega,1,omega**2,omega,1,omega**2,omega**2],[1,-1,1,1,1,-1,-1,1,1,1,-1,1],[1,-1,omega**2,1,1,-omega**2,-1,omega,omega**2,1,-omega,omega],[1,-1,omega,1,1,-omega,-1,omega**2,omega,1,-omega**2,omega**2]]
CharacterTable2Gamma6 = [[2,0,2,2,-1,0,0,2,-1,-1,0,-1],[2,0,2*omega**2,2,-1,0,0,2*omega,-omega**2,-1,0,-omega],[2,0,2*omega,2,-1,00,0,2*omega**2,-omega,-1,0,-omega**2]]
CharacterTable3Gamma6 = [[3,3,0,-1,3,0,-1,0,0,-1,0,0],[3,-3,0,-1,3,0,1,0,0,-1,0,0]]
CharacterTable4Gamma6 = []
CharacterTable5Gamma6 = []
CharacterTable6Gamma6 = [[6,0,0,-2,-3,0,0,0,0,1,0,0]]

CharacterTableGamma6 = [CharacterTable1Gamma6, CharacterTable2Gamma6, CharacterTable3Gamma6, CharacterTable4Gamma6, CharacterTable5Gamma6, CharacterTable6Gamma6]

Gamma6_Theory = [IrrepsGamma6, SIrrepsGamma6, TIrrepsGamma6, 6, ConjugacyClassesSizesGamma6, CharacterTableGamma6, "[72,44]", "Gamma_6"]

#==================================================
#Class of groups of the type: PermutationGroup x ZN
#==================================================

"S3xZ3 called S3Z3 here"
IrrepsS3Z3 = [["1_0","1_1","1_2","1_0'","1_1'","1_2'"],["2_0","2_1","2_2"],[]]

A1 = [[1,0],[0,-1]]
A2 = [[-1/2,-3**0.5/2],[-3**0.5/2,1/2]]

TIrreps1S3Z3 = [[1],[omega],[omega**2],[-1],[-omega],[-omega**2]]
SIrreps1S3Z3 = [[1],[1],[1],[-1],[-1],[-1]]
TIrreps2S3Z3 = [A1,MultScal(omega,A1),MultScal(omega**2,A1)] 
SIrreps2S3Z3 = [A2,A2,A2]
TIrreps3S3Z3 = []
SIrreps3S3Z3 = []

TIrrepsS3Z3 = [TIrreps1S3Z3, TIrreps2S3Z3, TIrreps3S3Z3]
SIrrepsS3Z3 = [SIrreps1S3Z3, SIrreps2S3Z3, SIrreps3S3Z3]

ConjugacyClassesSizesS3Z3 = [1,3,1,2,3,1,2,3,2]

CharacterTable1S3Z3 = [[1,1,1,1,1,1,1,1,1],[1,1,omega**2,1,omega**2,omega,omega**2,omega,omega],[1,1,omega,1,omega,omega**2,omega,omega**2,omega**2],[1,-1,1,1,-1,1,1,-1,1],[1,-1,omega**2,1,-omega**2,omega,omega**2,-omega,omega],[1,-1,omega,1,-omega,omega**2,omega,-omega**2,omega**2]]
CharacterTable2S3Z3 = [[2,0,2,-1,0,2,-1,0,-1],[2,0,2*omega**2,-1,0,2*omega,-omega**2,0,-omega],[2,0,2*omega,-1,0,2*omega**2,-omega,0,-omega**2]]
CharacterTable3S3Z3 = []

CharacterTableS3Z3 = [CharacterTable1S3Z3, CharacterTable2S3Z3, CharacterTable3S3Z3]

S3Z3_Theory = [IrrepsS3Z3, SIrrepsS3Z3, TIrrepsS3Z3, 6, ConjugacyClassesSizesS3Z3, CharacterTableS3Z3, "[18,3]", "S3xZ3"]

"A4xZ2 called A4Z2 here"
IrrepsA4Z2 = [["1_0","1_0'",'1_0"',"1_1","1_1'",'1_1"'],[],["3_0","3_1"]]

A1 = [[1,0,0],[0,omega,0],[0,0,omega**2]]
A2 = [[-1,2,2],[2,-1,2],[2,2,-1]]

TIrreps1A4Z2 = [[1],[omega],[omega**2],[-1],[-omega],[-omega**2]]
SIrreps1A4Z2 = [[1],[1],[1],[-1],[-1],[-1]]
TIrreps2A4Z2 = [] 
SIrreps2A4Z2 = []
TIrreps3A4Z2 = [A1,MultScal(-1,A1)]
SIrreps3A4Z2 = [MultScal(1/3,A2),MultScal(-1/3,A2)]

TIrrepsA4Z2 = [TIrreps1A4Z2, TIrreps2A4Z2, TIrreps3A4Z2]
SIrrepsA4Z2 = [SIrreps1A4Z2, SIrreps2A4Z2, SIrreps3A4Z2]

ConjugacyClassesSizesA4Z2 = [1,1,4,3,4,3,4,4]

CharacterTable1A4Z2 = [[1,1,1,1,1,1,1,1],[1,1,omega**2,1,omega**2,1,omega,omega],[1,1,omega,1,omega,1,omega**2,omega**2],[1,-1,1,1,-1,-1,1,-1],[1,-1,omega**2,1,-omega**2,-1,omega,-omega],[1,-1,omega,1,-omega,-1,omega**2,-omega**2]]
CharacterTable2A4Z2 = []
CharacterTable3A4Z2 = [[3,3,0,-1,0,-1,0,0],[3,-3,0,-1,0,1,0,0]]

CharacterTableA4Z2 = [CharacterTable1A4Z2, CharacterTable2A4Z2, CharacterTable3A4Z2]

A4Z2_Theory = [IrrepsA4Z2, SIrrepsA4Z2, TIrrepsA4Z2, 6, ConjugacyClassesSizesA4Z2, CharacterTableA4Z2,"[24,13]", "A4xZ2"]

"A4xZ_4 called A4Z4 here"
IrrepsA4Z4 = [["1_0","1_0'",'1_0"',"1_1","1_1'",'1_1"',"1_2","1_2'",'1_2"',"1_3","1_3'",'1_3"'],[],["3_0","3_1","3_2","3_3"]]

A1 = [[1,0,0],[0,omega,0],[0,0,omega**2]]
A2 = MultScal(1/3,[[-1,2,2],[2,-1,2],[2,2,-1]])

TIrreps1A4Z4 = [[1],[omega],[omega**2],[1j],[eta**7],[eta**11],[-1],[-omega],[-omega**2],[-1j],[-eta**7],[-eta**11]]
SIrreps1A4Z4 = [[1],[1],[1],[-1j],[-1j],[-1j],[-1],[-1],[-1],[1j],[1j],[1j]]
TIrreps2A4Z4 = [] 
SIrreps2A4Z4 = []
TIrreps3A4Z4 = [A1,MultScal(1j,A1),MultScal(-1,A1),MultScal(-1j,A1)]
SIrreps3A4Z4 = [A2,MultScal(-1j,A2),MultScal(-1,A2),MultScal(1j,A2)]

TIrrepsA4Z4 = [TIrreps1A4Z4, TIrreps2A4Z4, TIrreps3A4Z4]
SIrrepsA4Z4 = [SIrreps1A4Z4, SIrreps2A4Z4, SIrreps3A4Z4]

ConjugacyClassesSizesA4Z4 = [1,1,4,1,3,4,1,3,4,4,3,4,4,3,4,4]

CharacterTable1A4Z4 = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,omega,1,1,omega,1,1,omega**2,omega,1,omega**2,omega,1,omega**2,omega**2],[1,1,omega**2,1,1,omega**2,1,1,omega,omega**2,1,omega,omega**2,1,omega,omega],[1,1j,1,-1,1,1j,-1j,1j,1,-1,-1,1j,-1j,-1j,-1,-1j],[1,1j,omega,-1,1,-eta,-1j,1j,omega**2,-omega,-1,eta**11,eta,-1j,-omega**2,-eta**11],[1,1j,omega**2,-1,1,eta**11,-1j,1j,omega,-omega**2,-1,-eta,-eta**11,-1j,-omega,eta],[1,-1,1,1,1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1],[1,-1,omega,1,1,-omega,-1,-1,omega**2,omega,1,-omega**2,-omega,-1,omega**2,-omega**2],[1,-1,omega**2,1,1,-omega**2,-1,-1,omega,omega**2,1,-omega,-omega**2,-1,omega,-omega],[1,-1j,1,-1,1,-1j,1j,-1j,1,-1,-1,-1j,1j,1j,-1,1j],[1,-1j,omega,-1,1,eta,1j,-1j,omega**2,-omega,-1,-eta**11,-eta,1j,-omega**2,eta**11],[1,-1j,omega**2,-1,1,-eta**11,1j,-1j,omega,-omega**2,-1,eta,eta**11,1j,-omega,-eta]]
CharacterTable2A4Z4 = []
CharacterTable3A4Z4 = [[3,3,0,3,-1,0,3,-1,0,0,-1,0,0,-1,0,0],[3,3*1j,0,-3,-1,0,-3*1j,-1j,0,0,1,0,0,1j,0,0],[3,-3,0,3,-1,0,-3,1,0,0,-1,0,0,1,0,0],[3,-3*1j,0,-3,-1,0,3*1j,1j,0,0,1,0,0,-1j,0,0]]

CharacterTableA4Z4 = [CharacterTable1A4Z4, CharacterTable2A4Z4, CharacterTable3A4Z4]

A4Z4_Theory = [IrrepsA4Z4, SIrrepsA4Z4, TIrrepsA4Z4, 12, ConjugacyClassesSizesA4Z4, CharacterTableA4Z4, "[48,31]", "A4xZ4"]

"T'xZ_2 called TZ2 here"
IrrepsTZ2 = [["1_0","1_0'",'1_0"',"1_1","1_1'",'1_1"'],["2_0","2_0'",'2_0"',"2_1","2_1'",'2_1"'],["3_0","3_1"]]

A1 = [[omega,0],[0,1]]
A2 = [[-1j/3**0.5, -(2/3)**0.5*np.exp(1j*np.pi/12)],[(2/3)**0.5*np.exp(-1j*np.pi/12),1j/3**0.5]]
B1 = [[1,0,0],[0,omega,0],[0,0,omega**2]]
B2 = MultScal(1/3,[[-1,2,2],[2,-1,2],[2,2,-1]])

TIrreps1TZ2 = [[1],[omega],[omega**2],[-1],[-omega],[-omega**2]]
SIrreps1TZ2 = [[1],[1],[1],[-1],[-1],[-1]]
TIrreps2TZ2 = [A1, MultScal(omega,A1),MultScal(omega**2,A1),MultScal(-1,A1),MultScal(-omega,A1),MultScal(-omega**2,A1)] 
SIrreps2TZ2 = [A2, A2,A2,MultScal(-1,A2),MultScal(-1,A2),MultScal(-1,A2)] 
TIrreps3TZ2 = [B1, MultScal(-1,B1)]
SIrreps3TZ2 = [B2, MultScal(-1,B2)]

TIrrepsTZ2 = [TIrreps1TZ2, TIrreps2TZ2, TIrreps3TZ2]
SIrrepsTZ2 = [SIrreps1TZ2, SIrreps2TZ2, SIrreps3TZ2]

ConjugacyClassesSizesTZ2 = [1,1,4,6,1,4,6,1,4,4,4,4,4,4]

CharacterTable1TZ2 = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,omega**2,1,1,omega**2,1,1,omega,omega**2,omega,omega**2,omega,omega],[1,1,omega,1,1,omega,1,1,omega**2,omega,omega**2,omega,omega**2,omega**2],[1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,-1],[1,-1,omega**2,1,1,-omega**2,-1,-1,omega,omega**2,-omega,-omega**2,omega,-omega],[1,-1,omega,1,1,-omega,-1,-1,omega**2,omega,-omega**2,-omega,omega**2,-omega**2]]
CharacterTable2TZ2 = [[2,2,-omega,0,-2,-omega,0,-2,-omega**2,omega,-omega**2,omega,omega**2,omega**2],[2,2,-1,0,-2,-1,0,-2,-1,1,-1,1,1,1],[2,2,-omega**2,0,-2,-omega**2,0,-2,-omega,omega**2,-omega,omega**2,omega,omega],[2,-2,-omega,0,-2,omega,0,2,-omega**2,omega,omega**2,-omega,omega**2,-omega**2],[2,-2,-1,0,-2,1,0,2,-1,1,1,-1,1,-1],[2,-2,-omega**2,0,-2,omega**2,0,2,-omega,omega**2,omega,-omega**2,omega,-omega]]
CharacterTable3TZ2 = [[3,3,0,-1,3,0,-1,3,0,0,0,0,0,0],[3,-3,0,-1,3,0,1,-3,0,0,0,0,0,0]]

CharacterTableTZ2 = [CharacterTable1TZ2, CharacterTable2TZ2, CharacterTable3TZ2]

TZ2_Theory = [IrrepsTZ2, SIrrepsTZ2, TIrrepsTZ2, 6, ConjugacyClassesSizesTZ2, CharacterTableTZ2, "[48,32]", "T'xZ2"]

"S4xZ_3 called S4Z3 here"
IrrepsS4Z3 = [["1_0","1_1","1_2","1_0'","1_1'","1_2'"],["2_0","2_1","2_2"],["3_0","3_1","3_2","3_0'","3_1'","3_2'"]]

A1 = [[1,0],[0,-1]]
A2 = [[-1/2, 3**0.5/2],[3**0.5/2, 1/2]]
B1 = [[-1,0,0],[0,-1j,0],[0,0,1j]]
B2 = [[0,-1/2**0.5,-1/2**0.5],[-1/2**0.5,1/2,-1/2],[-1/2**0.5,-1/2,1/2]]  

TIrreps1S4Z3 = [[1],[omega],[omega**2],[-1],[-omega],[-omega**2]]
SIrreps1S4Z3 = [[1],[1],[1],[-1],[-1],[-1]]
TIrreps2S4Z3 = [A1, MultScal(omega,A1),MultScal(omega**2,A1)] 
SIrreps2S4Z3 = [A2,A2,A2]
TIrreps3S4Z3 = [B1, MultScal(omega,B1),MultScal(omega**2,B1),MultScal(-1,B1),MultScal(-omega,B1),MultScal(-omega**2,B1)]
SIrreps3S4Z3 = [B2,B2,B2,MultScal(-1,B2),MultScal(-1,B2),MultScal(-1,B2)]

TIrrepsS4Z3 = [TIrreps1S4Z3, TIrreps2S4Z3, TIrreps3S4Z3]
SIrrepsS4Z3 = [SIrreps1S4Z3, SIrreps2S4Z3, SIrreps3S4Z3]

ConjugacyClassesSizesS4Z3 = [1,6,1,8,3,6,6,1,8,3,6,6,8,3,6]

CharacterTable1S4Z3 = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,omega**2,1,1,omega**2,1,omega,omega**2,omega**2,omega,omega**2,omega,omega,omega],[1,1,omega,1,1,omega,1,omega**2,omega,omega,omega**2,omega,omega**2,omega**2,omega**2],[1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1],[1,-1,omega**2,1,1,-omega**2,-1,omega,omega**2,omega**2,-omega,-omega**2,omega,omega,-omega],[1,-1,omega,1,1,-omega,-1,omega**2,omega,omega,-omega**2,-omega,omega**2,omega**2,-omega**2]]
CharacterTable2S4Z3 = [[2,0,2,-1,2,0,0,2,-1,2,0,0,-1,2,0],[2,0,2*omega**2,-1,2,0,0,2*omega,-omega**2,2*omega**2,0,0,-omega,2*omega,0],[2,0,2*omega,-1,2,0,0,2*omega**2,-omega,2*omega,0,0,-omega**2,2*omega**2,0]]
CharacterTable3S4Z3 = [[3,1,3,0,-1,1,-1,3,0,-1,1,-1,0,-1,-1],[3,1,3*omega**2,0,-1,omega**2,-1,3*omega,0,-omega**2,omega,-omega**2,0,-omega,-omega],[3,1,3*omega,0,-1,omega,-1,3*omega**2,0,-omega,omega**2,-omega,0,-omega**2,-omega**2],[3,-1,3,0,-1,-1,1,3,0,-1,-1,1,0,-1,1],[3,-1,3*omega**2,0,-1,-omega**2,1,3*omega,0,-omega**2,-omega,omega**2,0,-omega,omega],[3,-1,3*omega,0,-1,-omega,1,3*omega**2,0,-omega,-omega**2,omega,0,-omega**2,omega**2]]

CharacterTableS4Z3 = [CharacterTable1S4Z3, CharacterTable2S4Z3, CharacterTable3S4Z3]

S4Z3_Theory = [IrrepsS4Z3, SIrrepsS4Z3, TIrrepsS4Z3, 12, ConjugacyClassesSizesS4Z3, CharacterTableS4Z3, "[72,42]","S4xZ3"]

"T':Z4 called SemiTZ4 here."
IrrepsSemiTZ4 = [["1_0","1_1","1_2","1_3"],["2_0","2_1","2_2","2_3","2_4","2_5"],["3_0","3_1","3_2","3_3"],["4_0","4_1"]]

TIrreps1SemiTZ4 = [[1],[-1],[1j],[-1j]]
SIrreps1SemiTZ4 = [[1],[-1],[-1j],[1j]]
TIrreps2SemiTZ4 = [[[1,0],[0,-1]],[[1j,0],[0,-1j]],[[xi,0],[0,xi**7]],[[xi**5,0],[0,xi**3]],[[xi,0],[0,xi**3]],[[xi**5,0],[0,xi**7]]]
SIrreps2SemiTZ4 = [[[-1/2,-1/2],[-3/2,1/2]],[[1j/2,-3j/2],[-1j/2,-1j/2]],[[1j/2**0.5,(-1+1j)/2],[(1+1j)/2,-1j/2**0.5]],[[-1j/2**0.5,(1+1j)/2],[(-1+1j)/2,1j/2**0.5]],[[-2**0.5/2,-(2**0.5+2)/2],[(2**0.5-2)/2,2**0.5/2]],[[2**0.5/2,1-2**0.5/2],[1+2**0.5/2,-2**0.5/2]]]
TIrreps3SemiTZ4 = [[[1,0,0],[0,1j,0],[0,0,-1j]],[[-1,0,0],[0,-1j,0],[0,0,1j]],[[-1,0,0],[0,1,0],[0,0,-1j]],[[-1,0,0],[0,1,0],[0,0,1j]]]
SIrreps3SemiTZ4 = [[[0,(1-1j)/2,(1+1j)/2],[(1+1j)/2,-1/2,1j/2],[(1-1j)/2,-1j/2,-1/2]],[[0,(-1-1j)/2,(-1+1j)/2],[(-1+1j)/2,1/2,1j/2],[(-1-1j)/2,-1j/2,1/2]],[[-1j/2,1/2,(-1-1j)/2],[-1/2,-1j/2,(1-1j)/2],[(1-1j)/2,(-1-1j)/2,0]],[[1j/2,1/2,(-1+1j)/2],[-1/2,1j/2,(1+1j)/2],[(1+1j)/2,(-1+1j)/2,0]]]
TIrreps4SemiTZ4 = [[[xi,0,0,0],[0,xi**3,0,0],[0,0,xi**5,0],[0,0,0,xi**7]],[[xi,0,0,0],[0,xi**3,0,0],[0,0,xi**5,0],[0,0,0,xi**7]]]
SIrreps4SemiTZ4 = [[[(2**0.5)/4,(1-1j)*(-(6**0.5)+(2**0.5))/8,(6**0.5)*1j/4,1j*(1-1j)*(-3*(2**0.5)+(6**0.5))/8],[1j*(1-1j)*(-(6**0.5)-(2**0.5))/8,-(2**0.5)/4,(1-1j)*((6**0.5)+3*(2**0.5))/8,-(6**0.5)*1j/4],[-(6**0.5)*1j/4,1j*(1-1j)*(-(6**0.5)+3*(2**0.5))/8,-(2**0.5)/4,(1-1j)*(-(2**0.5)+(6**0.5))/8],[(1-1j)*(-3*(2**0.5)-(6**0.5))/8,(6**0.5)*1j/4,1j*(1-1j)*((2**0.5)+(6**0.5))/8,(2**0.5)/4]],[[-(2**0.5)*1j/4,(1-1j)*(-(6**0.5)+3*(2**0.5))/8,(6**0.5)/4,1j*(1-1j)*(-(6**0.5)+(2**0.5))/8],[1j*(1-1j)*(-3*(2**0.5)-(6**0.5))/8,-(2**0.5)*1j/4,(1-1j)*(-(6**0.5)-(2**0.5))/8,(6**0.5)/4],[-(6**0.5)/4,1j*(1-1j)*(-(2**0.5)+(6**0.5))/8,(2**0.5)*1j/4,(1-1j)*(-3*(2**0.5)+(6**0.5))/8],[(1-1j)*((2**0.5)+(6**0.5))/8,-(6**0.5)/4,1j*(1-1j)*((6**0.5)+3*(2**0.5))/8,(2**0.5)*1j/4]]]

TIrrepsSemiTZ4 = [TIrreps1SemiTZ4, TIrreps2SemiTZ4, TIrreps3SemiTZ4, TIrreps4SemiTZ4]
SIrrepsSemiTZ4 = [SIrreps1SemiTZ4, SIrreps2SemiTZ4, SIrreps3SemiTZ4, SIrreps4SemiTZ4]

ConjugacyClassesSizesSemiTZ4 = [1,12,1,8,6,1,12,6,8,6,1,8,6,6,8,6]

CharacterTable1SemiTZ4 = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,-1,1,1,1,1,-1,-1,1,1,1,1,-1,-1,1,-1],[1,1j,-1,1,1,1,-1j,1j,-1,-1,-1,1,-1j,1j,-1,-1j],[1,-1j,-1,1,1,1,1j,-1j,-1,-1,-1,1,1j,-1j,-1,1j]]
CharacterTable2SemiTZ4 = [[2,0,2,-1,2,2,0,0,-1,2,2,-1,0,0,-1,0],[2,0,-2,-1,2,2,0,0,1,-2,-2,-1,0,0,1,0],[2,0,-2,-1,0,-2,0,2**0.5,1,0,2,1,-2**0.5,-2**0.5,-1,2**0.5],[2,0,-2,-1,0,-2,0,-2**0.5,1,0,2,1,2**0.5,2**0.5,-1,-2**0.5],[2,0,2,-1,0,-2,0,2**0.5*1j,-1,0,-2,1,2**0.5*1j,-2**0.5*1j,1,-2**0.5*1j],[2,0,2,-1,0,-2,0,-2**0.5*1j,-1,0,-2,1,-2**0.5*1j,2**0.5*1j,1,2**0.5*1j]]
CharacterTable3SemiTZ4 = [[3,-1,3,0,-1,3,-1,1,0,-1,3,0,1,1,0,1],[3,1,3,0,-1,3,1,-1,0,-1,3,0,-1,-1,0,-1],[3,1j,-3,0,-1,3,-1j,-1j,0,1,-3,0,1j,-1j,0,1j],[3,-1j,-3,0,-1,3,1j,1j,0,1,-3,0,-1j,1j,0,-1j]]
CharacterTable4SemiTZ4 = [[4,0,4,1,0,-4,0,0,1,0,-4,-1,0,0,-1,0],[4,0,-4,1,0,-4,0,0,-1,0,4,-1,0,0,1,0]]

CharacterTableSemiTZ4 = [CharacterTable1SemiTZ4, CharacterTable2SemiTZ4, CharacterTable3SemiTZ4, CharacterTable4SemiTZ4]

SemiTZ4_Theory = [IrrepsSemiTZ4, SIrrepsSemiTZ4, TIrrepsSemiTZ4, 8, ConjugacyClassesSizesSemiTZ4, CharacterTableSemiTZ4, "[96,66]","T':Z4"]

"T'xZ4 called TZ4 here."
IrrepsTZ4 = [["1_0","1_0'",'1_0"',"1_1","1_1'",'1_1"',"1_2","1_2'",'1_2"',"1_3","1_3'",'1_3"'],["2_0","2_0'",'2_0"',"2_1","2_1'",'2_1"',"2_2","2_2'",'2_2"',"2_3","2_3'",'2_3"'],["3_0","3_1","3_2","3_3"]]

A1 = [[omega,0],[0,1]]
A2 = [[-1j/3**0.5, -(2/3)**0.5*np.exp(1j*np.pi/12)],[(2/3)**0.5*np.exp(-1j*np.pi/12),1j/3**0.5]]
B1 = [[1,0,0],[0,omega,0],[0,0,omega**2]]
B2 = MultScal(1/3,[[-1,2,2],[2,-1,2],[2,2,-1]])

TIrreps1TZ4 = [[1],[omega],[omega**2],[1j],[1j*omega],[1j*omega**2],[-1],[-omega],[-omega**2],[-1j],[-1j*omega],[-1j*omega**2]]
SIrreps1TZ4 = [[1],[1],[1],[-1j],[-1j],[-1j],[-1],[-1],[-1],[1j],[1j],[1j]]
TIrreps2TZ4 = [A1, MultScal(omega,A1), MultScal(omega**2,A1), MultScal(1j,A1),MultScal(1j*omega,A1),MultScal(1j*omega**2,A1),MultScal(-1,A1),MultScal(-omega,A1),MultScal(-omega**2,A1),MultScal(-1j,A1),MultScal(-1j*omega,A1),MultScal(-1j*omega**2,A1)]
SIrreps2TZ4 = [A2,A2,A2,MultScal(-1j,A2),MultScal(-1j,A2),MultScal(-1j,A2),MultScal(-1,A2),MultScal(-1,A2),MultScal(-1,A2),MultScal(1j,A2),MultScal(1j,A2),MultScal(1j,A2)]
TIrreps3TZ4 = [B1,MultScal(1j,B1),MultScal(-1,B1),MultScal(-1j,B1)]
SIrreps3TZ4 = [B2,MultScal(-1j,B2),MultScal(-1,B2),MultScal(1j,B2)]

TIrrepsTZ4 = [TIrreps1TZ4, TIrreps2TZ4, TIrreps3TZ4]
SIrrepsTZ4 = [SIrreps1TZ4, SIrreps2TZ4, SIrreps3TZ4]

ConjugacyClassesSizesTZ4 = [1,1,4,1,6,1,4,1,6,1,4,4,4,6,1,4,4,4,6,1,4,4,4,4,4,4,4,4]

CharacterTable1TZ4 = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,omega,1,1,1,omega,1,1,1,omega**2,omega,omega,1,1,omega**2,omega,omega,1,1,omega**2,omega**2,omega,omega**2,omega**2,omega,omega**2,omega**2],[1,1,omega**2,1,1,1,omega**2,1,1,1,omega,omega**2,omega**2,1,1,omega,omega**2,omega**2,1,1,omega,omega,omega**2,omega,omega,omega**2,omega,omega],[1,1j,1,-1,1,1,1j,-1j,1j,1j,1,-1,1,-1,-1,1j,-1j,1j,-1j,-1j,-1,1,-1,-1j,1j,-1j,-1,-1j],[1,1j,omega,-1,1,1,eta**7,-1j,1j,1j,omega**2,-omega,omega,-1,-1,eta**11,eta,eta**7,-1j,-1j,-omega**2,omega**2,-omega,-eta**11,eta**11,eta,-omega**2,-eta**11],[1,1j,omega**2,-1,1,1,eta**11,-1j,1j,1j,omega,-omega**2,omega**2,-1,-1,eta**7,-eta**11,eta**11,-1j,-1j,-omega,omega,-omega**2,eta,eta**7,-eta**11,-omega,eta],[1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,-1,-1,-1,1,-1],[1,-1,omega,1,1,1,-omega,-1,-1,-1,omega**2,omega,omega,1,1,-omega**2,-omega,-omega,-1,-1,omega**2,omega**2,omega,-omega**2,-omega**2,-omega,omega**2,-omega**2],[1,-1,omega**2,1,1,1,-omega**2,-1,-1,-1,omega,omega**2,omega**2,1,1,-omega,-omega**2,-omega**2,-1,-1,omega,omega,omega**2,-omega,-omega,-omega**2,omega,-omega],[1,-1j,1,-1,1,1,-1j,1j,-1j,-1j,1,-1,1,-1,-1,-1j,1j,-1j,1j,1j,-1,1,-1,1j,-1j,1j,-1,1j],[1,-1j,omega,-1,1,1,eta,1j,-1j,-1j,omega**2,-omega,omega,-1,-1,-eta**11,eta**7,eta,1j,1j,-omega**2,omega**2,-omega,eta**11,-eta**11,eta**7,-omega**2,eta**11],[1,-1j,omega**2,-1,1,1,-eta**11,1j,-1j,-1j,omega,-omega**2,omega**2,-1,-1,eta,eta**11,-eta**11,1j,1j,-omega,omega,-omega**2,eta**7,eta,eta**11,-omega,eta**7]]
CharacterTable2TZ4 = [[2,2,-omega**2,2,0,-2,-omega**2,2,0,-2,-omega,-omega**2,omega**2,0,-2,-omega,-omega**2,omega**2,0,-2,-omega,omega,omega**2,-omega,omega,omega**2,omega,omega],[2,2,-1,2,0,-2,-1,2,0,-2,-1,-1,1,0,-2,-1,-1,1,0,-2,-1,1,1,-1,1,1,1,1],[2,2,-omega,2,0,-2,-omega,2,0,-2,-omega**2,-omega,omega,0,-2,-omega**2,-omega,omega,0,-2,-omega**2,omega**2,omega,-omega**2,omega**2,omega,omega**2,omega**2],[2,2j,-omega**2,-2,0,-2,-eta**11,-2j,0,-2j,-omega,omega**2,omega**2,0,2,eta,eta**11,eta**11,0,2j,omega,omega,-omega**2,eta**7,eta**7,-eta**11,-omega,eta],[2,2j,-1,-2,0,-2,-1j,-2j,0,-2j,-1,1,1,0,2,-1j,1j,1j,0,2j,1,1,-1,1j,1j,-1j,-1,-1j],[2,2j,-omega,-2,0,-2,eta,-2j,0,-2j,-omega**2,omega,omega,0,2,-eta**11,eta**7,eta**7,0,2j,omega**2,omega**2,-omega,eta**11,eta**11,eta,-omega**2,-eta**11],[2,-2,-omega**2,2,0,-2,omega**2,-2,0,2,-omega,-omega**2,omega**2,0,-2,omega,omega**2,-omega**2,0,2,-omega,omega,omega**2,omega,-omega,-omega**2,omega,-omega],[2,-2,-1,2,0,-2,1,-2,0,2,-1,-1,1,0,-2,1,1,-1,0,2,-1,1,1,1,-1,-1,1,-1],[2,-2,-omega,2,0,-2,omega,-2,0,2,-omega**2,-omega,omega,0,-2,omega**2,omega,-omega,0,2,-omega**2,omega**2,omega,omega**2,-omega**2,-omega,omega**2,-omega**2],[2,-2j,-omega**2,-2,0,-2,eta**11,2j,0,2j,-omega,omega**2,omega**2,0,2,eta**7,-eta**11,-eta**11,0,-2j,omega,omega,-omega**2,eta,eta,eta**11,-omega,eta**7],[2,-2j,-1,-2,0,-2,1j,2j,0,2j,-1,1,1,0,2,1j,-1j,-1j,0,-2j,1,1,-1,-1j,-1j,1j,-1,1j],[2,-2j,-omega,-2,0,-2,eta**7,2j,0,2j,-omega**2,omega,omega,0,2,eta**11,eta,eta,0,-2j,omega**2,omega**2,-omega,-eta**11,-eta**11,eta**7,-omega**2,eta**11]]
CharacterTable3TZ4 = [[3,3,0,3,-1,3,0,3,-1,3,0,0,0,-1,3,0,0,0,-1,3,0,0,0,0,0,0,0,0],[3,3j,0,-3,-1,3,0,-3j,-1j,3j,0,0,0,1,-3,0,0,0,1j,-3j,0,0,0,0,0,0,0,0],[3,-3,0,3,-1,3,0,-3,1,-3,0,0,0,-1,3,0,0,0,1,-3,0,0,0,0,0,0,0,0],[3,-3j,0,-3,-1,3,0,3j,1j,-3j,0,0,0,1,-3,0,0,0,-1j,3j,0,0,0,0,0,0,0,0]]

CharacterTableTZ4 = [CharacterTable1TZ4, CharacterTable2TZ4, CharacterTable3TZ4]


TZ4_Theory = [IrrepsTZ4, SIrrepsTZ4, TIrrepsTZ4, 12, ConjugacyClassesSizesTZ4, CharacterTableTZ4, "[96,69]", "T'xZ4"]

#===========================================================
#Class of groups of the type: (semi-)direct products of ZN's
#===========================================================

"(Z3:Z4)xZ3 isomorphic to 2D3xZ_3. Called Zs343 here"
IrrepsZs343 = [["1_0","1_1","1_2","1_0'","1_1'","1_2'","1_0^","1_1^","1_2^","1_0^'","1_1^'","1_2^'"],["2_0","2_1","2_2","2_0^","2_1^","2_2^"],[]]

A1 = [[1,0],[0,-1]]
A2 = [[1j,0],[0,-1j]]
B1 = MultScal(-1/2,[[1,1j*np.sqrt(3)],[-1j*np.sqrt(3),-1] ])
B2 =  MultScal(1j/2,[[1,1j*np.sqrt(3)],[-1j*np.sqrt(3),-1] ])

TIrreps1Zs343 = [[1],[omega],[omega**2],[-1],[-omega],[-omega**2],[eta],[eta**5],[-1j],[eta**7],[eta**11],[1j]]
SIrreps1Zs343 = [[1],[1],[1],[-1],[-1],[-1],[1j],[1j],[1j],[-1j],[-1j],[-1j]]
TIrreps2Zs343 = [A1, MultScal(omega,A1),MultScal(omega**2,A1),A2,MultScal(omega,A2),MultScal(omega**2, A2)] 
SIrreps2Zs343 = [B1,B1,B1,B2,B2,B2]
TIrreps3Zs343 = []
SIrreps3Zs343 = []

TIrrepsZs343 = [TIrreps1Zs343, TIrreps2Zs343, TIrreps3Zs343]
SIrrepsZs343 = [SIrreps1Zs343, SIrreps2Zs343, SIrreps3Zs343]

ConjugacyClassesSizesZs343 = [1,3,1,1,2,3,3,1,1,2,2,3,3,1,2,2,3,2]

CharacterTable1Zs343 = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,omega,1,1,omega,1,omega**2,omega,omega,1,omega**2,omega,omega**2,omega**2,omega,omega**2,omega**2],[1,1,omega**2,1,1,omega**2,1,omega,omega**2,omega**2,1,omega,omega**2,omega,omega,omega**2,omega,omega],[1,-1,1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,1],[1,-1,omega,1,1,-omega,-1,omega**2,omega,omega,1,-omega**2,-omega,omega**2,omega**2,omega,-omega**2,omega**2],[1,-1,omega**2,1,1,-omega**2,-1,omega,omega**2,omega**2,1,-omega,-omega**2,omega,omega,omega**2,-omega,omega],[1,-1j,omega,-1,1,-eta**7,1j,omega**2,-omega,omega,-1,eta**5,eta**7,-omega**2,omega**2,-omega,-eta**5,-omega**2],[1,-1j,omega**2,-1,1,eta**5,1j,omega,-omega**2,omega**2,-1,-eta**7,-eta**5,-omega,omega,-omega**2,eta**7,-omega],[1,-1j,1,-1,1,-1j,1j,1,-1,1,-1,-1j,1j,-1,1,-1,1j,-1],[1,1j,omega,-1,1,eta**7,-1j,omega**2,-omega,omega,-1,-eta**5,-eta**7,-omega**2,omega**2,-omega,eta**5,-omega**2],[1,1j,omega**2,-1,1,-eta**5,-1j,omega,-omega**2,omega**2,-1,eta**7,eta**5,-omega,omega,-omega**2,-eta**7,-omega],[1,1j,1,-1,1,1j,-1j,1,-1,1,-1,1j,-1j,-1,1,-1,-1j,-1]]
CharacterTable2Zs343 = [[2,0,2,2,-1,0,0,2,2,-1,-1,0,0,2,-1,-1,0,-1],[2,0,2*omega,2,-1,0,0,2*omega**2,2*omega,-omega,-1,0,0,2*omega**2,-omega**2,-omega,0,-omega**2],[2,0,2*omega**2,2,-1,0,0,2*omega,2*omega**2,-omega**2,-1,0,0,2*omega,-omega,-omega**2,0,-omega],[2,0,2,-2,-1,0,0,2,-2,-1,1,0,0,-2,-1,1,0,1],[2,0,2*omega,-2,-1,0,0,2*omega**2,-2*omega,-omega,1,0,0,-2*omega**2,-omega**2,omega,0,omega**2],[2,0,2*omega**2,-2,-1,0,0,2*omega,-2*omega**2,-omega**2,1,0,0,-2*omega,-omega,omega**2,0,omega]]
CharacterTable3Zs343 = []

CharacterTableZs343 = [CharacterTable1Zs343, CharacterTable2Zs343, CharacterTable3Zs343]

Zs343_Theory = [IrrepsZs343, SIrrepsZs343, TIrrepsZs343, 12, ConjugacyClassesSizesZs343, CharacterTableZs343, "[36,6]", "(Z3:Z4)xZ3"]

"Z7:Z6 called Zs76 here. This group also has a dim-6 irrep."
IrrepsZs76 = [["1_0","1_1","1_2","1_3","1_4","1_5"],[],[],[],[],["6"]]

TIrreps1Zs76 = [[1],[-1],[omega**2],[-omega**2],[omega],[-omega]]
SIrreps1Zs76 = [[1],[-1],[1],[-1],[1],[-1]]
TIrreps2Zs76 = []
SIrreps2Zs76 = []
TIrreps3Zs76 = [] 
SIrreps3Zs76 = []
TIrreps4Zs76 = []
SIrreps4Zs76 = []
TIrreps5Zs76 = []
SIrreps5Zs76 = []
TIrreps6Zs76 = [[[-1,0,0,0,0,0],[0,omega**2,0,0,0,0],[0,0,omega,0,0,0],[0,0,0,-omega,0,0],[0,0,0,0,-omega**2,0],[0,0,0,0,0,1]]]
SIrreps6Zs76 = [[[ 1/6,  0.40668889+1.70436465e-01j,-0.40668889+1.70436465e-01j,0.43525413+7.06985822e-02j,-0.42282538-1.25152466e-01j,0.4299028-9.81225088e-02j],
        [ 0.40668889-1.70436465e-01j,-1/6,0.39507823-1.95851049e-01j,0.43885788+4.29907958e-02j,0.29992799-3.23245492e-01j,0.10302748-4.28753756e-01j],
        [-0.40668889-1.70436465e-01j,0.39507823+1.95851049e-01j,-1/6,0.12997481+4.21368001e-01j,0.41405029+1.51679928e-01j,-0.27885383-3.41591840e-01j],
        [ 0.43525413-7.06985822e-02j,0.43885788-4.29907958e-02j,0.12997481-4.21368001e-01j,1/6,0.02792713+4.40073312e-01j,-0.20035179+3.92814975e-01j],
        [-0.42282538+1.25152466e-01j,0.29992799+3.23245492e-01j,0.41405029-1.51679928e-01j,0.02792713-4.40073312e-01j,1/6,0.43203241+8.82747947e-02j],
        [ 0.4299028 +9.81225088e-02j,0.10302748+4.28753756e-01j,-0.27885383+3.41591840e-01j,-0.20035179-3.92814975e-01j,0.43203241-8.82747947e-02j, -1/6]]]

TIrrepsZs76 = [TIrreps1Zs76, TIrreps2Zs76, TIrreps3Zs76, TIrreps4Zs76, TIrreps5Zs76, TIrreps6Zs76]
SIrrepsZs76 = [SIrreps1Zs76, SIrreps2Zs76, SIrreps3Zs76, SIrreps4Zs76, SIrreps5Zs76, SIrreps6Zs76]

ConjugacyClassesSizesZs76 = [1,7,7,6,7,7,7]

CharacterTable1Zs76 = [[1,1,1,1,1,1,1],[1,-1,1,1,-1,1,-1],[1,1,omega,1,omega,omega**2,omega**2],[1,-1,omega,1,-omega,omega**2,-omega**2],[1,1,omega**2,1,omega**2,omega,omega],[1,-1,omega**2,1,-omega**2,omega,-omega]]
CharacterTable2Zs76 = []
CharacterTable3Zs76 = []
CharacterTable4Zs76 = []
CharacterTable5Zs76 = []
CharacterTable6Zs76 = [[6,0,0,-1,0,0,0]]

CharacterTableZs76 = [CharacterTable1Zs76,CharacterTable2Zs76,CharacterTable3Zs76,CharacterTable4Zs76,CharacterTable5Zs76,CharacterTable6Zs76]

Zs76_Theory = [IrrepsZs76, SIrrepsZs76, TIrrepsZs76, 6, ConjugacyClassesSizesZs76, CharacterTableZs76, "[42,1]", "Z7:Z6"]

"((Z4xZ2):Z2):Z3 called Zs4223 here."
IrrepsZs4223 = [["1_0","1_1","1_2","1_3","1_4","1_5"],["2_0","2_1","2_2","2_3","2_4","2_5"],["3","3'"]]

z = 0.21132487-0.788675135j 
z4 = 0.78867513+1.36602540e+00j
z4bis =  0.21132487-3.66025404e-01j
B1 = [[-1,-1-3**0.5*1j,1-3**0.5*1j],[-1+3**0.5*1j,-1,1+3**0.5*1j],[1+3**0.5*1j,1-3**0.5*1j,-1]]
B2 = [[1,1+3**0.5*1j,1-3**0.5*1j],[1-3**0.5*1j,1,1+3**0.5*1j],[1+3**0.5*1j,1-3**0.5*1j,1]]

TIrreps1Zs4223 = [[1],[-1],[-omega],[-omega**2],[omega**2],[omega]]
SIrreps1Zs4223 = [[1],[-1],[-1],[-1],[1],[1]]
TIrreps2Zs4223 = [[[eta**11,0],[0,eta**7]],[[eta,0],[0,eta**5]],[[eta**5,0],[0,eta**9]],[[eta**3,0],[0,eta**7]],[[eta**3,0],[0,eta**11]],[[eta,0],[0,eta**9]]]
SIrreps2Zs4223 = [[[-3**0.5/3,z],[np.conj(z),3**0.5/3]],[[-3**0.5/3,1j*z],[np.conj(1j*z),3**0.5/3]],[[-3**0.5/3,1j*np.conj(z)],[-1j*z,3**0.5/3]],[[3**0.5/3,z4],[z4bis,-3**0.5/3]],[[-3**0.5/3,-z],[-np.conj(z),3**0.5/3]],[[3**0.5/3,z],[np.conj(z),-3**0.5/3]]]
TIrreps3Zs4223 = [[[eta**4,0,0],[0,eta**8,0],[0,0,1]],[[eta**6,0,0],[0,eta**2,0],[0,0,eta**10]]]
SIrreps3Zs4223 = [MultScal(1/3,B1),MultScal(1/3,B2)]

TIrrepsZs4223 = [TIrreps1Zs4223, TIrreps2Zs4223, TIrreps3Zs4223]
SIrrepsZs4223 = [SIrreps1Zs4223, SIrreps2Zs4223, SIrreps3Zs4223]

ConjugacyClassesSizesZs4223 = [1,1,4,6,1,4,6,1,4,4,4,4,4,4]

CharacterTable1Zs4223 = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,-1],[1,-1,omega**2,1,1,-omega**2,-1,-1,omega,omega**2,-omega,-omega**2,omega,-omega],[1,-1,omega,1,1,-omega,-1,-1,omega**2,omega,-omega**2,-omega,omega**2,-omega**2],[1,1,omega,1,1,omega,1,1,omega**2,omega,omega**2,omega,omega**2,omega**2],[1,1,omega**2,1,1,omega**2,1,1,omega,omega**2,omega,omega**2,omega,omega]]
CharacterTable2Zs4223 = [[2,2j,-1,0,-2,-1j,0,-2j,-1,1,-1j,1j,1,1j],[2,-2j,-1,0,-2,1j,0,+2j,-1,1,1j,-1j,1,-1j],[2,-2j,-omega**2,0,-2,-eta**5,0,2j,-omega,omega**2,eta**7,eta**5,omega,-eta**7],[2,2j,-omega,0,-2,-eta**7,0,-2j,-omega**2,omega,eta**5,eta**7,omega**2,-eta**5],[2,2j,-omega**2,0,-2,eta**5,0,-2j,-omega,omega**2,-eta**7,-eta**5,omega,eta**7],[2,-2j,-omega,0,-2,eta**7,0,2j,-omega**2,omega,-eta**5,-eta**7,omega**2,eta**5]]
CharacterTable3Zs4223 = [[3,3,0,-1,3,0,-1,3,0,0,0,0,0,0],[3,-3,0,-1,3,0,1,-3,0,0,0,0,0,0]]

CharacterTableZs4223 = [CharacterTable1Zs4223, CharacterTable2Zs4223, CharacterTable3Zs4223]

Zs4223_Theory = [IrrepsZs4223, SIrrepsZs4223, TIrrepsZs4223, 12, ConjugacyClassesSizesZs4223, CharacterTableZs4223, "[48,33]", "((Z4xZ2):Z2):Z3"]

"(Z3xZ3):Z6 called Zs336 here. This group also has dim-6 irrep."
IrrepsZs336 = [["1_0","1_1","1_2","1_3","1_4","1_5"],["2", "2'",'2"'],[],[],[],["6"]]

A1 = [[1,0],[0,-1]]
A2 = MultScal(-1/2,[[1,1j*3**0.5],[-1j*3**0.5,-1]])
z = 1/2*np.exp(2*1j*np.pi/3)

TIrreps1Zs336 = [[1],[-1],[omega**2],[-omega**2],[omega],[-omega]]
SIrreps1Zs336 = [[1],[-1],[1],[-1],[1],[-1]]
TIrreps2Zs336 = [A1,MultScal(omega**2,A1),MultScal(omega,A1)] 
SIrreps2Zs336 = [A2,A2,A2]
TIrreps3Zs336 = []
SIrreps3Zs336 = []
TIrreps4Zs336 = []
SIrreps4Zs336 = []
TIrreps5Zs336 = []
SIrreps5Zs336 = []
TIrreps6Zs336 = [[[-1,0,0,0,0,0],[0,omega**2,0,0,0,0],[0,0,omega,0,0,0],[0,0,0,-omega,0,0],[0,0,0,0,-omega**2,0],[0,0,0,0,0,1]]]
SIrreps6Zs336 = [[[0,np.conj(z),-z,np.conj(z),1/2,0],[z,0,np.conj(z),1/2,0,1/2],[-np.conj(z),z,0,0,-np.conj(z),z],[z,1/2,0,0,-z,-1/2],[ 1/2,0,-z,-np.conj(z),0,np.conj(z)],[0,1/2,np.conj(z),-1/2,z,0]]]

TIrrepsZs336 = [TIrreps1Zs336, TIrreps2Zs336, TIrreps3Zs336, TIrreps4Zs336, TIrreps5Zs336, TIrreps6Zs336]
SIrrepsZs336 = [SIrreps1Zs336, SIrreps2Zs336, SIrreps3Zs336, SIrreps4Zs336, SIrreps5Zs336, SIrreps6Zs336]

ConjugacyClassesSizesZs336 = [1,9,3,6,2,9,3,6,9,6]

CharacterTable1Zs336 = [[1,1,1,1,1,1,1,1,1,1],[1,-1,1,1,1,-1,1,1,-1,1],[1,1,omega,1,1,omega,omega**2,omega,omega**2,omega**2],[1,-1,omega,1,1,-omega,omega**2,omega,-omega**2,omega**2],[1,1,omega**2,1,1,omega**2,omega,omega**2,omega,omega],[1,-1,omega**2,1,1,-omega**2,omega,omega**2,-omega,omega]]
CharacterTable2Zs336 = [[2,0,2,-1,2,0,2,-1,0,-1],[2,0,2*omega,-1,2,0,2*omega**2,-omega,0,-omega**2],[2,0,2*omega**2,-1,2,0,2*omega,-omega**2,0,-omega]]
CharacterTable3Zs336 = []
CharacterTable4Zs336 = []
CharacterTable5Zs336 = []
CharacterTable6Zs336 = [[6,0,0,0,-3,0,0,0,0,0]]

CharacterTableZs336 = [CharacterTable1Zs336, CharacterTable2Zs336, CharacterTable3Zs336, CharacterTable4Zs336, CharacterTable5Zs336, CharacterTable6Zs336]


Zs336_Theory = [IrrepsZs336, SIrrepsZs336, TIrrepsZs336, 6, ConjugacyClassesSizesZs336, CharacterTableZs336, "[54,5]","(Z3xZ3):Z6"]

"Z13:Z6 called Zs136"
IrrepsZs136 = [["1_0","1_1","1_2","1_3","1_4","1_5"],[],[],[],[],["6_0","6_1"]]

TIrreps1Zs136 = [[1],[-1],[omega**2],[-omega**2],[omega],[-omega]]
SIrreps1Zs136 = [[1],[-1],[1],[-1],[1],[-1]]
TIrreps2Zs136 = []
SIrreps2Zs136 = []
TIrreps3Zs136 = []
SIrreps3Zs136 = []
TIrreps4Zs136 = []
SIrreps4Zs136 = []
TIrreps5Zs136 = []
SIrreps5Zs136 = []
TIrreps6Zs136 = [[[1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,omega,0,0,0],[0,0,0,-omega,0,0],[0,0,0,0,omega**2,0],[0,0,0,0,0,-omega**2]],
[[1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,omega,0,0,0],[0,0,0,-omega,0,0],[0,0,0,0,omega**2,0],[0,0,0,0,0,-omega**2]]]
SIrreps6Zs136 = S_new_list = [[[0.21713,-0.17414j,-0.04830-0.35797j,-0.40407-0.40929j,-0.04830+0.35797j,0.40407-0.40929j],[0.17414j,-0.21713,0.40407+0.40929j,0.04830+0.35797j,-0.40407+0.40929j,0.04830-0.35797j],[-0.04830+0.35797j,0.40407-0.40929j,0.21713,-0.17414j,-0.04830-0.35797j,-0.40407-0.40929j],[-0.40407+0.40929j,0.04830-0.35797j,0.17414j,-0.21713,0.40407+0.40929j,0.04830+0.35797j],[-0.04830-0.35797j,-0.40407-0.40929j,-0.04830+0.35797j,0.40407-0.40929j,0.21713,-0.17414j],[0.40407+0.40929j,0.04830+0.35797j,-0.40407+0.40929j,0.04830-0.35797j,0.17414j,-0.21713]],[[-0.38380,-0.57514j,-0.29357-0.38006j,-0.04615+0.16791j,-0.29357+0.38006j,0.04615+0.16791j],[0.57514j,0.38380,0.04615-0.16791j,0.29357+0.38006j,-0.04615-0.16791j,0.29357-0.38006j],[-0.29357+0.38006j,0.04615+0.16791j,-0.38380,-0.57514j,-0.29357-0.38006j,-0.04615+0.16791j],[-0.04615-0.16791j,0.29357-0.38006j,0.57514j,0.38380,0.04615-0.16791j,0.29357+0.38006j],[-0.29357-0.38006j,-0.04615+0.16791j,-0.29357+0.38006j,0.04615+0.16791j,-0.38380,-0.57514j],[0.04615-0.16791j,0.29357+0.38006j,-0.04615-0.16791j,0.29357-0.38006j,0.57514j,0.38380]]]

TIrrepsZs136 = [TIrreps1Zs136, TIrreps2Zs136, TIrreps3Zs136, TIrreps4Zs136, TIrreps5Zs136, TIrreps6Zs136]
SIrrepsZs136 = [SIrreps1Zs136, SIrreps2Zs136, SIrreps3Zs136, SIrreps4Zs136, SIrreps5Zs136, SIrreps6Zs136]

ConjugacyClassesSizesZs136 = [1,13,13,6,13,13,6,13]

CharacterTable1Zs136 = [[1,1,1,1,1,1,1,1],[1,-1,1,1,-1,1,1,-1],[1,1,omega,1,omega,omega**2,1,omega**2],[1,-1,omega,1,-omega,omega**2,1,-omega**2],[1,1,omega**2,1,omega**2,omega,1,omega],[1,-1,omega**2,1,-omega**2,omega,1,-omega]]
CharacterTable2Zs136 = []
CharacterTable3Zs136 = []
CharacterTable4Zs136 = []
CharacterTable5Zs136 = []
CharacterTable6Zs136 = [[6,0,0,(-1+13**0.5)/2,0,0,(-1-13**0.5)/2,0],[6,0,0,(-1-13**0.5)/2,0,0,(-1+13**0.5)/2,0]]

CharacterTableZs136 = [CharacterTable1Zs136, CharacterTable2Zs136, CharacterTable3Zs136, CharacterTable4Zs136, CharacterTable5Zs136, CharacterTable6Zs136]

Zs136_Theory = [IrrepsZs136, SIrrepsZs136, TIrrepsZs136, 6, ConjugacyClassesSizesZs136, CharacterTableZs136, "[78,1]","Z13:Z6"]

"(((Z4xZ4):Z2):Z3) called Zs4423 here."
IrrepsZs4423 = [["1_0","1_1","1_2","1_3","1_4","1_5"],[],["3_0","3_1"],[],[],["6_0","6_1"]]

TIrreps1Zs4423 = [[1],[-1],[omega**2],[-omega**2],[omega],[-omega]]
SIrreps1Zs4423 = [[1],[-1],[1],[-1],[1],[-1]]
TIrreps2Zs4423 = []
SIrreps2Zs4423 = []
TIrreps3Zs4423 = [[[1,0,0],[0,omega,0],[0,0,omega**2]],[[-1,0,0],[0,-omega,0],[0,0,-omega**2]]]
SIrreps3Zs4423 = [[[-1/3,2*omega/3,2*omega**2/3],[2*omega**2/3,-1/3,2*omega/3],[2*omega/3,2*omega**2/3,-1/3]],[[1/3,-2*omega/3,-2*omega**2/3],[-2*omega**2/3,1/3,-2*omega/3],[-2*omega/3,-2*omega**2/3,1/3]]]
TIrreps4Zs4423 = []
SIrreps4Zs4423 = []
TIrreps5Zs4423 = []
SIrreps5Zs4423 = []
TIrreps6Zs4423 = [[[1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,omega,0,0,0],[0,0,0,omega**2,0,0],[0,0,0,0,-omega,0],[0,0,0,0,0,-omega**2]],[[1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,omega,0,0,0],[0,0,0,omega**2,0,0],[0,0,0,0,-omega,0],[0,0,0,0,0,-omega**2]]]
SIrreps6Zs4423 = [[[1/3,0,-1/6+1j*3**0.5/6,-1/6-1j*3**0.5/6,3**0.5/6-1j/2,-3**0.5/6-1j/2],[0,-1/3,-3**0.5/6+1j/2,3**0.5/6+1j/2,1/6-1j*3**0.5/6,1/6+1j*3**0.5/6],[-1/6-1j*3**0.5/6,-3**0.5/6-1j/2,1/3,-1/6+1j*3**0.5/6,0,3**0.5/6-1j/2],[-1/6+1j*3**0.5/6,3**0.5/6-1j/2,-1/6-1j*3**0.5/6,1/3,-3**0.5/6-1j/2,0],[3**0.5/6+1j/2,1/6+1j*3**0.5/6,0,-3**0.5/6+1j/2,-1/3,1/6-1j*3**0.5/6],[-3**0.5/6+1j/2,1/6-1j*3**0.5/6,3**0.5/6+1j/2,0,1/6+1j*3**0.5/6,-1/3]],[[-1/3,-2j/3,1/6-1j*3**0.5/6,1/6+1j*3**0.5/6,-3**0.5/6-1j/6,3**0.5/6-1j/6],[2j/3,1/3,3**0.5/6+1j/6,-3**0.5/6+1j/6,-1/6+1j*3**0.5/6,-1/6-1j*3**0.5/6],[1/6+1j*3**0.5/6,3**0.5/6-1j/6,-1/3,1/6-1j*3**0.5/6,-2j/3,-3**0.5/6-1j/6],[1/6-1j*3**0.5/6,-3**0.5/6-1j/6,1/6+1j*3**0.5/6,-1/3,3**0.5/6-1j/6,-2j/3],[-3**0.5/6+1j/6,-1/6-1j*3**0.5/6,2j/3,3**0.5/6+1j/6,1/3,-1/6+1j*3**0.5/6],[3**0.5/6+1j/6,-1/6+1j*3**0.5/6,-3**0.5/6+1j/6,2j/3,-1/6-1j*3**0.5/6,1/3]]]

TIrrepsZs4423 = [TIrreps1Zs4423, TIrreps2Zs4423, TIrreps3Zs4423, TIrreps4Zs4423, TIrreps5Zs4423, TIrreps6Zs4423]
SIrrepsZs4423 = [SIrreps1Zs4423, SIrreps2Zs4423, SIrreps3Zs4423, SIrreps4Zs4423, SIrreps5Zs4423, SIrreps6Zs4423]

ConjugacyClassesSizesZs4423 = [1,4,16,6,3,16,12,16,6,16]

CharacterTable1Zs4423 = [[1,1,1,1,1,1,1,1,1,1],[1,-1,1,1,1,-1,-1,1,1,-1],[1,1,omega,1,1,omega,1,omega**2,1,omega**2],[1,-1,omega,1,1,-omega,-1,omega**2,1,-omega**2],[1,1,omega**2,1,1,omega**2,1,omega,1,omega],[1,-1,omega**2,1,1,-omega**2,-1,omega,1,-omega]]
CharacterTable2Zs4423 = []
CharacterTable3Zs4423 = [[3,3,0,-1,3,0,-1,0,-1,0],[3,-3,0,-1,3,0,1,0,-1,0]]
CharacterTable4Zs4423 = []
CharacterTable5Zs4423 = []
CharacterTable6Zs4423 = [[6,0,0,2,-2,0,0,0,-2,0],[6,0,0,-2,-2,0,0,0,2,0]]

CharacterTableZs4423 = [CharacterTable1Zs4423, CharacterTable2Zs4423, CharacterTable3Zs4423, CharacterTable4Zs4423, CharacterTable5Zs4423, CharacterTable6Zs4423]

Zs4423_Theory = [IrrepsZs4423, SIrrepsZs4423, TIrrepsZs4423, 6, ConjugacyClassesSizesZs4423, CharacterTableZs4423, "[96,72]","(((Z4xZ4):Z2):Z3)"]

#=======================
#Class of special groups
#=======================

"2D_3 the binary dihedral group. We call it BD3 here"
IrrepsBD3 = [[str(1),str(1)+"'",str(1)+"^",str(1)+"^'"],[str(2),str(2)+"^"],[] ]

A1 = [[1,0],[0,-1]]
A2 = [[-1/2,3**0.5/2],[3**0.5/2,1/2]]

TIrreps1BD3 = [[1],[-1],[-1j],[1j]]
SIrreps1BD3 = [[1],[-1],[1j],[-1j]]
TIrreps2BD3 = [A1,MultScal(-1j,A1)] 
SIrreps2BD3 = [A2,MultScal(1j,A2)]
TIrreps3BD3 = []
SIrreps3BD3 = []

TIrrepsBD3 = [TIrreps1BD3, TIrreps2BD3, TIrreps3BD3]
SIrrepsBD3 = [SIrreps1BD3, SIrreps2BD3, SIrreps3BD3]

ConjugacyClassesSizesBD3 = [1,3,1,2,3,2]

CharacterTable1BD3 = [[1,1,1,1,1,1],[1,-1,1,1,-1,1],[1,-1j,-1,1,1j,-1],[1,1j,-1,1,-1j,-1]]
CharacterTable2BD3 = [[2,0,2,-1,0,-1],[2,0,-2,-1,0,1]]
CharacterTable3BD3 = []

CharacterTableBD3 = [CharacterTable1BD3, CharacterTable2BD3, CharacterTable3BD3]

BD3_Theory = [IrrepsBD3, SIrrepsBD3, TIrrepsBD3, 4, ConjugacyClassesSizesBD3, CharacterTableBD3, "[12,1]", "2D3"]

"2O the binary octahedral group. We call it BO here. This group also has a dim-4 irrep."
IrrepsBO = [[str(1),str(1)+"'"],[str(2), str(2)+"^", str(2)+"^'"],[str(3),str(3)+"'"],["4^"]]

A1 = [[-1/2,-3**0.5/2],[-3**0.5/2,1/2]]
A2 = MultScal(1j/2**0.5,[[1,1],[1,-1]])
B1 = [[1,0,0],[0,xi**2,0],[0,0,xi**6]]
B2 = MultScal(1/2, [[0,2**0.5,2**0.5],[2**0.5,-1,1],[2**0.5,1,-1]])
C1 = [[xi**5,0,0,0],[0,xi**3,0,0],[0,0,xi**7,0],[0,0,0,xi]]
C2 = MultScal(1j/(2*2**0.5),[[1,1,3**0.5,-3**0.5],[1,-1,-3**0.5,-3**0.5],[3**0.5,-3**0.5,1,1],[-3**0.5,-3**0.5,1,-1]])

TIrreps1BO = [[1],[-1]]
SIrreps1BO = [[1],[-1]]
TIrreps2BO = [[[1,0],[0,-1]],[[xi**3,0],[0,xi**5]],[[xi,0],[0,xi**7]]] 
SIrreps2BO = [A1,A2,A2]
TIrreps3BO = [B1,MultScal(-1,B1)]
SIrreps3BO = [B2,MultScal(-1,B2)]
TIrreps4BO = [C1]
SIrreps4BO = [C2]

TIrrepsBO = [TIrreps1BO, TIrreps2BO, TIrreps3BO, TIrreps4BO]
SIrrepsBO = [SIrreps1BO, SIrreps2BO, SIrreps3BO, SIrreps4BO]

ConjugacyClassesSizesBO = [1,1,6,6,6,8,8,12]

CharacterTable1BO = [[1,1,1,1,1,1,1,1],[1,1,-1,-1,1,1,1,-1]]
CharacterTable2BO = [[2,2,0,0,2,-1,-1,0],[2,-2,-2**0.5,2**0.5,0,1,-1,0],[2,-2,2**0.5,-2**0.5,0,1,-1,0]]
CharacterTable3BO = [[3,3,1,1,-1,0,0,-1],[3,3,-1,-1,-1,0,0,1]]
CharacterTable4BO = [[4,-4,0,0,0,-1,1,0]]

CharacterTableBO = [CharacterTable1BO, CharacterTable2BO, CharacterTable3BO, CharacterTable4BO]

BO_Theory = [IrrepsBO, SIrrepsBO, TIrrepsBO, 8, ConjugacyClassesSizesBO, CharacterTableBO, "[48,28]", "2O"]

"GL(2,3) called GL2_3 here. This group also has a dim-4 irrep."
IrrepsGL2_3 = [[str(1), str(1)+"'"], [str(2),str(2)+"'",str(2)+'"'],[str(3), str(3)+"'"],["4"]]

A1 = [[-1/2,-3**0.5/2],[-3**0.5/2,1/2]]
A2 = MultScal(1/2**0.5,[[1,1],[1,-1]])
A3 = MultScal(1/2**0.5,[[-1,1],[1,1]])
B1 = [[1,0,0],[0,1j,0],[0,0,-1j]]
B2 = MultScal(1/2,[[0,2**0.5,2**0.5],[2**0.5,-1,1],[2**0.5,1,-1]])
C1 = [[xi**5,0,0,0],[0,xi**3,0,0],[0,0,xi**7,0],[0,0,0,xi]]
C2 = MultScal(1/(2*2**0.5),[[-1,3**0.5,1,3**0.5],[3**0.5,-1,3**0.5,1],[1,3**0.5,1,-3**0.5],[3**0.5,1,-3**0.5,1]])

TIrreps1GL2_3 = [[1],[-1]]
SIrreps1GL2_3 = [[1],[-1]]
TIrreps2GL2_3 = [[[1,0],[0,-1]],[[xi**3,0],[0,xi]],[[xi**7,0],[0,xi**5]]] 
SIrreps2GL2_3 = [A1,A2,A3]
TIrreps3GL2_3 = [B1,MultScal(-1,B1)]
SIrreps3GL2_3 = [B2,MultScal(-1,B2)]
TIrreps4GL2_3 = [C1]
SIrreps4GL2_3 = [C2]

TIrrepsGL2_3 = [TIrreps1GL2_3, TIrreps2GL2_3, TIrreps3GL2_3, TIrreps4GL2_3]
SIrrepsGL2_3 = [SIrreps1GL2_3, SIrreps2GL2_3, SIrreps3GL2_3, SIrreps4GL2_3]

ConjugacyClassesSizesGL2_3 = [1,1,6,6,6,8,8,12]

CharacterTable1GL2_3 = [[1,1,1,1,1,1,1,1],[1,1,-1,-1,1,1,1,-1]]
CharacterTable2GL2_3 = [[2,2,0,0,2,-1,-1,0],[2,-2,2**0.5*1j,-2**0.5*1j,0,-1,1,0],[2,-2,-2**0.5*1j,2**0.5*1j,0,-1,1,0]]
CharacterTable3GL2_3 = [[3,3,1,1,-1,0,0,-1],[3,3,-1,-1,-1,0,0,1]]
CharacterTable4GL2_3 = [[4,-4,0,0,0,1,-1,0]] 

CharacterTableGL2_3 = [CharacterTable1GL2_3, CharacterTable2GL2_3, CharacterTable3GL2_3, CharacterTable4GL2_3]

GL2_3_Theory = [IrrepsGL2_3, SIrrepsGL2_3, TIrrepsGL2_3, 8, ConjugacyClassesSizesGL2_3, CharacterTableGL2_3, "[48,29]", "GL(2,3)"]

"Δ(96) called Delta96 here."
IrrepsDelta96 = [["1","1'"],["2"],["3_0","3_1","3_2","3_3","3_4","3_5"],[],[],["6"]]

TIrreps1Delta96 = [[1],[-1]]
SIrreps1Delta96 = [[1],[-1]]
TIrreps2Delta96 = [[[1,0],[0,-1]]]
SIrreps2Delta96 = [[[-1/2,-1j*3**0.5/2],[1j*3**0.5/2,1/2]]]
TIrreps3Delta96 = [[[-1,0,0],[0,-1j,0],[0,0,1j]],[[1,0,0],[0,-1j,0],[0,0,1j]],[[-xi,0,0],[0,xi,0],[0,0,1j]],[[-xi,0,0],[0,xi,0],[0,0,-1j]],[[-1j,0,0],[0,-1/xi,0],[0,0,1/xi]],[[1j,0,0],[0,-1/xi,0],[0,0,1/xi]]]
SIrreps3Delta96 = [[[0,-1,-1],[-1/2,1/2,-1/2],[-1/2,-1/2,1/2]],[[0,1,1],[1/2,-1/2,1/2],[1/2,1/2,-1/2]],[[-1/2,1/2,-1j/2],[1/2,-1/2,-1j/2],[1j,1j,0]],[[1/2,-1/2,1j/2],[-1/2,1/2,1j/2],[-1j,-1j,0]],[[0,-1j,-1j],[1j/2,-1/2,1/2],[1j/2,1/2,-1/2]],[[0,1j,1j],[-1j/2,1/2,-1/2],[-1j/2,-1/2,1/2]]]
TIrreps4Delta96 = []
SIrreps4Delta96 = []
TIrreps5Delta96 = []
SIrreps5Delta96 = []
TIrreps6Delta96 = [[[1/xi,0,0,0,0,0],[0,-1/xi,0,0,0,0],[0,0,1,0,0,0],[0,0,0,-1,0,0],[0,0,0,0,xi,0],[0,0,0,0,0,-xi]]]
SIrreps6Delta96 = [[[0,0,-1j/2,1j/2,(1-1j)/(2*2**0.5),(1-1j)/(2*2**0.5)],[0,0,-1j/2,1j/2,(-1+1j)/(2*2**0.5),(-1+1j)/(2*2**0.5)],[1j/2,1j/2,0,0,(1+1j)/(2*2**0.5),(-1-1j)/(2*2**0.5)],[-1j/2,-1j/2,0,0,(1+1j)/(2*2**0.5),(-1-1j)/(2*2**0.5)],[(1+1j)/(2*2**0.5),-(1+1j)/(2*2**0.5),(1-1j)/(2*2**0.5),(1-1j)/(2*2**0.5),0,0],[(1+1j)/(2*2**0.5),-(1+1j)/(2*2**0.5),(-1+1j)/(2*2**0.5),(-1+1j)/(2*2**0.5),0,0]]]

TIrrepsDelta96 = [TIrreps1Delta96, TIrreps2Delta96, TIrreps3Delta96, TIrreps4Delta96, TIrreps5Delta96, TIrreps6Delta96]
SIrrepsDelta96 = [SIrreps1Delta96, SIrreps2Delta96, SIrreps3Delta96, SIrreps4Delta96, SIrreps5Delta96, SIrreps6Delta96]

ConjugacyClassesSizesDelta96 = [1,12,32,6,3,12,12,3,3,12]

CharacterTable1Delta96 = [[1,1,1,1,1,1,1,1,1,1],[1,-1,1,1,1,-1,-1,1,1,-1]]
CharacterTable2Delta96 = [[2,0,-1,2,2,0,0,2,2,0]]
CharacterTable3Delta96 = [[3,1,0,-1,3,-1,1,-1,-1,-1],[3,-1,0,-1,3,1,-1,-1,-1,1],[3,-1,0,1,-1,1j,1,-1-2j,-1+2j,-1j],[3,1,0,1,-1,-1j,-1,-1-2j,-1+2j,1j],[3,-1,0,1,-1,-1j,1,-1+2j,-1-2j,1j],[3,1,0,1,-1,1j,-1,-1+2j,-1-2j,-1j]]
CharacterTable4Delta96 = []
CharacterTable5Delta96 = []
CharacterTable6Delta96 = [[6,0,0,-2,-2,0,0,2,2,0]]

CharacterTableDelta96 = [CharacterTable1Delta96, CharacterTable2Delta96, CharacterTable3Delta96, CharacterTable4Delta96, CharacterTable5Delta96, CharacterTable6Delta96]

Delta96_Theory = [IrrepsDelta96, SIrrepsDelta96, TIrrepsDelta96, 8, ConjugacyClassesSizesDelta96, CharacterTableDelta96, "[96,64]", "Δ(96)"]

"Δ(384) called Delta384 here."
IrrepsDelta384 = [["1","1'"],["2"],["3_0","3_1","3_2","3_3","3_4","3_5","3_6","3_7","3_8","3_9","3_10","3_11","3_12","3_13"],[],[],["6_0","6_1","6_2","6_3","6_4","6_5","6_6"]]

TIrreps1Delta384 = [[1],[-1]]
SIrreps1Delta384 = [[1],[-1]]
TIrreps2Delta384 = [[[1,0],[0,-1]]]
SIrreps2Delta384 = [[[-1/2,1j*3**0.5/2],[-1j*3**0.5/2,1/2]]]
TIrreps3Delta384 = [[[-1,0,0],[0,1j,0],[0,0,-1j]],[[1,0,0],[0,1j,0],[0,0,-1j]],[[1j,0,0],[0,xi,0],[0,0,-xi]],[[-1j,0,0],[0,xi,0],[0,0,-xi]],[[-1j,0,0],[0,xi**3,0],[0,0,-xi**3]],[[1j,0,0],[0,xi**3,0],[0,0,-xi**3]],[[xi**3,0,0],[0,(-xi)**0.5,0],[0,0,-(-xi)**0.5]],[[-xi**3,0,0],[0,(-xi)**0.5,0],[0,0,-(-xi)**0.5]],[[-xi**3,0,0],[0,xi**0.5,0],[0,0,-xi**0.5]],[[xi**3,0,0],[0,xi**0.5,0],[0,0,-xi**0.5]],[[xi,0,0],[0,(-xi**3)**0.5,0],[0,0,-(-xi**3)**0.5]],[[-xi,0,0],[0,(-xi**3)**0.5,0],[0,0,-(-xi**3)**0.5]],[[-xi,0,0],[0,(xi**3)**0.5,0],[0,0,-(xi**3)**0.5]],[[xi,0,0],[0,(xi**3)**0.5,0],[0,0,-(xi**3)**0.5]]]
SIrreps3Delta384 = [[[0,1,1],[1/2,1/2,-1/2],[1/2,-1/2,1/2]],[[0,-1,-1],[-1/2,-1/2,1/2],[-1/2,1/2,-1/2]],[[0,-1,-1],[-1/2,-1/2,1/2],[-1/2,1/2,-1/2]],[[0,1,1],[1/2,1/2,-1/2],[1/2,-1/2,1/2]],[[0,-1,-1],[-1/2,-1/2,1/2],[-1/2,1/2,-1/2]],[[0,1,1],[1/2,1/2,-1/2],[1/2,-1/2,1/2]],[[0,-1,-1],[-1/2,1/2,-1/2],[-1/2,-1/2,1/2]],[[0,1,1],[1/2,-1/2,1/2],[1/2,1/2,-1/2]],[[0,-1,-1],[-1/2,1/2,-1/2],[-1/2,-1/2,1/2]],[[0,1,1],[1/2,-1/2,1/2],[1/2,1/2,-1/2]],[[0,-1,-1],[-1/2,1/2,-1/2],[-1/2,-1/2,1/2]],[[0,1,1],[1/2,-1/2,1/2],[1/2,1/2,-1/2]],[[0,-1,-1],[-1/2,1/2,-1/2],[-1/2,-1/2,1/2]],[[0,1,1],[1/2,-1/2,1/2],[1/2,1/2,-1/2]]]
TIrreps4Delta384 = []
SIrreps4Delta384 = []
TIrreps5Delta384 = []
SIrreps5Delta384 = []
TIrreps6Delta384 = [[[1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,xi**7,0,0,0],[0,0,0,xi**3,0,0],[0,0,0,0,xi,0],[0,0,0,0,0,xi**5]],[[1j,0,0,0,0,0],[0,-1j,0,0,0,0],[0,0,xi**(-1/2),0,0,0],[0,0,0,-xi**(-1/2),0,0],[0,0,0,0,xi**(-3/2),0],[0,0,0,0,0,-xi**(-3/2)]],[[1j,0,0,0,0,0],[0,-1j,0,0,0,0],[0,0,xi**(3/2),0,0,0],[0,0,0,-xi**(3/2),0,0],[0,0,0,0,xi**(1/2),0],[0,0,0,0,0,-xi**(1/2)]],[[1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,xi**(3/2),0,0,0],[0,0,0,-xi**(3/2),0,0],[0,0,0,0,xi**(-3/2),0],[0,0,0,0,0,-xi**(-3/2)]],[[1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,xi**(-1/2),0,0,0],[0,0,0,-xi**(-1/2),0,0],[0,0,0,0,xi**(1/2),0],[0,0,0,0,0,-xi**(1/2)]],[[xi,0,0,0,0,0],[0,-xi,0,0,0,0],[0,0,xi**(1/2),0,0,0],[0,0,0,-xi**(1/2),0,0],[0,0,0,0,xi**(-3/2),0],[0,0,0,0,0,-xi**(-3/2)]],[[xi**(-1),0,0,0,0,0],[0,-xi**(-1),0,0,0,0],[0,0,xi**(3/2),0,0,0],[0,0,0,-xi**(3/2),0,0],[0,0,0,0,xi**(-1/2),0],[0,0,0,0,0,-xi**(-1/2)]]]
SIrreps6Delta384 = [[[0,0,-0.3535533906+0.3535533906j,0.3535533906-0.3535533906j,-0.5,-0.5],[0,0,0.3535533906+0.3535533906j,-0.3535533906-0.3535533906j,-0.5j,-0.5j],[-0.3535533906-0.3535533906j,0.3535533906-0.3535533906j,0,0,0.3535533906+0.3535533906j,-0.3535533906-0.3535533906j],[0.3535533906+0.3535533906j,-0.3535533906+0.3535533906j,0,0,0.3535533906+0.3535533906j,-0.3535533906-0.3535533906j],[-0.5,0.5j,0.3535533906-0.3535533906j,0.3535533906-0.3535533906j,0,0],[-0.5,0.5j,-0.3535533906+0.3535533906j,-0.3535533906+0.3535533906j,0,0]],[[0,0,-0.5j,-0.4619397663-0.1913417162j,-0.1913417162+0.4619397663j,0.1913417162-0.4619397663j],[0,0,-0.5j,-0.4619397663-0.1913417162j,0.1913417162-0.4619397663j,-0.1913417162+0.4619397663j],[0.5j,0.5j,0,0,-0.4619397663-0.1913417162j,-0.4619397663-0.1913417162j],[-0.4619397663+0.1913417162j,-0.4619397663+0.1913417162j,0,0,0.5j,0.5j],[-0.1913417162-0.4619397663j,0.1913417162+0.4619397663j,-0.4619397663+0.1913417162j,-0.5j,0,0],[0.1913417162+0.4619397663j,-0.1913417162-0.4619397663j,-0.4619397663+0.1913417162j,-0.5j,0,0]],[[0,0,-0.5j,-0.1913417162+0.4619397663j,-0.4619397663-0.1913417162j,0.3535533906-0.3535533906j],[0,0,-0.5j,-0.1913417162+0.4619397663j,0.4619397663+0.1913417162j,-0.3535533906+0.3535533906j],[0.5j,0.5j,0,0,-0.1913417162+0.4619397663j,0.3535533906+0.3535533906j],[-0.1913417162-0.4619397663j,-0.1913417162-0.4619397663j,0,0,0.5j,0.4619397663+0.1913417162j],[-0.4619397663+0.1913417162j,0.4619397663-0.1913417162j,-0.1913417162-0.4619397663j,-0.5j,0,0],[0.3535533906+0.3535533906j,-0.3535533906-0.3535533906j,0.3535533906-0.3535533906j,0.4619397663-0.1913417162j,0,0]],[[0,0,0.3535533906-0.3535533906j,0.4619397663-0.1913417162j,-0.1913417162-0.4619397663j,0.1913417162+0.4619397663j],[0,0,-0.3535533906+0.3535533906j,-0.4619397663+0.1913417162j,-0.1913417162-0.4619397663j,0.1913417162+0.4619397663j],[0.3535533906+0.3535533906j,-0.3535533906-0.3535533906j,0,0,0.1913417162-0.4619397663j,0.1913417162-0.4619397663j],[0.4619397663+0.1913417162j,-0.4619397663-0.1913417162j,0,0,0.5j,0.5j],[-0.1913417162+0.4619397663j,-0.1913417162+0.4619397663j,0.1913417162+0.4619397663j,-0.5j,0,0],[0.1913417162-0.4619397663j,0.1913417162-0.4619397663j,0.1913417162+0.4619397663j,-0.5j,0,0]],[[0,0,-0.3535533906+0.3535533906j,-0.1913417162-0.4619397663j,-0.4619397663+0.1913417162j,0.5j],[0,0,0.3535533906-0.3535533906j,0.1913417162+0.4619397663j,-0.4619397663+0.1913417162j,0.5j],[-0.3535533906-0.3535533906j,0.3535533906+0.3535533906j,0,0,0.4619397663+0.1913417162j,-0.3535533906+0.3535533906j],[-0.1913417162+0.4619397663j,0.1913417162-0.4619397663j,0,0,0.5j,-0.4619397663-0.1913417162j],[-0.4619397663-0.1913417162j,-0.4619397663-0.1913417162j,0.4619397663-0.1913417162j,-0.5j,0,0],[-0.5j,-0.5j,-0.3535533906-0.3535533906j,-0.4619397663+0.1913417162j,0,0]],[[0,0,0.5j,-0.1913417162-0.4619397663j,0.4619397663-0.1913417162j,0.5j],[0,0,0.5j,-0.1913417162-0.4619397663j,-0.4619397663+0.1913417162j,-0.5j],[-0.5j,-0.5j,0,0,-0.1913417162-0.4619397663j,-0.5],[-0.1913417162+0.4619397663j,-0.1913417162+0.4619397663j,0,0,-0.5j,-0.4619397663-0.1913417162j],[0.4619397663+0.1913417162j,-0.4619397663-0.1913417162j,-0.1913417162+0.4619397663j,0.5j,0,0],[-0.5j,0.5j,-0.5,-0.4619397663+0.1913417162j,0,0]],[[0,0,-0.5j,0.4619397663+0.1913417162j,0.5j,0.5j],[0,0,-0.5j,0.4619397663+0.1913417162j,-0.5j,-0.5j],[0.5j,0.5j,0,0,-0.5,0.5],[0.4619397663-0.1913417162j,0.4619397663-0.1913417162j,0,0,-0.1913417162-0.4619397663j,0.1913417162+0.4619397663j],[-0.5j,0.5j,-0.5,-0.1913417162+0.4619397663j,0,0],[-0.5j,0.5j,0.5,0.1913417162-0.4619397663j,0,0]]]

TIrrepsDelta384 = [TIrreps1Delta384, TIrreps2Delta384, TIrreps3Delta384, TIrreps4Delta384, TIrreps5Delta384, TIrreps6Delta384]
SIrrepsDelta384 = [SIrreps1Delta384, SIrreps2Delta384, SIrreps3Delta384, SIrreps4Delta384, SIrreps5Delta384, SIrreps6Delta384]

ConjugacyClassesSizesDelta384 = [1,24,128,6,6,3,24,24,24,6,6,6,3,6,3,24,24,24,24,6,3,3,3,3]

CharacterTable1Delta96 = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,-1,1,1,1,1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,1]]
CharacterTable2Delta96 = [[2,0,-1,2,2,2,0,0,0,2,2,2,2,2,2,0,0,0,0,2,2,2,2,2]]
CharacterTable3Delta96 = [[3,1,0,-1,3,3,-1,1,1,-1,-1,-1,-1,-1,3,-1,-1,-1,1,-1,-1,-1,3,-1],[3,-1,0,-1,3,3,1,-1,-1,-1,-1,-1,-1,-1,3,1,1,1,-1,-1,-1,-1,3,-1],[3,-1,0,1,-1,3,1j,1,-1,1,-1-2j,1,-1+2j,1,-1,-1j,-1j,1j,1,-1+2j,-1-2j,-1+2j,-1,-1-2j],[3,1,0,1,-1,3,-1j,-1,1,1,-1-2j,1,-1+2j,1,-1,1j,1j,-1j,-1,-1+2j,-1-2j,-1+2j,-1,-1-2j],[3,-1,0,1,-1,3,-1j,1,-1,1,-1+2j,1,-1-2j,1,-1,1j,1j,-1j,1,-1-2j,-1+2j,-1-2j,-1,-1+2j],[3,1,0,1,-1,3,1j,-1,1,1,-1+2j,1,-1-2j,1,-1,-1j,-1j,1j,-1,-1-2j,-1+2j,-1-2j,-1,-1+2j],[3,1,0,1-xi+xi**3,1,-1,xi**3,-1j,-1,1+xi-xi**3,-1j,-1+xi+xi**3,-2*xi-xi**2,-1-xi-xi**3,-1-2j,-xi,xi,-xi**3,1j,1j,xi**2+2*xi**3,2*xi-xi**2,-1+2j,xi**2-2*xi**3],[3,-1,0,1-xi+xi**3,1,-1,-xi**3,1j,1,1+xi-xi**3,-1j,-1+xi+xi**3,-2*xi-xi**2,-1-xi-xi**3,-1-2j,xi,-xi,xi**3,-1j,1j,xi**2+2*xi**3,2*xi-xi**2,-1+2j,xi**2-2*xi**3],[3,1,0,1+xi-xi**3,1,-1,-xi**3,-1j,-1,1-xi+xi**3,-1j,-1-xi-xi**3,2*xi-xi**2,-1+xi+xi**3,-1-2j,xi,-xi,xi**3,1j,1j,xi**2-2*xi**3,-2*xi-xi**2,-1+2j,xi**2+2*xi**3],[3,-1,0,1+xi-xi**3,1,-1,xi**3,1j,1,1-xi+xi**3,-1j,-1-xi-xi**3,2*xi-xi**2,-1+xi+xi**3,-1-2j,-xi,xi,-xi**3,-1j,1j,xi**2-2*xi**3,-2*xi-xi**2,-1+2j,xi**2+2*xi**3],[3,1,0,1+xi-xi**3,1,-1,xi,1j,-1,1-xi+xi**3,1j,-1+xi+xi**3,xi**2-2*xi**3,-1-xi-xi**3,-1+2j,-xi**3,xi**3,-xi,-1j,-1j,2*xi-xi**2,xi**2+2*xi**3,-1-2j,-2*xi-xi**2],[3,-1,0,1+xi-xi**3,1,-1,-xi,-1j,1,1-xi+xi**3,1j,-1+xi+xi**3,xi**2-2*xi**3,-1-xi-xi**3,-1+2j,xi**3,-xi**3,xi,1j,-1j,2*xi-xi**2,xi**2+2*xi**3,-1-2j,-2*xi-xi**2],[3,1,0,1-xi+xi**3,1,-1,-xi,1j,-1,1+xi-xi**3,1j,-1-xi-xi**3,xi**2+2*xi**3,-1+xi+xi**3,-1+2j,xi**3,-xi**3,xi,-1j,-1j,-2*xi-xi**2,xi**2-2*xi**3,-1-2j,2*xi-xi**2],[3,-1,0,1-xi+xi**3,1,-1,xi,-1j,1,1+xi-xi**3,1j,-1-xi-xi**3,xi**2+2*xi**3,-1+xi+xi**3,-1+2j,-xi**3,xi**3,-xi,1j,-1j,-2*xi-xi**2,xi**2-2*xi**3,-1-2j,2*xi-xi**2]]
CharacterTable4Delta96 = []
CharacterTable5Delta96 = []
CharacterTable6Delta96 = [[6,0,0,-2,-2,6,0,0,0,-2,2,-2,2,-2,-2,0,0,0,0,2,2,2,-2,2],[6,0,0,0,-2,-2,0,0,0,0,2,-2*xi-2*xi**3,-2-2*xi-2*xi**3,2*xi+2*xi**3,2,0,0,0,0,2,-2+2*xi+2*xi**3,-2+2*xi+2*xi**3,2,-2-2*xi-2*xi**3],[6,0,0,0,-2,-2,0,0,0,0,2,2*xi+2*xi**3,-2+2*xi+2*xi**3,-2*xi-2*xi**3,2,0,0,0,0,2,-2-2*xi-2*xi**3,-2-2*xi-2*xi**3,2,-2+2*xi+2*xi**3],[6,0,0,2*xi-2*xi**3,-2,-2,0,0,0,-2*xi+2*xi**3,-2,0,2-2*xi+2*xi**3,0,2,0,0,0,0,-2,2-2*xi+2*xi**3,2+2*xi-2*xi**3,2,2+2*xi-2*xi**3],[6,0,0,-2*xi+2*xi**3,-2,-2,0,0,0,2*xi-2*xi**3,-2,0,2+2*xi-2*xi**3,0,2,0,0,0,0,-2,2+2*xi-2*xi**3,2-2*xi+2*xi**3,2,2-2*xi+2*xi**3],[6,0,0,-2,2,-2,0,0,0,-2,2j,2,2j,2,-2-4j,0,0,0,0,-2j,-2j,2j,-2+4j,-2j],[6,0,0,-2,2,-2,0,0,0,-2,-2j,2,-2j,2,-2+4j,0,0,0,0,2j,2j,-2j,-2-4j,2j]]

CharacterTableDelta384 = [CharacterTable1Delta96, CharacterTable2Delta96, CharacterTable3Delta96, CharacterTable4Delta96, CharacterTable5Delta96, CharacterTable6Delta96]

Delta384_Theory = [IrrepsDelta384, SIrrepsDelta384, TIrrepsDelta384, 16, ConjugacyClassesSizesDelta384, CharacterTableDelta384, "[384,568]","Δ(384)"]


#=========================
#FullGroupTheory and Tests
#=========================
FullGroupTheory = [Gamma2_Theory, Gamma3_Theory, DoubleGamma3_Theory, Gamma4_Theory, DoubleGamma4_Theory, Gamma5_Theory, 
                   DoubleGamma5_Theory, Gamma6_Theory, S3Z3_Theory, A4Z2_Theory, A4Z4_Theory, TZ2_Theory, S4Z3_Theory, SemiTZ4_Theory, TZ4_Theory,
                   Zs343_Theory, Zs76_Theory, Zs4223_Theory, Zs336_Theory, 
                   Zs136_Theory, Zs4423_Theory, 
                   BD3_Theory, BO_Theory, 
                   GL2_3_Theory, Delta96_Theory, Delta384_Theory
                   ]

FullGroupTheoryTest(FullGroupTheory)
