# Algorithmic-construction-of-VVMFs


## Overview
The code in this repository aims to compute the Vector-Valued Modular Forms (VVMF) of any finite modular groups. This type of a function appears naturally in modular symmetric theories. The ability to construct them is crucial in order to make predictions. They are characterized by a representation of the finite modular group (a quotient group of SL(2,Z)) and a positive integer called the weight.

The usual method to compute VVMFs it to build the first ones and do successive tensor products of them to get the higher weights. This assumes to be able to compute the first VVMFs and to have a criterion to know if a set of functions can produce the others. The computation of the lowest-weight VVMFs relies on general result on the mathematical of VVMFs while the criterion on a representation-theoric property namely the faithfulness of the representations of these VVMFs.

This code is associated with the following scientific publication:

- **Paper:** [Algorithmic Construction of Vector-Valued Modular Forms and Application to an A4xZ4-invariant Model of Leptons]
- **Authors:** [Quentin Bonnefoy, Michel Rausch de Traubenberg and Victor Saulquin]
- **DOI:** [DOI or link, if available]



## Usage

The main file of the repository is VVMFBuilder.py. If you want to compute the VVMFs of a tabulated group for e.g. A4 modify the first line of the file to specify the group you are interested in. The list of the tabulated groups lies in the file GroupTheoryDictionnary.py. For our example, the line becomes: group_theory=FullGroupTheory[1] since A4 is the second tabulated group in GroupTheoryDictionnary (indexation begins at 0 in Python).

Then, the code needs the Clebsch-Gordan (CG) coefficients of your group. They are located in the folder ClebschGordanLibrary. Change the second line to tell the code in which file of the folder the CG coefficients of the group are and what is the function to use in this file. This becomes: from ClebschGordanLibrary.Gamma3_CG import ProductGamma3 as Product

Finally, tell the code which weight you want. The code can be launched. It will first identify the minimal set thank to MinimalSet.py and describe it is in the consol. In order to build the VMMFs of the minimal set, it will ask you to copy a line L=[r1,r2,...] and paste it in the dedicated Mathematica file. These files can be found in the folder MathematicaFiles. Paste the line at the beginning of the file it will tell you, execute it and copy the list of coefficients at the end in the Python console. 

Then, thanks to MinimalVVMF.py, it will first compute what the overall constants of the entries of the VVMFs are since they are not predicted by the formulae encoded in the Mathematica files. They are printed in the console and file builds the minimal VVMFs. Finally, VVMFBuilder.py will automatically do the tensor products to reach the requested weight, eliminate the linear dependent VVMFs, normalize them and print the result in the console. These are the VVMFs of your group at the weight you asked.

If you now want to change the weight, you just have to modify the corresponding line of VVMFBuilder.py. However, since the minimal set has already been built, the Boolean FirstUse can be settle to zero to skip the first steps and directly perform the tensor products. Pay attention that FirstUse = True when you reopen the Python file or when you change the group under study.

In case the group you want to focus on is not already tabulated, add it to GroupTheoryDictionnary.py to collect the representation data needed. The needed information and the format are described at the beginning of this file. For the CG coefficients, create a file GroupName_CG in the folder ClebschGordanLibrary. Enter the coefficients in this file. It must contain a fonction ProductGroupName whose arguments are two multiplets and the output is a list of multipets corresponding to the decomposition of the product of the two initial ones. A multiplet is a list, the first entry being the name of the irreducible representation and the second the list of its components.

More details are available in the associated paper.
## Citation

If you use this repository in your research, please cite the associated article and its associated scientific publication.

**Associated publication:**

[Full reference will be added upon publication.]

A `CITATION.cff` file is also provided in this repository.

## License

This project is distributed under the MIT license. See the `LICENSE` file for details.

## Contact

For questions or suggestions, please contact:

victor.saulquin@iphc.cnrs.fr
