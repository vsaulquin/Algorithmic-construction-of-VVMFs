# Algorithmic-construction-of-VVMFs


## Overview
The code in this repository aims to compute the Vector-Valued Modular Forms (VVMFs) of any finite modular group (i.e., any finite quotient group of SL(2,Z)). This type of functions appears naturally in supersymmetric modular-invariant field theories. The ability to construct them is crucial in order to work out the actual predictions of a given model. They are characterized by a representation of the finite modular group and a positive integer called the weight.
 
The usual method to compute VVMFs it to build the first, low-weight ones, then to perform their successive tensor products to get the higher-weight VVMFs. This assumes that one is able to compute the first VVMFs, and to check whether a given set of VVMFs produces all other by tensor products. In the present code, the computation of the lowest-weight VVMFs relies on general results of the mathematical theory of VVMFs, while the criterion corresponds to a property of the representations featured in the set of the VVMFs, namely their faithfulness.

This code is associated with the following scientific publication:

- **Paper:** [Algorithmic Construction of Vector-Valued Modular Forms and Application to an A4xZ4-invariant Model of Leptons]
- **Authors:** [Quentin Bonnefoy, Michel Rausch de Traubenberg and Victor Saulquin]
- **DOI:** [DOI or link, if available]



## Usage

The main file of the repository is VVMFBuilder.py. If you want to compute the VVMFs of a tabulated group, e.g. A4, modify the first line of the file to specify which group you are interested in. The list of the tabulated groups can be found in the file GroupTheoryDictionnary.py. For our example, the line becomes: group_theory=FullGroupTheory[1] since A4 is the second tabulated group in GroupTheoryDictionnary (indexation begins at 0 in Python).
 
Then, the code needs the Clebsch-Gordan (CG) coefficients of your group. They are located in the folder ClebschGordanLibrary. Change the second line of VVMFBuilder.py so that it displays which file of that folder contains the CG coefficients of the group, and what is the master function to use in this file. For our A4 example, this reads: from ClebschGordanLibrary.Gamma3_CG import ProductGamma3 as Product. As explained below, only the structure of the arguments and of the output of this master function are rigid, the precise details of the implementation of the tensor products are free.
 
Finally, indicate in the third line for which weight you want the code to build the VVMFs.
 
The code can then be launched. It will first identify the minimal set thanks to MinimalSet.py and describe it in the console. In order to build the VMMFs of the minimal set, you will be asked to copy a line L=[r1,r2,...] and paste it in a specific Mathematica file. These files can be found in the folder MathematicaFiles. Paste the line at the beginning of the required file, execute it and copy the final list of coefficients at the end in the Python console.
 
Then, thanks to MinimalVVMF.py, the code will first fix overall constants in the entries of the VVMFs, that are left unspecified in the formulae encoded in the Mathematica files, and print them in the console. VVMFBuilder.py then builds the minimal set VVMFs, automatically perform their tensor products up to the requested weight, eliminate the linearly-dependent VVMFs, normalize the resulting ones and print the result in the console. These are the q-expansions of the VVMFs for the group at the weight that you requested.
 
If you now want to obtain results for another weight, you only have to modify the corresponding line of VVMFBuilder.py. However, since the minimal set has already been built, the Boolean FirstUse can be set to zero in order to skip the first steps and directly perform the tensor products. Pay attention to set back FirstUse = True when you reopen the Python file or when you change the group under study.
 
In case the group you want to focus on is not already tabulated, add all needed representation-theoretic data to GroupTheoryDictionnary.py. The necessary information and the syntax are described at the beginning of this file. For the CG coefficients, create a file GroupName_CG in the folder ClebschGordanLibrary. Enter the coefficients in this file. It must contain a fonction ProductGroupName whose arguments are two multiplets and the output is a list of multiplets corresponding to the decomposition of the product of the two initial ones. A multiplet is a list, the first entry being the name of the irreducible representation and the second the list of its components.

More details can be found in the associated paper.
## Citation

If you use this repository in your research, please cite the associated article.

**Associated publication:**

[Full reference will be added upon publication.]

A `CITATION.cff` file is also provided in this repository.

## License

This project is distributed under the MIT license. See the `LICENSE` file for details.

## Contact

For questions or suggestions, please contact:

victor.saulquin@iphc.cnrs.fr
