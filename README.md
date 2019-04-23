**@Author: `Ali Hashmi`**

credits: `Sham Tlili` for the strain-rate measure. 
###### The associated Mathematica functions for strain-rate measure are translations of her Matlab code

`J.M` (stack exchange) for `getEigenSystem` implementation 
###### I used it rather than the built-in `EigenSystem` to remain true to Sham's code

___

The script `pyramidal-KLT.m` is based on a pyramidal implementaion of Lucas Kanade motion tracking algorithm. It can be used
to obtain flow-fields at a given instance using images. The features to be tracked in the flow are taken as points inside a masked
region. 

`VectorFieldOperations.m` contains functions for computing velocity gradients and strain-rate from velocity-fields

`pyramidal-KLT-allmethods.m` is a merger of the two script files


<img src="https://github.com/alihashmiii/flow-fields/blob/master/for%20Readme/plot.png">
