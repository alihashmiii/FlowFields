**@Author: `Ali Hashmi`**
___
credits: `Sham Tlili` for the strain-rate measure. 
###### The associated Mathematica functions for strain-rate measure are translations of her Matlab code

`J.M` (Stack Exchange) for the `getEigenSystem` implementation 
##### Please note that you can also use the built-in `EigenSystem` function.
___

The script `pyramidal-KLT.m` is based on a pyramidal implementaion of Lucas Kanade motion tracking algorithm. It can be used
to obtain flow-field at a given instance using consecutive images. The features to be tracked in the flow are taken as points inside a masked region. 

`VectorFieldOperations.m` contains functions for computing velocity gradients and strain-rate from velocity-fields

`pyramidal-KLT-allmethods.m` is a merger of the two script files

![alt-text](https://github.com/alihashmiii/FlowFields/blob/master/for%20Readme/flow_upload_image_github.png)
