# Inverse dielectric
**Inverse dielectric** is a modern Fortran library providing a common framework for GW codes, specifically for the most common operations involving the inverse dielectric function.

**Supported operations:**
* _Anisotropic averaging for **q**=**Γ**_: the library contains functions for the computation of the inverse dielectric function at **q**=**Γ** for 3D systems. This last by numerically computing the anisotropic average in a small region (i.e. a small version of the reciprocal cell) centered around  **Γ** using a 131<sup>st</sup>-order Lebedev grid. 
* _Body inversion_: the code provides functions for inverting the body of the dielectric matrix either using CPUs or GPUs. 


