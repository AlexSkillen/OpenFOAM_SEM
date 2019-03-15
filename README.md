# Synthetic eddy method for OpenFOAM

This code provides an inlet conditions for LES, DES, etc. It is based on the
synthetic eddy method, with the full algorithm described in: 

A. Skillen, A. Revell and T. Craft. "Accuracy and Efficiency Improvements in
Synthetic Eddy Methods". International Journal of Heat and Fluid Flow. (2016) 
DOI: 10.1016/j.ijheatfluidflow.2016.09.008

This code has been tested on OpenFOAM v1806. 


## How to ...

### Compile

Clone the repository and compile with `wmake` e.g.:

```bash
git clone https://github.com/AlexSkillen/OpenFOAM_SEM.git
cd OpenFOAM_SEM
wmake
```

By default the compiled library will be in `FOAM_USER_LIBBIN`. 

### Run

Add to your `system/controlDict`

```
libs ("libSEM.so");
```

Finally add to your velocity input conditions the following subdictionary
(assuming inlet is called `inlet`):

```
boundaryField
{
    ...

    inlet
    {
        type    turbulentInletSEM;
        value   uniform (0.0 0.0 0.0);
        UIn     uniform (0.497899  0.0 0.0); 
        RIn     uniform (0.0436187 -0.000158973 0.0 1.11767e-05 0.0 0.0158168);
        sigma   uniform (0.04 0.04 0.04);
    }
}
```

There are three user inputs: `UIn` is the prescribed mean velocity at the
inlet. `RIn` is the prescribed Reynolds stress tensor at the inlet. Sigma is
the prescribed lengthscale of the turbulence (which may be anisotropic). All 
inputs can be inhomogeneous if desired, by setting them to nonuniform Lists.

Sigma can usually be estimated from engineering judgement, or from a prior RANS 
calculation as k^{3/2} / epsilon, where k and epsilon are the turbulent 
kinetic energy, and dissipation rate respectively. For an anisotropic 
length-scale, the siszes can be estimated from an RSM closure, with sigma 
being set to (uu^{3/2}/epsilon, vv^{3/2}/epsilon, ww^{3/2}/epsilon). 
In any case, it is very important to manually clip the lengthscale
returned by any empirical estimate to prevent excessively small or large eddies. 
For example, eddies should not be smaller than the local LES filter width, or
maximum cell dimension. Eddies should also not be larger than the geometry for
internal flows (typically, for a channel flow, sigma < 0.4 delta seems to work
well).

