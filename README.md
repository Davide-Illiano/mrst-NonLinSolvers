# mrst-NonLinSolvers
## Description
In this module we study two different linearization schemes: the Newton method, quadratically but 
only locally convergent, and the L-scheme, linearly and globally convergent. The former is already available in MRST,
the latter is implemented in the mrst-NonLinSolvers module given here. A particular case of fully 
coupled multiphase flow and reactive transport is here investigated. We will observe as, due to the complexity of 
the problems, the Newton method fails to converge. The L-scheme, converging also for these cases, 
is here improved thanks to the Anderson acceleration. 
This technique can drastically reduce the numbers of iterations required by the linearly convergent scheme.

In this module we find two separate folders, in one we investigate the Richards equation, in the other we consider also the transportation of
an external component dissolved into the water phase. 
For both of the problems we investigate different solving techiniques. The equations are linearized thanks to the Newton method,
using the authomatic differentiation techiniques already implemented in mrst. We use also an alternative fix point solver, based on the L-scheme. 
This linearization scheme is linearly but globally converegnt and more robust than the Newton method. We observed as the latter fails to converge
in particularly complex domains and configurations. 
Unfortunatelly, the L-shcme, being only linearly convergent, results slowere, in term of numbers of iterations, than the Newton method. 
Thus, we implemented the Anderson acceleration tool, which can drastically improve the convergence rate of linearly convergent schemes.
The obtained improved L-scheme required a number of iterations comparable with the one needed by the Newton method.

We included 5 differente example:
* Richards equation on unsaturated domain
* Richards equation on variably saturated domain
* Richards and transport equations, fully coupled together, on unsaturated domain
* 2 examples regarding Richards and transport equations, fully coupled together, on variably saturated domain

The examples here investigated regards a particular case of multiphase flow in variably saturated porous media. Anyhow, we beleive the 
codes here included, can help in the implementation of the L-scheme also for different sets of equations.

This module was largely based on:

Davide Illiano, Iuliu Sorin Pop, Florin Adrian Radu: Iterative schemes for surfactant transport in porous media, 
arXiv:1906.00224, 2019.

### Prerequisites

* MRST (Tested version: 2019a)
* MATLAB (Tested version: R2019a)


