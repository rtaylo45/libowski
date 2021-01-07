# libowski
![The Dude](https://images2.minutemediacdn.com/image/upload/c_crop,h_1191,w_2118,x_41,y_0/f_auto,q_auto,w_1100/v1554931540/shape/mentalfloss/61708-gramercy_pictures.jpg "The Dude")

This is a c++ library used to solve PDE's using exponential time differencing that is being build as part of my PhD work. The purpose is to solve the set of nuclear burnup equations along with mass transport. The species transport equation is shown below.

<p align="center">
<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Bequation%2A%7D%0A++++%5Cfrac%7B%5Cpartial+%5Crho_%7Bi%7D%7D%7B%5Cpartial+t%7D+%2B+%5Cnabla+%5Ccdot%28%5Crho_%7Bi%7D%28r%2Ct%29v%29+%2B+%5Cnabla+%5Ccdot+j_%7Bi%7D%28r%2Ct%29+%3D+R_%7Bi%7D%28r%2Ct%29%0A%5Cend%7Bequation%2A%7D%0A" 
alt="\begin{equation*}
    \frac{\partial \rho_{i}}{\partial t} + \nabla \cdot(\rho_{i}(r,t)v) + \nabla \cdot j_{i}(r,t) = R_{i}(r,t)
\end{equation*}
">
</p>

The first term represents the change in species with time, the second is the flux through the volume surface the last is the change in concentration with a volumetric source term. For molten salt reactors this equation becomes:

<p align="center">
<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A%5Cbegin%7Bequation%2A%7D%0A%5Cbegin%7Bsplit%7D%0A++++%5Cfrac%7B%5Cpartial+%5Crho_%7Bi%7D%7D%7B%5Cpartial+t%7D%0A++++%2B+%5Cnabla+%5Ccdot+%5Crho_%7Bi%7D%28r%2Ct%29%5Cboldsymbol%7Bv%7D%0A++++%2B+%5Cnabla+%5Ccdot+j_%7Bi%7D%28r%29%0A++++%26%3D%0A++++%5Csum_%7Bj%3D1%7D%5E%7BN%7D%5Cfrac%7BM_%7Bi%7D%7D%7BM_%7Bj%7D%7D%5Cbigg%28b_%7Bj%5Crightarrow+i%7D%5Clambda_%7Bj%7D+%2B+%0A++++%5Csum_%7Bk%3D1%7D%5E%7BK%7D%5Cgamma_%7Bj%5Crightarrow+i%2Ck%7D%5Csigma_%7Bk%2Cj%7D%28r%29%5Cphi%28r%2Ct%29+%5Cbigg%29%5Crho_%7Bj%7D%28r%2Ct%29%5C%5C%0A++++%26-+%5Cbigg%28%5Clambda_%7Bi%7D+%2B+%5Cphi%28r%2Ct%29%5Csum_%7Bk%3D1%7D%5E%7BK%7D+%5Csigma_%7Bk%2Ci%7D%28r%29%5Cbigg%29%5Crho_%7Bi%7D%28r%2Ct%29.%0A%5Cend%7Bsplit%7D%0A%5Cend%7Bequation%2A%7D%0A%5Cend%7Balign%2A%7D%0A" 
alt="\begin{align*}
\begin{equation*}
\begin{split}
    \frac{\partial \rho_{i}}{\partial t}
    + \nabla \cdot \rho_{i}(r,t)\boldsymbol{v}
    + \nabla \cdot j_{i}(r)
    &=
    \sum_{j=1}^{N}\frac{M_{i}}{M_{j}}\bigg(b_{j\rightarrow i}\lambda_{j} + 
    \sum_{k=1}^{K}\gamma_{j\rightarrow i,k}\sigma_{k,j}(r)\phi(r,t) \bigg)\rho_{j}(r,t)\\
    &- \bigg(\lambda_{i} + \phi(r,t)\sum_{k=1}^{K} \sigma_{k,i}(r)\bigg)\rho_{i}(r,t).
\end{split}
\end{equation*}
\end{align*}
">
</p>

where N is the number of nuclides in the system and K is the number of neutron induced reactions for a specific isotope.  The first term on the left hand side represent generation from decay of nuclide j with the second being generation from neutron induced reactions. The third and fourth term includes losses from decay and transmutation reactions. Other source terms can be included from chemical reactions, both linear and nonlinear. These can include chemical reactions, phase migration and surface reactions. After the control volume method is applied to the integral transport equation,

<p align="center">
<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Bequation%2A%7D%0A%5Cbegin%7Bsplit%7D%0A++++%5Cfrac%7B%5Cpartial+%5Coverline%7B%5Crho%7D_%7Bi%7D%7D%7B%5Cpartial+t%7D%0A++++%3D+%5Cfrac%7B-1%7D%7BV%7D%5Cint_%7BV%7D%5Cbigg%28%5Cnabla+%5Ccdot+%5Crho_%7Bi%7D%28r%2Ct%29%5Cboldsymbol%7Bv%7D%0A++++%2B+%26%5Cnabla+%5Ccdot+j_%7Bi%7D%28r%2Ct%29%5Cbigg%29dV%0A++++%2B%0A++++%5Csum_%7Bj%3D1%7D%5E%7BN%7D%5Cfrac%7BM_%7Bi%7D%7D%7BM_%7Bj%7D%7D%5Cbigg%28b_%7Bj%5Crightarrow+i%7D%5Clambda_%7Bj%7D+%2B+%0A++++%5Csum_%7Bk%3D1%7D%5E%7BK%7D%5Cgamma_%7Bj%5Crightarrow+i%2Ck%7D%5Coverline%7B%5Csigma%7D_%7Bk%2Cj%7D%5Coverline%7B%5Cphi%7D+%5Cbigg%29%5Coverline%7B%5Crho%7D_%7Bj%7D%28t%29%5C%5C%0A++++%26-+%5Cbigg%28%5Clambda_%7Bi%7D+%2B+%5Coverline%7B%5Cphi%7D%5Csum_%7Bk%3D1%7D%5E%7BK%7D+%5Coverline%7B%5Csigma%7D_%7Bk%2Ci%7D%5Cbigg%29%5Coverline%7B%5Crho%7D_%7Bi%7D%28t%29.%0A%5Cend%7Bsplit%7D%0A%5Cend%7Bequation%2A%7D%0A" 
alt="\begin{equation*}
\begin{split}
    \frac{\partial \overline{\rho}_{i}}{\partial t}
    = \frac{-1}{V}\int_{V}\bigg(\nabla \cdot \rho_{i}(r,t)\boldsymbol{v}
    + &\nabla \cdot j_{i}(r,t)\bigg)dV
    +
    \sum_{j=1}^{N}\frac{M_{i}}{M_{j}}\bigg(b_{j\rightarrow i}\lambda_{j} + 
    \sum_{k=1}^{K}\gamma_{j\rightarrow i,k}\overline{\sigma}_{k,j}\overline{\phi} \bigg)\overline{\rho}_{j}(t)\\
    &- \bigg(\lambda_{i} + \overline{\phi}\sum_{k=1}^{K} \overline{\sigma}_{k,i}\bigg)\overline{\rho}_{i}(t).
\end{split}
\end{equation*}
">
</p>

The transport flux in x and y are approximated using a second order upwind TVD scheme for convection and a second order central differencing scheme for diffusion. 

The governing equations become are a set of first order ODEs as follows.

<p align="center">
<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Bequation%2A%7D%0A++++%5Cfrac%7Bd%5Cboldsymbol%7B%5Crho%7D%7D%7Bdt%7D+%3D+%5Cboldsymbol%7BL%5Crho%7D+%2B+%5Cboldsymbol%7BN%7D%28t%2C%5Cboldsymbol%7B%5Crho%7D%29.%0A%5Cend%7Bequation%2A%7D%0A" 
alt="\begin{equation*}
    \frac{d\boldsymbol{\rho}}{dt} = \boldsymbol{L\rho} + \boldsymbol{N}(t,\boldsymbol{\rho}).
\end{equation*}
">
</p>

L is the transition matrix. It containts the coefficients for burnup as well as the convection/diffusion source terms. The nonlinear source terms are collected in N. This equation is in vector matrix notation. Rho is a vector of species concentrations. We are solving for the concentrations in each cell volume, so the size of our system is the number of species X number of cells. Solving the system of ODE's is done witht matrix exponential time differencing. This solution is 

<p align="center">
<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cbegin%7Bequation%2A%7D%0A++++%5Cboldsymbol%7B%5Crho%7D%28t_%7Bn%7D+%2B+%5CDelta+t%29+%3D+e%5E%7B%5CDelta+t%5Cboldsymbol%7BL%7D%7D%5Cboldsymbol%7B%5Crho%7D%28t_%7Bn%7D%29+%2B+e%5E%7B%5CDelta+t%5Cboldsymbol%7BL%7D%7D%5Cint_%7B0%7D%5E%7B%5CDelta+t%7De%5E%7B-%5Cboldsymbol%7BL%7D%5Ctau%7D%5Cboldsymbol%7BN%7D%28t_%7Bn%7D+%2B+%5Ctau%2C%5Cboldsymbol%7B%5Crho%7D%28t_%7Bn%7D+%2B+%5Ctau%29%29d%5Ctau.%0A%5Cend%7Bequation%2A%7D%0A" 
alt="\begin{equation*}
    \boldsymbol{\rho}(t_{n} + \Delta t) = e^{\Delta t\boldsymbol{L}}\boldsymbol{\rho}(t_{n}) + e^{\Delta t\boldsymbol{L}}\int_{0}^{\Delta t}e^{-\boldsymbol{L}\tau}\boldsymbol{N}(t_{n} + \tau,\boldsymbol{\rho}(t_{n} + \tau))d\tau.
\end{equation*}
">
</p>


This formula is exact and ETD methods approximate the integral in the expression. The problem lies in evaluating the exponential of the transition matrix. Currently there are 6 algorithms for solving the matrix exponential. One based on a truncated Taylor series (`Taylor`), two on the Pade approximation (`pade-method1` and `pade-method2`) and three on Cauchys integral formula (`CRAM`, `parabolic` and `hyperbolic`). The Krylov substace method is also avaliable as a preprocessing step.

# Dependencies
The eigen3 library was used extensivly in the creation of libowski and is included in the download. Additionally a number of different error functions were required for some of the initial conditions in the unit test. Some of these fucntions are not defined in the standard C++ library, therefor the Faddevva Package is included in `lib/utils` folder in the main directory. 

# Usage
Pull the code to your favorite place and create a `build` folder in the main `libowski` directory. `CD` into the `build` folder and run `cmake ..` to generate the make file. Then run `make all` to build the project. The only dependence is the eigen3 library. It is used for all matrix types and linear algebra functions. It is included in the download.

Testing is done using `ctest`. To run the unit test run `make test` in your build directory. Some unit test require files that have not been uploaded into the reop and so some test will fail. I need to fix this going forward.

This probject can be built with MPI if cmake detects it with the `Find MPI` command. Only the Cauchy algorithms are built to run in parallel and can run on up to 8 procs for CRAM and 16 for the parabolic and hyperbolic solvers. 
