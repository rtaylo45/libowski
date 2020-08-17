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

The first term represents the change in species with time, the second is the flux through the volume surface the last is the change in concentration with a volumetric source term. The flux through the surface is comprised of a convective flux and diffusive. Volumetric source terms include those found in nuclear reactors such as decay, transmutation and generation from fission. Other source terms come from chemical reactions, both linear and nonlinear. These can include chemical reactions, phase migration and surface reactions. After the control volume method is applied to the integral transport equation,

<p align="center">
<img src="http://www.sciweavers.org/tex2img.php?eq=%20%20%20%20%5Cfrac%7B%5Cpartial%20%5Coverline%7BC%7D_%7Bi%7D%7D%7B%5Cpartial%20t%7D%20%3D%20-%20%5Cfrac%7B1%7D%7BV%7D%5Ciint_S%20n%20%5Ccdot%20F_%7Bi%7D%20%5C%2CdS%20%2B%20%5Csum%20%5Coverline%7BC%7D_%7BV%2Ci%7D%0A&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="    \frac{\partial \overline{C}_{i}}{\partial t} = - \frac{1}{V}\iint_S n \cdot F_{i} \,dS + \sum \overline{C}_{V,i}" width="251" height="51" />
</p>

Next the method of lines is used to discretize the spatial varialbes. For a fixed control volume in 2D cartesian gemoetry, the integral transport equation for species i is shown below.

<p align="center">
<img src="http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7Bd%20%5Coverline%7BC%7D_%7Bi%7D%7D%7Bd%20t%7D%20%3D%20-%20%5Cfrac%7B1%7D%7BV%7D%5Cbigg%5B%20F_%7Bx%7DA_%7Bx%7D%20%2B%20F_%7By%7DA_%7By%7D%5Cbigg%5D%20%2B%20%5Csum%20%5Coverline%7BC%7D_%7BV%2Ci%7D%0A&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\frac{d \overline{C}_{i}}{d t} = - \frac{1}{V}\bigg[ F_{x}A_{x} + F_{y}A_{y}\bigg] + \sum \overline{C}_{V,i}" width="278" height="49" />
</p>

The flux in x and y are approximated using a first order upwind differencing scheme for convection and a second order central differencing scheme for diffusion. 

The governing equations become are a set of first order ODEs as follows.

<p align="center">
<img src="http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7Bd%20C%7D%7Bd%20t%7D%20%3D%20AC%20%2B%20F%28C%2Ct%29%0A&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\frac{d C}{d t} = AC + F(C,t)" width="144" height="43" />
</p>

A is the transition matrix. It containts the coefficients for burnup as well as the convection/diffusion source terms. The nonlinear source terms are collected in F(c,t). This equation is in vector matrix notation. C is a vector of species concentrations. We are solving for the concentrations in each cell volume, so the size of our system is the number of species X number of cells. Solving the system of ODE's is done witht matrix exponential time differencing. This solution is 

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

This formula is exact and ETD methods approximate the integral in the expression. The problem lies in evaluating the exponential of the transition matrix. Currently the only matrix exponential solver in libowski is based transforming the matrix exponential with Cauchy's integral formula into a contour integral that winds the spectrum A. This is done with either CRAM or a ration approximation derived using a contour on the left hand side of the complex plane. This might not be the best, other matrix exponential solvers are currently in the works. 

# Usage
Pull the code to your favorite place and create a `build` folder in the main `libowski` directory. `CD` into the `build` folder and run `cmake ..` to generate the make file. Then run `make all` to build the project. The only dependence is the eigen3 library. It is used for all matrix types and linear algebra functions. It is included in the download.

Testing is done using `ctest`. To run the unit test run `make test` in your build directory.
