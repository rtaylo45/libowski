# libowski
A numberical solver even the Dude can git behind
![The Dude](https://images2.minutemediacdn.com/image/upload/c_crop,h_1191,w_2118,x_41,y_0/f_auto,q_auto,w_1100/v1554931540/shape/mentalfloss/61708-gramercy_pictures.jpg "The Dude")

This is a matrix expotential library that is being build as part of my PhD work. The purpose is to solve the set of nuclear burnup equations along with species transport. The governing equations are a set of linear first order ODEs as follows.

N'(t) = A*N

A is the transition matrix. It containts the coefficients for burnup as well as the convection source terms 

# Usage
Pull the code to your favorite place and create a `build` folder in the main `libowski` directory. `CD` into the `build` folder and run `cmake ..` to generate the make file. Then run `make all` to build the project. The only dependence is the eigen3 library. It is used for all matrix types and linear algebra functions. It is included in the download.

Testing is done using `ctest`. To run the unit test run `make test` in your build directory.
