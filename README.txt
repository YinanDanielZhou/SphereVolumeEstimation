To compile the program, run 

g++ *.cpp -o proj

then to run the program, run

./proj 2

Here the argument 2 is the number of dimensions.
The monte carlo method support dimension 2 to 7, takes long time on higher dimensions.
The cube based method support dimension 2 and 3, and takes very long time on both. 

To run the program with fixed number of sample points (cubes) N = 10^6, run

./proj 2 fixedN
