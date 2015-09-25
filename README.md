# roulette

This is a MATLAB implementation of parameter sampling methods for Markov Jump Processes (aka Continuous Time Markov Chains), using random truncations of the state-space. The method used to truncate the space is called Russian Roulette.

The experiments folder contains sample files on how to call the sampler, and examples of model descriptions and observation files. More detailed instructions can be found in experiments/predPreyGibbsRT.m.

The expmv folder is taken from http://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector and is copyright (c) 2010, Nick Higham and Awad Al-Mohy. It is the implementation of the method described in:
A. H. Al-Mohy and N. J. Higham, Computing the action of the matrix exponential, with an application to exponential integrators, SIAM Journal on Scientific Computing
