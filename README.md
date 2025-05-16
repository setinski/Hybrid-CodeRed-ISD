# Hybrid-CodeRed-ISD

In this project, we analyze the concrete complexity of hybrid ISD algorithms, which combine the ISD framework with code reduction techniques, as well as their non-hybrid counterparts used to solve the codeword-finding problem.

## Project Dependencies

- Standard C++ libraries

## Running the Experiments

To run the experiments and obtain both the empirical results and the theoretical predictions of the success probability, proceed as follows:

1. Modify the parameter `NUM_BITS` in `binop.hpp` to suit your needs. It is currently hard-coded to `1280`, which corresponds to the code length.
2. Compile and run `experiment.cpp`.

By default, the code dimension `k` is set to half the code length. The weight parameter `w2` is hard-coded to `3`, and the number of merging levels is set to `2`. The collision parameters of the Stern-based and Wagner-based algorithms are tuned to generate a list at each level in amortized time 1.

To change any of these parameters, modify them in the `main` function of `experiment.cpp`.
