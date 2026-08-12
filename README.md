# Supplementary R Code for "Marginal Analysis of Clustered Interval-Censored Competing-Risks Data with Informative Cluster Size"

The R implementation of the proposed method, a simulated example dataset, and an example analysis script are provided


Files included:

1. analysis.R :
   
       Main R function implementing the proposed marginal analysis procedure based on the generalized odds-rate transformation model
       using B-spline-based sieve estimation.

2. example_dataset.csv :
   
       A simulated example dataset illustrating the required input data structure.

3. run example_dataset.R :

       Example R script demonstrating how to load the proposed method, read the simulated example dataset, fit the model, and obtain
       regression coefficient estimates.

Before running the example, please make sure that analysis.R, example_dataset.csv, and run example_dataset.R are in the same folder.


Notes

The provided dataset is simulated and is intended only to illustrate the implementation of the proposed method.
It is not intended to reproduce the full Monte Carlo simulation results reported in the manuscript.
