# Censored Shifted Wald model

This repository contains files related to my Stan implementation of a hierarchical Bayesian extension of the censored shifted Wald model (Miller et al., 2018).

List of files:
- swald.stan: this is the Stan script that implements the censored shifted Wald.   
- sce.R: R script that shows how to extract censored shifted Wald estimates in the context of the size-congruity effect.
- estimates.csv: a set of parameter estimates from sce.R. Gives censored SW parameters for each participant derived from the distribution of RTs in both congruent (1) and incongruent (2) trials.
- estimates.jasp: a JASP file showing one way to use the censored shifted Wald estimates in inference.
- talk.pdf: my Virtual MathPsych 2023 talk.
