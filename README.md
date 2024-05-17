# A partial Benders decomposition (PBD) method for solving choice-based competitive facility location problems (CBCFLP)

## Reference
The PBD method is presented in Legault, R. & Frejinger, E. *A model-free approach for solving choice-based competitive facility location problems using simulation and submodularity*. arXiv:2203.11329 (2022). See our [arXiv preprint](https://arxiv.org/abs/2203.11329):
```
@article{legault2022model,
  title={A model-free approach for solving choice-based competitive facility location problems using simulation and submodularity},
  author={Legault, Robin and Frejinger, Emma},
  journal={arXiv preprint arXiv:2203.11329},
  year={2022}
}
```

## ARTICLE ABSTRACT
This paper considers facility location problems in which a firm entering a market seeks to open facilities on a subset of candidate locations so as to maximize its expected market share, assuming that customers choose the available alternative that maximizes a random utility function. We introduce a deterministic equivalent reformulation of this stochastic problem as a maximum covering location problem with an exponential number of demand points, each of which is covered by a different set of candidate locations. Estimating the prevalence of these preference profiles through simulation generalizes a sample average approximation method from the literature and results in a maximum covering location problem of manageable size. To solve it, we develop a partial Benders reformulation in which the contribution to the objective of the least influential preference profiles is aggregated and bounded by submodular cuts. This set of profiles is selected by a knee detection method that seeks to identify the best tradeoff between the fraction of the demand that is retained in the master problem and the size of the model. We develop a theoretical analysis of our approach and show that the solution quality it provides for the original stochastic problem, its computational performance, and the automatic profile-retention strategy it exploits are directly connected to the entropy of the preference profiles in the population. Computational experiments on existing and new benchmark sets indicate that our approach dominates the classical sample average approximation method on large instances of the competitive facility location problem, can outperform the best heuristic method from the literature under the multinomial logit model, and achieves state-of-the-art results under the mixed multinomial logit model. We characterize a broader class of problems, which includes assortment optimization, to which the solving methodology and the analyses developed in this paper can be extended.

## Directory Structure
The detailed results of our experiments are provided in folder `results`. These results can be reproduced by running scripts `SB/scripts/run_experiments_MNL.jl` and `SB/scripts/run_experiments_MMNL.jl`. The MNL instances are provided in folder `instances`. The MMNL instances have not been uploaded to the repository due to their large size, but they can be generated locally by running script `SB/scripts/generate_instances_MMNL.jl`. The other scripts of folder `scripts` can be executed to reproduce the illustrative examples, the figures, and the tables presented in the paper. Our algorithm is implemented in the `PBD` function in file `SB/src/PBD_method.jl`.
