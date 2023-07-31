# The simulation-based hybrid submodular (SHS) method for solving choice-based competitive facility location problems (CBCFLP)

## Reference
The SHS method is described in Legault, R. & Frejinger, E. A model-free approach for solving choice-based competitive facility location problems using simulation and submodularity. arXiv:2203.11329 (2022). See our [arXiv preprint](https://arxiv.org/abs/2203.11329):
```
@article{legault2022simulation,
  title={A Simulation Approach for Competitive Facility Location with Random Utility Maximizing Customers},
  author={Legault, Robin and Frejinger, Emma},
  journal={arXiv preprint arXiv:2203.11329},
  year={2022}
}
```

## ARTICLE ABSTRACT
This paper considers facility location problems in which a firm entering a market seeks to open a set of available locations so as to maximize its expected market share, assuming that customers choose the alternative that maximizes a random utility function. We introduce a novel deterministic equivalent reformulation of this probabilistic model and, extending the results of previous studies, show that its objective function is submodular under any random utility maximization model. This reformulation characterizes the demand based on a finite set of preference profiles. Estimating their prevalence through simulation generalizes a sample average approximation method from the literature and results in a maximum covering problem for which we develop a new branch-and-cut algorithm. The proposed method takes advantage of the submodularity of the objective value to replace the least influential preference profiles by an auxiliary variable that is bounded by submodular cuts. This set of profiles is selected by a knee detection method. We provide a theoretical analysis of our approach and show that its computational performance, the solution quality it provides, and the efficiency of the knee detection method it exploits are directly connected to the entropy of the preference profiles in the population. Computational experiments on existing and new instances indicate that our approach dominates the classical sample average approximation method, can outperform the best heuristic method from the literature under the multinomial logit model, and achieves state-of-the-art results under the mixed multinomial logit model.

## Directory Structure
The detailed results of our experiments are provided in folder `results`. These results can be reproduced by running scripts `SB/scripts/run_experiments_MNL.jl` and `SB/scripts/run_experiments_MMNL.jl`. The MNL instances are provided in folder `instances`. The MMNL instances have not been uploaded to the repository due to their large size, but they can be generated locally by running script `SB/scripts/generate_instances_MMNL.jl`. The other scripts of folder `scripts` can be executed to reproduce the illustrative examples, the figures, and the tables presented in the paper. Our algorithm is implemented in the `SHS` function in file `SB/source/SHS_method.jl`.
