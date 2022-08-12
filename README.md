# Subgradient step-size update rules in Lagrangian relaxations of chance-constrained optimization models
This project is based on the paper **Statistical performance of subgradient step-size update rules in Lagrangian relaxations of chance-constrained optimization models** by C. Ritter and B. Singh. 
There, using different computational experiments and statistical tests on two example instances of chance-constrained optimization problems, we investigate whether particular choices of the step-size rules are better suited for solving the Lagrangian relaxation via the subgradient method for these kinds of models in practice.
We present the results of those tests and several additional investigations related to bounding techniques.
We refer to the paper for further information on methods and conclusions.
A preprint of the manuscript that is currently under peer review is available at xxx.
## Repository content
The repository contains the following content:
- `Code` contains all the GAMS code one needs to run the algorithms presented in the paper. Please refer to the subfolders for a specific instance of the desired model. Starter files for the different algorithms are included. `main_for_naive.gms` for the naive solution method, `main_for_lr1.gms` up to `main_for_lr6.gms` for the Lagrangian relaxation algorithm (Alg. 1) using step-size rule one to six, `SAA.gms` as part 1 of Algorithm 2 and `fixed_problem.gms` as part 2 of Algorithm 2. Analogously, `SAA_hi.gms` as part 1 of Algorithm 3 and `fixed_problem_hi.gms` as part 2 of Algorithm 3. We also have a starter file `MainHeuristic.gms` for Algorithm 4. For Model I, that is also included and put out in the Lagrangian relaxation. 
- `Results` contains the Excel tables and figures visualizing the results, some of which are also included in the paper. We differentiate between the `Example_Instances` and `Samples` for each model instance. The `Example_Instances` include the 16 instances we presented separately in the paper for Model I and Model II, respectively. The `Samples` include all 20 batches of each scenario-size regime for Model I and Model II, respectively. For Model II, the example instances can be found in a subfolder of `Model II` in `Code`. For Model I, the example instances are part of the 20 batches of each regime. Which of those we used as the example instances can be read in `Example_instance_Model_I.txt` in the `Model I` subfolder of `Code`.
## Requirements to run code
The code uses the optimization solver Gurobi and the modeling language GAMS which are therefore required to run it.