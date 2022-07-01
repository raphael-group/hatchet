# check

This command of HATCHet verifies that dependencies that are unique to specific parts of HATCHet are correctly installed. The relevant commands with dependencies are `count-reads`, `phase-snps`, and `compute-cn`.

All checks can be run simultaneously via `hatchet check`, or an individual command can be checked via, e.g., `hatchet check compute-cn`.

The check for `compute-cn` runs the step on a set of small data files (.bbc/.seg) pre-packaged with HATCHet, and is a quick way to verify if your solver is working correctly.
If you are unable to run this command, it likely indicates a licensing issue with default (Gurobi) solver. To use alternative solvers, see the
[Using a different Pyomo-supported solver](README.html#usingasolver_other) section of the README for more details.

## Input

This command takes no inputs.

## Output

This command produces debugging output from the solver as well as other dependency checks. For `compute-cn`, look for the message towards the end of the output:

```
# Your current solver ... seems to be working correctly
```

to verify if your selected solver works correctly.
