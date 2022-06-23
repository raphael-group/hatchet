# check

This command of HATCHet runs the `compute_cn` step on a set of small data files (.bbc/.seg) pre-packaged with HATCHet, and is a quick way to verify if your solver is working correctly.
If you are unable to run this command, it likely indicates a licensing issue with default (Gurobi) solver. To use alternative solvers, see the
[Using a different Pyomo-supported solver](README.html#usingasolver_other) section of the README for more details).

Other checks are planned in the future that will check if other 3rd party dependencies are working correctly. You will still only need to run the `check` command to run a series of these diagnostic checks.

## Input

This command takes no inputs.

## Output

This command produces debugging output from the solver. Look for the message towards the end of the output:

```
# Your current solver ... seems to be working correctly
```

to verify if your selected solver works correctly.
