# SNEK
Quality of Life Functions for CobraPy. 


## snek.py

More specialized functions to investigate CobraPy models and to look at fluxes. For more information on available functions just do

```
import snek
dir(snek)
``` 

## snek_utils.py

Simple wrappers to important cobra functions that do additional checks:
* ensuring correct spelling,
* check for logical errors.

Includes following functions:

``` 
sensitive_optimize() # Returns an None object instead of dictionary with old fluxes during infeasible optimization.
set_bounds()         # Ensures that no spelling errors occur during definition of bounds.
set_objective()      # checks for unusual bounds on objective function. Ensures no spelling errors occur.
``` 
