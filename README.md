# Introduction

Snek provides quality of Life Functions for CobraPy.


## 1. Core Module

There are three reasons why this module is helpful for CobraPy users:

### 1.1 Ensuring Correct Spelling

When setting bounds of a reaction or the objective function with the main CobraPy
package, no error is raised when a spelling error occurs.

```
Reaction.lower_bound = 5
Reaction.lower_buond = 5 # no error raised, but the bound is not updated
```

This can happen frequently and is hard to debug. Thus following functions ensure
correct spelling:

* ```set_bounds()```
* ```set_objective()```

### 1.2 Checking for Logical Errors

CobraPy gives the user a lot of freedom. However, here we implemented some basic
checks that ensure that no logical errors are made and prints warning if this is true.

* ```sensitive_optimize()```
* ```set_objective()```

### 1.3 Easy Programmatic Access

Certain functionalities are somewhat "hidden" in CobraPy. Here we implemented some
basic functions that easily access this hidden functionalities.

* ```get_objective()```
* ```get_solver()```
* ```find_biomass_reaction()```


## 2. Elements Module

This moduel contains several functions for the elemental analysis of molecules
and the calculation of molecular weights from sum formulas. For more information
on available functions just do

```
import snek
dir(snek.elements)
```
or visit the documentation.

## 3. Analysis Module

More specialized functions to investigate CobraPy models and to look at fluxes. For more information on available functions just do

```
import snek
dir(snek.analysis)
```
or visit the documentation.
