# Introduction

Snek provides quality of Life Functions for CobraPy. <br>
There are several reasons why this module is helpful for CobraPy users:

## 1. Core Module

### 1.1 Ensuring Correct Spelling

When setting bounds of a reaction or the objective function with the main CobraPy
package, no error is raised when a spelling error occurs.

```
Reaction.lower_bound = 5
Reaction.lower_buond = 5 # no error raised, but the bound is not updated
```

This can happen frequently and is hard to debug. Thus following functions ensure
correct spelling:

* ```snek.set_bounds()```
* ```snek.set_objective()```

### 1.2 Checking for Unintended Results

CobraPy gives the user a lot of freedom. However, here we implemented some basic
checks that ensure that no unintended errors are made and prints warning if this is true.

* ```snek.sensitive_optimize()```
* ```snek.set_objective()```

### 1.3 Easy Programmatic Access

Certain functionalities are somewhat "hidden" in CobraPy. Here we implemented some
basic functions that easily access this hidden functionalities.

* ```snek.get_objective()```
* ```snek.get_solver()```

### 1.4 New Functions

We added some new functions to easily search for reactions and metabolites as well.

* ```snek.find_biomass_reaction()```
* ```snek.find_reaction()```
* ```snek.find_metabolite()```

## 2. Elements Module

This module contains several functions for the elemental analysis of molecules
and the calculation of molecular weights from sum formulas. Following functions
are currently available:

* ```snek.elements.count_atom()```
* ```snek.elements.unique_elements()```
* ```snek.elements.formula_to_dictionary()```
* ```snek.elements.dictionary_to_formula()```
* ```snek.elements.element_flux_coefficient()```
* ```snek.elements.element_fluxes()```
* ```snek.elements.molecular_weight()```
* ```snek.elements.get_pDNA_formula()```
* ```snek.elements.get_protein_formula()```


## 3. Analysis Module

More specialized functions to investigate CobraPy models and to look at fluxes.
Following functions are currently available.

* ```snek.analysis.get_constrained_reactions()```
* ```snek.analysis.in_bounds()```
* ```snek.analysis.in_flux()```
* ```snek.analysis.investigate_network()```
* ```snek.analysis.investigate_network_solution()```
* ```snek.analysis.out_bounds()```
* ```snek.analysis.out_flux()```
