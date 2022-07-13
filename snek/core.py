'''
The core module contains simple wrappers to important cobra functions that do additional checks:

* ensuring correct spelling,
* check for logical errors.

The functions can also be called with ``snek.<function>``.
'''

import warnings
import cobra
from optlang.glpk_interface import Model as glpkModel

def set_objective(model,reaction,direction='max'):
    '''
    Updates the objective function of a model. More importantly, also checks if the objective reaction has unusual bounds.
    I.e. something else than (lower bound, upper bound)

    * -1000/1000 for maximization and minimization,
    *     0/1000 for maximization, or
    * -1000/   0 for minimization.

    Parameters
    ----------
        model       : cobra.core.model.Model
                      CobraPy Model
        reaction    : Str
                      Objective reaction ID.
        direction   : Str, default='max'
                      Objective direction, either ``'min'`` or ``'max'``.

    Returns
    -------
        model       : cobra.core.model.Model
                      Updated CobraPy Model.
    '''

    if model.reactions.get_by_id(reaction).upper_bound == 1000 and model.reactions.get_by_id(reaction).lower_bound == -1000:
        pass
    elif direction == 'min' and (model.reactions.get_by_id(reaction).lower_bound != -1000 or model.reactions.get_by_id(reaction).upper_bound != 0):
        warnings.warn('The objective reaction has a lower bound at {:.2f} and a upper bound at {:.2f}. This could lead to problems during optimization.'.format(model.reactions.get_by_id(reaction).lower_bound,model.reactions.get_by_id(reaction).upper_bound))
    elif direction == 'max' and (model.reactions.get_by_id(reaction).lower_bound != 0 or model.reactions.get_by_id(reaction).upper_bound != 1000):
        warnings.warn('The objective reaction has a lower bound at {:.2f} and a upper bound at {:.2f}. This could lead to problems during optimization.'.format(model.reactions.get_by_id(reaction).lower_bound,model.reactions.get_by_id(reaction).upper_bound))

    model.objective = reaction
    model.objective_direction = direction
    return model

def set_bounds(model,reaction,lb=None,ub=None):
    '''
    Updates the bounds of a specific reaction. If None, the bounds are not changed.

    Parameters
    ----------
        model      : cobra.core.model.Model
                     CobraPy Model.
        reaction   : Str
                     Reaction ID
        lb         : Float, default=None
                     Lower bound.
        ub         : Float, default=None
                     Upper bound.

    Returns
    -------
        model      : cobra.core.model.Model
                     Updated CobraPy Model.
    '''

    if lb is None:
        lb = model.reactions.get_by_id(reaction).lower_bound
    if ub is None:
        ub = model.reactions.get_by_id(reaction).upper_bound

    model.reactions.get_by_id(reaction).bounds = float(lb),float(ub)
    return model

def sensitive_optimize(model,pFBA=False):
    '''
    This function does some additional checks before model optimization.

    #. In contrast to the original implementation where fluxes can still be extracted from
       infeasible solutions (e.g. ``solution[reaction_id] = <float>``, even if
       ``solution.status == 'infeasible')``, this function raises a
       ValueError if the solver is infeasible.

    #. Additionally, this function checks if the objective reaction has unusual bounds.
       I.e., if other bounds than (lower bound, upper bound)

        * -1000/1000 for maximization and minimization,
        *     0/1000 for maximization, or
        * -1000/   0 for minimization.

       are present, a warning is printed.

    #. Finally, this function warns about inconsistent results that can occur when
       doing pFBA with the GLPK solver.

    Parameters
    ----------
        model    : cobra.core.model.Model
                   CobraPy Model.
        pFBA : Bool, default=False
            Wether or not to do parsimonious FBA.

    Returns
    -------
        solution : cobra.core.solution.Solution
                   Optimized solution if ``solution.status != 'infeasible'``.
    '''

    objective_reaction_id = get_objective(model)
    objective_direction = model.objective_direction
    model = set_objective(model,objective_reaction_id,objective_direction)

    if pFBA:
        if isinstance(model.solver,glpkModel):
            warnings.warn('You are performing pFBA with the GLPK Solver. This can lead to inconsistent results.')
        solution = cobra.flux_analysis.parsimonious.pfba(model)
    else:
        solution = model.optimize()
    if solution.status == 'infeasible':
        raise ValueError('Solution is infeasible')
    else:
        return solution

def find_biomass_reaction(model):
    '''
    Searches 'Biomass' string in reaction names or reaction IDs.
    Prints reaction IDs of found reactions.
    Not sensitive to capitalization.

    Parameters
    ----------
    model : cobra.core.model.Model
        CobraPy Model to be searched.
    '''

    for reaction in model.reactions:
        if 'BIOMASS' in reaction.name.upper() or 'BIOMASS' in reaction.id.upper():
            print(reaction.id)

def get_objective(model):
    '''
    Returns the reaction ID of a single reaction objective of a CobraPy model.

    Parameters
    ----------
    model : cobra.core.model.Model
        CobraPy Model.

    Returns
    -------
    reaction_id : Str
        Reaction id of the objective reaction.
    '''

    reaction_id = str(model.objective.expression).split('*')[1].split(' ')[0]

    return reaction_id
