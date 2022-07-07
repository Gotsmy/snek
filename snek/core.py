import warnings
'''
The functions listed here are wrappers to important cobra functions that do additional checks:
* ensuring correct spelling,
* check for logical errors.
'''

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

def sensitive_optimize(model):
    '''
    In contrast to the original implementation where fluxes can still be extracted from
    infeasible solutions (e.g. ``solution[reaction_id] = <float>``, even if ``solution.status == 'infeasible')``,
    this implementation returns a None object if the solver is infeasible.

    Parameters
    ----------
        model    : cobra.core.model.Model
                   CobraPy Model.

    Returns
    -------
        solution : cobra.core.solution.Solution or None
                   Optimized solution if ``solution.status != 'infeasible'``, otherwise None.
    '''

    solution = model.optimize()
    if solution.status == 'infeasible':
        return None
    else:
        return solution
