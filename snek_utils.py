import cobra

def set_objective(model,reaction,direction='max'):
    '''
    Updates the objective function of a model.
    -
    Input:
        model       cobrapy model
        reaction    objective ration ID
        direction   objective direction. default = 'max'
    -
    Output:
        model       updated cobrapy model
    '''
    
    model.objective = reaction
    model.objective_direction = direction
    return model

def set_bounds(model,reaction,lb=None,ub=None):
    '''
    Updates the bounds of a specific reaction.
    -
    Input: 
        model       cobrapy model
        reaction    String. Reaction ID
        lb          Number. Lower bound, default = None
        ub          Number. Upper bound, default = None
    -
    Output:
        model       updated cobrapy model
    '''

    if lb is None:
        lb = model.reactions.get_by_id(reaction).lower_bound
    if ub is None:
        ub = model.reactions.get_by_id(reaction).upper_bound
    
    model.reactions.get_by_id(reaction).bounds = float(lb),float(ub)
    return model

def sensitive_optimize(model):
    '''
    In contrast to the original implementation where fluxes can still be extracted from infeasible solutions (e.g. solution[reaction_id] = float, even if solution.status == 'infeasible'), this implementation returns a None object if the solver is infeasible.
    -
    Input:
        model               cobrapy model
    -
    Output:
        solution or None    Optimized solution if solution.status != 'infeasible', otherwise None.
    '''
    
    solution = model.optimize()
    if solution.status == 'infeasible':
        return None
    else:
        return solution