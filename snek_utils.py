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
    
    model.reactions.get_by_id(reaction).bounds = lb,ub
    return model