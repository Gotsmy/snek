import cobra
import numpy as np

def count_atom(formula,atom):
    """
    Counts the atoms in a chemical formula. 
    Works with multiple-letter atom names and floats. 
    Case sensitive!
    INPUT  
        formula    Chemical Formula, e.g. H2O
        atom       Atom name, e.g. H
    RETURNS 
        Number of atoms as float, e.g. 2.0
    """
    found = 0
    while formula.find(atom) != -1:
        idx = formula.find(atom)+len(atom)-1
        tmp = ''
        for f in formula[idx+1:]:
            if f.isupper() == False:
                tmp += f
            else:
                break
        if len(tmp) == 0:
            nr = '1'
            found = 1
            break
        
        elif tmp[0].islower() == False:
            nr = tmp
            found = 1
            break
        else:
            found = 0
            formula = formula[idx+1:]
            
    if found == 0:
        nr = '0'
    return float(nr)

def in_bounds(model):
    """ returns all exchange reactions with lower_bounds < 0  as dictionary """
    in_bnds ={}
    for ex in model.exchanges:
        if ex.lower_bound != 0:
            in_bnds[ex.id] = ex.lower_bound
    return in_bnds

def out_bounds(model):
    """ returns all exchange reactions with lower_bounds > 0  as dictionary """
    out_bnds ={}
    for ex in model.exchanges:
        if ex.upper_bound != 0:
            out_bnds[ex.id] = ex.upper_bound
    return out_bnds

def io_bounds(model,info=True):
    """ returns all exchange reactions that have non 0 input and their output bounds """
    io_bnds ={}
    if info == True:
        print('[lower_bound,upper_bound]')
    for ex in model.exchanges:
        if ex.lower_bound != 0:
            io_bnds[ex.id] = [ex.lower_bound,ex.upper_bound]
    return io_bnds

def in_flux(model):
    """ returns all exchange reactions with flux < 0  as dictionary """
    in_flux = {}
    solution = model.optimize()
    for ex in model.exchanges:
        if solution[ex.id] < 0:
            in_flux[ex.id] = solution[ex.id]
    return in_flux

def out_flux(model):
    """ returns all exchange reactions with flux > 0  as dictionary """
    out_flux = {}
    solution = model.optimize()
    for ex in model.exchanges:
        if solution[ex.id] > 0:
            out_flux[ex.id] = solution[ex.id]
    return out_flux

def c_flux_coefficient(reaction_id,model):
    """ returns the amount of C consumed/produced in a non-mass-balanced reaction """
    carbon_flux_coefficient = 0
    for meta in model.reactions.get_by_id(reaction_id).metabolites:
        carbon_flux_coefficient += model.reactions.get_by_id(reaction_id).metabolites[meta]*count_atom(meta.formula,'C')
    return carbon_flux_coefficient

def carbon_fluxes(model):
    """ Returns dictionary of non-zero carbon fluxes during FBA.
    Carbon fluxes are only non-zero when a reaction is not mass balanced. """
    c_fluxes = {}
    solution = model.optimize()
    for r in model.reactions:
        rate = c_flux_coefficient(r.id,model)
        if rate != 0 and solution.fluxes[r.id] != 0:
            c_fluxes[r.id] = solution.fluxes[r.id]*rate
    #do a sanity check
    control = 0
    for c_flux in c_fluxes:
        control += c_fluxes[c_flux]
    if round(control,5) != 0:
        print("There is a problem! Fluxes don't add up.")
        print('Sum :',control)
        return
    elif round(control,10) != 0:
        print("There might be a problem! Fluxes don't add up.")
        print('Sum :',control)
    return c_fluxes

def pfba_carbon_fluxes(model):
    """ Returns dictionary of non-zero carbon fluxes during pFBA.
    Carbon fluxes are only non-zero when a reaction is not mass balanced. """
    c_fluxes = {}
    solution = cobra.flux_analysis.pfba(model)
    for r in model.reactions:
        rate = c_flux_coefficient(r.id,model)
        if rate != 0 and solution.fluxes[r.id] != 0:
            c_fluxes[r.id] = solution.fluxes[r.id]*rate
    #do a sanity check
    control = 0
    for c_flux in c_fluxes:
        control += c_fluxes[c_flux]
    if round(control,5) != 0:
        print("There is a problem! Fluxes don't add up.")
        print('Sum :',control)
        return
    elif round(control,10) != 0:
        print("There might be a problem! Fluxes don't add up.")
        print('Sum :',control)
    return c_fluxes

def michaelis_menten(c,vmax,km):
    """returns negative v. be aware of units!"""
    if c <= 0:
        v = 0
    else:
        v = -vmax * c / (km + c)
    return v

def lex_dfba(model,compounds,y_zero,time,objectives,objectives_direction,dynamic_constraints):
    
    """ \nLexicographic Dynamic FBA\n
    INPUT
    \t compounds  : list of reactions to track with size n
    \t y_zero     : list of start concentrations of reaction metabolites with size n
    \t time       : list with t0, tmax, and dt
    \t objectives : list of reactions to optimize with size m
    \t objectives_direction : list of direction of reactions to optimize ('max' or 'min') with size m
    \t dynamic_constraints  : list of lists in form of [['reaction_id',vmax,km],...]
    
    RETURNS
    \t dict of tracked reaction concentrations
    \t dict of tracked reaction fluxes"""
   

    y = {}
    f = {}
 
    compounds = compounds # compound exchange reactions
    y[0] = y_zero # compound concentrations
    t0,tmax,dt = time # start time, max time, time step size
    objectives = objectives # list of reactions for objective function
    objectives_direction = objectives_direction # direction of objetives
    
    stat = ''
    for n,t in enumerate(np.arange(t0,tmax,dt)):
        with model:
            for dc in dynamic_constraints:
                model.reactions.get_by_id(dc[0]).lower_bound = michaelis_menten(y[n][compounds.index(dc[0])],dc[1],dc[2])
            try:
                lex_constraints = cobra.util.add_lexicographic_constraints(model, objectives, objectives_direction)
            except:
                stat = 'infeasible'
            solution = model.optimize()
            if model.solver.status == 'infeasible' or stat == 'infeasible':
                break
            tmp_y = []
            tmp_f = []
            for nr, c in enumerate(compounds):
                tmp_y.append(y[n][nr]+solution.fluxes[c]*dt*y[n][0])
                tmp_f.append(solution.fluxes[c])
            y[n+1]=tmp_y
            f[n+1]=tmp_f
    return y,f
