import cobra
import numpy as np
import pandas as pd
from . import snek_utils

def count_atom(formula,element):
    """
    Counts the number of atoms of an element in a chemical formula. 
    Works with multiple-letter element names and floats. 
    Case sensitive!
    -
    Input:
        formula    Str. Chemical formula, e.g. 'H2O'.
        element    Str. Element name, e.g. 'H'.
    -
    Output:
        nr         Num. Number of atoms as float, e.g. 2.0
    """
    found = 0
    while formula.find(element) != -1:
        idx = formula.find(element)+len(element)-1
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

def unique_elements(formula):
    '''
    Get a list of unique elements from a chemical formula.
    -
    Input:
        formula    Str. Chemical formula, e.g. 'H2O'.
    -
    Output:
        list_elements    List. List of unique elements, e.g. ['H','O'].
    '''
    
    list_elements = []
    
    tmp = ''
    for letter in formula:
        if not letter.isnumeric():
            if letter.isupper():
                if len(tmp) > 0 and tmp not in list_elements:
                    list_elements.append(tmp)
                tmp = ''
                tmp += letter
            else:
                tmp += letter
            
    if len(tmp) > 0 and tmp not in list_elements:
        list_elements.append(tmp)
    
    return list_elements

def element_composition(formula):
    '''
    Get the elemental composition of a molecule as dictionary. 
    Keys correspond to elements, values to their number of appearance.
    -
    Input:
        formula    Str. Chemical formula, e.g. 'H2O'.
    -
    Output:
        composition    Dict. Dictionary with unique elements as keys
                       and number of atoms per element as items,
                       e.g. {'H': 2.0, 'O': 1.0}.

    '''

    list_elements = unique_elements(formula)
    composition  = {}
    for element in list_elements:
        composition[element] = count_atom(formula,element)

    return composition

def get_constrained_reactions(model):
    '''
    Returns data frames of all reactions that have constrained flux bounds.
    I.e. something else than (lower bound, upper bound)
    * -1000/1000,
    *     0/1000, or
    * -1000/   0.
    -
    Input:
        model    CobraPy model.
    -
    Output:
        constrained    DataFrame. Pandas data frame of constrained reactions with name, 
                       id, lower bound, upper bound as columns.
        blocked        DataFrame. Pandas data frame of constrained reactions with name, 
                       id, lower bound, upper bound as columns.
    '''
    
    constrained = []
    blocked = []
    for reaction in model.reactions:
        if reaction.bounds == (0,1000):
            pass
        elif reaction.bounds == (-1000,1000):
            pass
        elif reaction.bounds == (-1000,0):
            pass
        elif reaction.bounds == (0,0):
            blocked.append([reaction.name,reaction.id,*reaction.bounds])
        else:
            constrained.append([reaction.name,reaction.id,*reaction.bounds])    
    blocked = pd.DataFrame(blocked,columns=['name','id','lower_bound','upper_bound'])
    constrained = pd.DataFrame(constrained,columns=['name','id','lower_bound','upper_bound'])
    
    return constrained, blocked

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

def in_flux(model,pFBA=True):
    """ 
    Returns all exchange reactions with flux < 0  as dictionary.
    -
    Input:
        model    cobrapy model.
        pFBA     If True a pFBA is performed to get fluxes. default = True.
    Output:
        in_flux  dictionary with all exchange reactions with flux < 0.
    
    """
    
    in_flux = {}
    if pFBA:
        solution = cobra.flux_analysis.pfba(model)
    else:
        solution = snek_utils.sensitive_optimize(model)
    for ex in model.exchanges:
        if solution[ex.id] < 0:
            in_flux[ex.id] = solution[ex.id]
    return in_flux

def out_flux(model,pFBA=True):
    """ 
    Returns all exchange reactions with flux > 0  as dictionary.
    -
    Input:
        model    cobrapy model.
        pFBA     If True a pFBA is performed to get fluxes. default = True.
    Output:
        out_flux dictionary with all exchange reactions with flux > 0.
    
    """
    
    out_flux = {}
    if pFBA:
        solution = cobra.flux_analysis.pfba(model)
    else:
        solution = snek_utils.sensitive_optimize(model)
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

def investigate_reaction(model,reaction_id,silent_metabolites=None):
    '''
    Returns metabolites which are part of the reaction.
    Exempt for metabolite ids listed in silent_metabolites 
    (default: [}).
    '''
    if silent_metabolites == None:
        silent_metabolites = []
    metabolite_list = []
    for m in model.reactions.get_by_id(reaction_id).metabolites:
        if m.id not in silent_metabolites:
            metabolite_list.append(m.id)
    
    return metabolite_list

def investigate_metabolite(model,metabolite_id,silent_reactions=None):
    '''
    Returns reactions in which this metabolite is present.
    Exempt for reaction ids listed in silent_reactions
    (default: []).
    '''
    if silent_reactions == None:
        silent_reactions = []
    reaction_list = []
    for r in model.metabolites.get_by_id(metabolite_id).reactions:
        if r.id not in silent_reactions:
            reaction_list.append(r.id)
    
    return reaction_list

def mimic_excel():
    'Just creates some letters, Excel style.'
    # https://stackoverflow.com/questions/63875471/enumerate-with-letters-instead-of-numbers
    for i in range(0, 26):
        yield chr(i + 65)

    i, j = [0, 0]

    for j in range(0, 26):
        for i in range(0, 26):
            yield "{}{}".format(chr(j + 65), chr(i + 65))
            
def investigate_network_solution(model, solution, start_reaction, depth, silent_reactions = None, silent_metabolites = None):
    '''
    See the enviroment of the Metabolic Network around a reaction of interest.
    
    Input:
        model                 cobrapy model
        solution              a cobrapy model solution
        start_reaction        reaction_id of reaction of interest
        depth                 reaction depth to investigate
        silent_reactions      reactions you are not interested in for whatever reason
                              default: []
        silent_metabolites    metabolites you are not interested in for whatever reason
                              default: ['h_c','atp_c','h2o_c','pi_c','h_p','adp_c','o2_p','h2o_c','h2o_p']
    Output:
        None - Prints reactions and metabolites in the network with non-zero flux in the solution object.
    '''

    if silent_reactions == None:
        silent_reactions = []
    if silent_metabolites == None:
        silent_metabolites = ['h_c','atp_c','h2o_c','pi_c','h_p','adp_c','o2_p','h2o_c','h2o_p']

    reaction_list    = [start_reaction]
    metabolite_list  = []
    for i in range(1,depth+1):
        new_reaction_list = []
        for n,r in zip(mimic_excel(),reaction_list):
            if r in silent_reactions:
                pass
            elif solution[r] == 0:
                pass
            else:
                print('{:1}.{:3} {:17} ({:.5f})'.format(i,n,r,solution[r]))
                metabolite_list = investigate_reaction(model,r,silent_metabolites)
                silent_reactions.append(r)
                for m in metabolite_list:
                    print('      > {:15} ({:5})'.format(m,model.reactions.get_by_id(r).metabolites[model.metabolites.get_by_id(m)]))
                    new_reaction_list += investigate_metabolite(model,m,silent_reactions)
        reaction_list = set(new_reaction_list)
    return

def investigate_network(model,start_reaction,depth,silent_reactions = None,silent_metabolites = None):
    '''
    See the enviroment of the Metabolic Network around a reaction of interest.
    
    Input:
        model                 cobrapy model
        start_reaction        reaction_id of reaction of interest
        depth                 reaction depth to investigate
        silent_reactions      reactions you are not interested in for whatever reason
                              default: []
        silent_metabolites    metabolites you are not interested in for whatever reason
                              default: ['h_c','atp_c','h2o_c','pi_c','h_p','adp_c','o2_p','h2o_c','h2o_p']
    Output:
        None - Prints reactions and metabolites in the Network.    '''

    if silent_reactions == None:
        silent_reactions = []
    if silent_metabolites == None:
        silent_metabolites = ['h_c','atp_c','h2o_c','pi_c','h_p','adp_c','o2_p','h2o_c','h2o_p']

    reaction_list    = [start_reaction]
    metabolite_list  = []
    for i in range(1,depth+1):
        new_reaction_list = []
        for n,r in zip(mimic_excel(),reaction_list):
            if r in silent_reactions:
                pass
            else:
                print('{:1}.{:3} {}'.format(i,n,r))
                metabolite_list = investigate_reaction(model,r,silent_metabolites)
                silent_reactions.append(r)
                for m in metabolite_list:
                    print('      > {:15} ({:5})'.format(m,model.reactions.get_by_id(r).metabolites[model.metabolites.get_by_id(m)]))
                    new_reaction_list += investigate_metabolite(model,m,silent_reactions)
        reaction_list = set(new_reaction_list)
    return
