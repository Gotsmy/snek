'''
This module contains several functions for the analysis of CobraPy models.
'''

import cobra
import numpy as np
import pandas as pd
from .core import sensitive_optimize
from .elements import *

def get_constrained_reactions(model):
    '''
    Returns data frames of all reactions that have constrained flux bounds.
    I.e. something else than (lower bound, upper bound)

    * -1000/1000,
    *     0/1000, or
    * -1000/   0.

    Parameters
    ----------
        model   : cobra.core.model.Model
                  CobraPy model.

    Returns
    -------
        constrained  :  pandas.DataFrame
                        Pandas data frame of constrained reactions with name,
                        id, lower bound, upper bound as columns.
        blocked      :  pandas.DataFrame
                        Pandas data frame of constrained reactions with name,
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
    """
    Returns all exchange reactions with lower_bounds < 0  as dictionary.

    Parameters
    ----------
    model : cobra.core.model.Model
            CobraPy Model.

    Returns
    -------
    in_bnds : Dict
              Dictionary of all exchange reactions that have lower bounds lower than 0.

    """
    in_bnds ={}
    for ex in model.exchanges:
        if ex.lower_bound != 0:
            in_bnds[ex.id] = ex.lower_bound
    return in_bnds

def out_bounds(model):
    """
    Returns all exchange reactions with lower_bounds > 0  as dictionary.

    Parameters
    ----------
    model : cobra.core.model.Model
            CobraPy Model.

    Returns
    -------
    out_bnds : Dict
              Dictionary of all exchange reactions that have lower bounds greater than 0.

    """
    out_bnds ={}
    for ex in model.exchanges:
        if ex.upper_bound != 0:
            out_bnds[ex.id] = ex.upper_bound
    return out_bnds

def io_bounds(model,info=True):
    """Returns all exchange reactions that have non 0 input and their output bounds

    Parameters
    ----------
    model : cobra.core.model.Model
            CobraPy Model.
    info :  Bool, default=True
            If True prints the name of variables in output list.

    Returns
    -------
    io_bnds : Dict
              Dictionary of exchange reactions with non-zero lower bound in the form of
              ``io_bnds[<reaction_id>] = [<lower bound>,<upper bound>]``.

    """
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

    Parameters
    ----------
        model    : cobra.core.model.Model
                   Cobrapy model.
        pFBA     : Bool, default=True
                   If True a pFBA is performed to get fluxes, else FBA.

    Returns
    -------
        in_flux : Dict
                  Dictionary with all exchange reactions with flux < 0.

    """

    in_flux = {}
    if pFBA:
        solution = cobra.flux_analysis.pfba(model)
    else:
        solution = sensitive_optimize(model)
    for ex in model.exchanges:
        if solution[ex.id] < 0:
            in_flux[ex.id] = solution[ex.id]
    return in_flux

def out_flux(model,pFBA=True):
    """
    Returns all exchange reactions with flux > 0  as dictionary.

    Parameters
    ----------
        model  : cobra.core.model.Model
                 Cobrapy model.
        pFBA   : Bool, default=True
                 If True a pFBA is performed to get fluxes.

    Returns
    -------
        out_flux : Dict
                   Dictionary with all exchange reactions with flux > 0.

    """

    out_flux = {}
    if pFBA:
        solution = cobra.flux_analysis.pfba(model)
    else:
        solution = sensitive_optimize(model)
    for ex in model.exchanges:
        if solution[ex.id] > 0:
            out_flux[ex.id] = solution[ex.id]
    return out_flux

def c_flux_coefficient(reaction_id,model):
    """
    Returns the amount of C consumed/produced in a non-mass-balanced reaction

    Parameters
    ----------
    reaction_id : Str
                  Cobra reaction id as string.
    model       : cobra.core.model.Model
                  CobraPy model.

    Returns
    -------
    carbon_flux_coefficient : Float
                              C flux of the reaction (only non-0 when not mass balanced).

    """
    carbon_flux_coefficient = 0
    for meta in model.reactions.get_by_id(reaction_id).metabolites:
        carbon_flux_coefficient += model.reactions.get_by_id(reaction_id).metabolites[meta]*count_atom(meta.formula,'C')
    return carbon_flux_coefficient

def carbon_fluxes(model):
    """
    Returns dictionary of non-zero carbon fluxes during FBA.
    Carbon fluxes are only non-zero when a reaction is not mass balanced.

    Parameters
    ----------
    model : cobra.core.model.Model
            CobraPy model.

    Returns
    -------
    c_fluxes : Dict
               Dictionary of non-zero carbon fluxes during FBA.

    """
    c_fluxes = {}
    solution = sensitive_optimize(model)
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
    Carbon fluxes are only non-zero when a reaction is not mass balanced.

    Parameters
    ----------
    model : cobra.core.model.Model
            CobraPy model.

    Returns
    -------
    c_fluxes : Dict
               Dictionary of non-zero carbon fluxes during pFBA.
    """
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
    """
    Simple implementation of a Michaelis Menten enzyme kinetic.
    v = - vmax*c/(k_M+c)
    Attention: returns negative fluxes.
    Attention: be aware of units.


    Parameters
    ----------
    c : float
        Concentration.
    vmax : float
        Maximum reaction rate.
    km : float
        K_M as defined in Michaelis Menten kinetics.

    Returns
    -------
    v : float
        Reaction rate, (negative sign!).

    """
    if c <= 0:
        v = 0
    else:
        v = -vmax * c / (km + c)
    return v

def lex_dfba(model,compounds,y_zero,time,objectives,objectives_direction,dynamic_constraints):
    """
    Lexicographic Dynamic FBA

    Parameters
    ----------
    compounds  : list
        List of reactions to track with size n
    y_zero     : list
        List of start concentrations of reaction metabolites with size n
    time       : list
        List with t0, tmax, and dt
    objectives : list
        List of reactions to optimize with size m
    objectives_direction : list
        List of direction of reactions to optimize ('max' or 'min') with size m
    dynamic_constraints  : list
        List of lists in form of ``[['reaction_id',vmax,km],...]``

    Returns
    -------
    y : dict
        Dictionary of tracked reaction concentrations
    f : dict
        Dictionary of tracked reaction fluxes
        """


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

    Parameters
    ----------


    Returns
    -------

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

    Parameters
    ----------

    Returns
    -------

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

    Parameters
    ----------
    model : cobra.core.model.Model
        Cobrapy model.
    solution : cobra.core.Solution
        A CobraPy model solution.
    start_reaction : Str
        Reaction_id of reaction of interest
    depth : Int
        Reaction depth to investigate
    silent_reactions : list ``default=[]``
        Reactions you are not interested in for whatever reason
    silent_metabolites : list, ``default=['h_c','atp_c','h2o_c','pi_c','h_p','adp_c','o2_p','h2o_c','h2o_p']``
        Metabolites you are not interested in for whatever reason (usually because they are present in almost all reactions).


    Returns
    -------
    None : None
        Prints reactions and metabolites in the network with non-zero flux in the solution object.
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

    Parameters
    ----------
    model : cobra.core.model.Model
        Cobrapy model.
    start_reaction : Str
        Reaction_id of reaction of interest
    depth : Int
        Reaction depth to investigate
    silent_reactions : list ``default=[]``
        Reactions you are not interested in for whatever reason
    silent_metabolites : list, ``default=['h_c','atp_c','h2o_c','pi_c','h_p','adp_c','o2_p','h2o_c','h2o_p']``
        Metabolites you are not interested in for whatever reason (usually because they are present in almost all reactions).

    Returns
    -------
    None : None
        Prints reactions and metabolites in the Network.    '''

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
