'''
This moduel contains several functions for the elemental analysis of molecules and
the calculation of molecular weights from sum formulas.
'''

import pkg_resources
import numpy as np
import pandas as pd
from .core import sensitive_optimize

def count_atom(formula,element):
    """
    Counts the number of atoms of an element in a chemical formula.
    Works with multiple-letter element names and floats.
    Case sensitive!

    Parameters
    ----------
        formula  :  Str
                    Chemical formula, e.g. ``'H2O'``.
        element  :  Str
                    Element name, e.g. ``'H'``.

    Returns
    -------
        nr       :  float
                    Number of atoms as float, e.g. ``2.0``.
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

    Parameters
    ----------
        formula   : Str
                    Chemical formula, e.g. ``'H2O'``.

    Returns
    -------
        list_elements :  List
                         List of unique elements, e.g. ``['H','O']``.
    '''

    list_elements = []

    tmp = ''
    for letter in formula:
        if not letter.isnumeric() and not letter=='.':
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

    Parameters
    ----------
        formula  :  Str
                    Chemical formula, e.g. ``'H2O'``.

    Returns
    -------
        composition :   Dict
                        Dictionary with unique elements as keys
                        and number of atoms per element as items,
                        e.g. ``{'H': 2.0, 'O': 1.0}``.

    '''

    list_elements = unique_elements(formula)
    composition  = {}
    for element in list_elements:
        composition[element] = count_atom(formula,element)

    return composition

def molecular_weight(sum_formula):
    '''
    Calculates molecular weight from sum formula. The mass table is derived from the CobraPy package [1]_.

    Parameters
    ----------
    sum_formula : Str
        String of chemical sum formula. Capitalization has to be correct.

    Returns
    -------
    molecular_weight : Float
        Float of molecular weight in g/mol.

    Notes
    -----
    .. [1] Ebrahim, A., Lerman, J. A., Palsson, B. O., & Hyduke, D. R. (2013).
        COBRApy: constraints-based reconstruction and analysis for python.
        BMC systems biology, 7(1), 1-6.
        https://doi.org/10.1186/1752-0509-7-74
    '''

    elements = element_composition(sum_formula)

    # the molecular masses of the elements is downloaded from CobraPy
    # https://github.com/opencobra/cobrapy/blob/devel/src/cobra/core/formula.py
    location = pkg_resources.resource_stream(__name__,'data/mol_mass_table.csv')
    mol_mass_table = pd.read_csv(location)

    molecular_weight = 0
    for element in elements:
        molecular_weight += mol_mass_table[mol_mass_table['Symbol']==element]['AtomicMass'].values[0]*elements[element]

    return molecular_weight

def sum_formula_pDNA(sequence,return_dic=False):
    '''
    Returns the chemical sum formula of the pDNA molecule.
    Attention: this only works for circular, double stranded DNA!

    Parameters
    ----------
    sequence : Str
        String of pDNA sequence (only ``['A','T','C','G']`` are allowed).
    return_dic : Bool, default=False
        If True, additionally a dictionary with elements and counts in returned.

    Returns
    -------
    sum_formula : Str
        Chemical sum formula of the pDNA.
    sum_formula_dic : Dict
        Dictionary with elements and counts, only returned if return_dic=True.
    '''

    # H2O was subtracted from all dNMPs
    dAMP_delta_H2O = 'C10H12N5O5P'
    dTMP_delta_H2O = 'C10H13N2O7P'
    dCMP_delta_H2O =  'C9H12N3O6P'
    dGMP_delta_H2O = 'C10H12N5O6P'

    base_names = ['A','T','C','G']
    base_formulas_fw = [dAMP_delta_H2O,dTMP_delta_H2O,dCMP_delta_H2O,dGMP_delta_H2O]
    base_formulas_rv = [dTMP_delta_H2O,dAMP_delta_H2O,dGMP_delta_H2O,dCMP_delta_H2O]

    assert len(sequence) == sum([sequence.count(i) for i in base_names])

    sum_formula_dic = {'C':0,'H':0,'N':0,'O':0,'P':0}

    # forward strand
    for name, formula in zip(base_names,base_formulas_fw):
        tmp_count = sequence.count(name)
        tmp_elements = element_composition(formula)
        for element in tmp_elements:
                sum_formula_dic[element] += tmp_count * tmp_elements[element]
    # reverse strand
    for name, formula in zip(base_names,base_formulas_rv):
        tmp_count = sequence.count(name)
        tmp_elements = element_composition(formula)
        for element in tmp_elements:
                sum_formula_dic[element] += tmp_count * tmp_elements[element]

    sum_formula = ''
    for element in sum_formula_dic:
        sum_formula += element
        sum_formula += str(round(sum_formula_dic[element]))

    if return_dic:
        return sum_formula,sum_formula_dic
    else:
        return sum_formula

def element_flux_coefficient(model,element,reaction):
    """
    Returns the amount of an user-defined chemical element consumed/produced in a CobraPy reaction.

    Parameters
    ----------
    model : cobra.core.model.Model
        CobraPy model.
    element : Str
        Chemical element as string, e.g. 'C'.
    reaction : Str
        Cobra reaction id as string.

    Returns
    -------
    element_flux_coefficient : Float
        Elemental flux of the reaction (only non-0 when not mass balanced).
    """

    element_flux_coefficient = 0
    for meta in model.reactions.get_by_id(reaction).metabolites:
        stoichiometric_coefficient = model.reactions.get_by_id(reaction).metabolites[meta]
        number_of_atoms = count_atom(meta.formula,element)
        element_flux_coefficient += stoichiometric_coefficient * number_of_atoms

    return element_flux_coefficient

def element_fluxes(model,element,solution=None):
    """
    Returns dictionary of non-zero element fluxes.
    Element fluxes are only non-zero when a reaction is not mass balanced,
    typically exchange, drain, and sink reactions.

    If no solution object is provided a pFBA is performed.
    If a solution object is provided metabolite fluxes are extracted from there.
    We recommend to do this analysis only solutions from pFBA.

    Parameters
    ----------
    model : cobra.core.model.Model
            CobraPy model.
    element : Str
        Chemical element as string, e.g. 'C'.
    solution : cobra.core.solution.Solution, default=None
        CobraPy solution object.

    Returns
    -------
    element_fluxes : pandas.DataFrame
               Pandas dataframe of non-zero element fluxes. The sign of fluxes
               corresponds to the sign of fluxes in the solution object.
               The index of the data frame corresponds to reaction ids,
               the first column to the absolute element flux in mmol/(g biomass h),
               and the second column to the relative element flux in %.
    """

    if isinstance(solution,type(None)):
        solution = sensitive_optimize(model,pFBA=True)

    # calculate fluxes
    element_fluxes = {}
    for r in model.reactions:
        rate = element_flux_coefficient(model,element,r.id)
        if rate != 0 and solution.fluxes[r.id] != 0:
            element_fluxes[r.id] = [solution.fluxes[r.id]*rate]

    # raise error if there are no fluxes in the solution
    if len(element_fluxes) == 0:
        raise ValueError(f'No {element} fluxes were found.')

    # do a sanity check
    control = 0
    for flux in element_fluxes:
        control += element_fluxes[flux][0]
    assert np.isclose(control,0,atol=model.tolerance,rtol=0), 'Fluxes do not add up. Are all non-boundary reactions mass balanced?'

    # create nice dataframe
    element_fluxes = pd.DataFrame(element_fluxes).T
    abs_col_name = f'{element}_flux'
    element_fluxes.columns = [abs_col_name]
    element_fluxes = element_fluxes.sort_values(abs_col_name,ascending=False)
    total_flux = np.sum(element_fluxes[element_fluxes[abs_col_name] > 0].values)
    element_fluxes.insert(1,abs_col_name+'%',element_fluxes[abs_col_name]/total_flux*100)

    return element_fluxes
