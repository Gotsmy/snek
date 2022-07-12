'''
This moduel contains several functions for the elemental analysis of molecules and
the calculation of molecular weights from sum formulas.
'''

import pkg_resources
import pandas as pd

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
    Calculates molecular weight from sum formula. The mass table is derived from the enviPat package [1].

    Parameters
    ----------
    sum_formula : Str
        String of chemical sum formula. Capitalization has to be correct.

    Returns
    -------
    molecular_weight : Float
        Float of molecular weight in g/mol.

    [1] Loos, M., Gerber, C., Corona, F., Hollender, J., Singer, H. (2015). Accelerated isotope fine structure calculation using pruned transition trees, Analytical Chemistry 87(11), 5738-5744. 
    '''

    elements = element_composition(sum_formula)

    # the molecular masses of the elements is downloaded from
    # https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee
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

    # leading strand
    for name, formula in zip(base_names,base_formulas_fw):
        tmp_count = sequence.count(name)
        tmp_elements = element_composition(formula)
        for element in tmp_elements:
                sum_formula_dic[element] += tmp_count * tmp_elements[element]
    # lagging strand
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
