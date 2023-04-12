import os
import pytest
import numpy as np
import cobra
from Bio import SeqUtils
import snek

@pytest.fixture
def model():
    location = os.path.dirname(os.path.abspath(__file__))
    model = cobra.io.read_sbml_model(os.path.join(location,'data','textbook.xml'))
    return model

# check snek.elements.count_atom
def test_count_atom():
    nr = snek.elements.count_atom('C16H15N3.99','N')
    assert nr == 3.99

# check snek.elements.unique_elements
def test_unique_elements():
    elements = snek.elements.unique_elements('C16H15N3.99')
    assert len(elements) == 3

# check snek.elements.formula_to_dictionary
def test_formula_to_dictionary():
    dictionary = snek.elements.formula_to_dictionary('C16H15N3.99')
    assert dictionary == {'C':16,'H':15,'N':3.99}

# check snek.elements.element_formula
def test_dictionary_to_formula():
    formula = snek.elements.dictionary_to_formula({'C':16,'H':15,'N':3.99})
    assert formula == 'C16H15N3.99'

# check snek.elements.molecular_weight
def test_molecular_weight():
    formula = 'C6H12O6'
    dummy = cobra.Metabolite()
    dummy.formula = formula
    assert np.isclose(dummy.formula_weight,snek.elements.molecular_weight(formula))

# check snek.elements.get_pDNA_formula
def test_get_pDNA_formula():
    # test if my two functions work as intended and compare results to the Bio.SeqUtils.molecular_weight function.
    # the big difference: with my functions I also get the sum formula of the pDNA
    sequence = 'ATCG'*100
    test_sum_formula, test_sum_formula_dic = snek.elements.get_pDNA_formula(sequence,return_dict=True)
    test_mw_my_implementation = snek.elements.molecular_weight(test_sum_formula)
    test_mw_bio = SeqUtils.molecular_weight(sequence,seq_type='DNA',double_stranded=True,circular=True,monoisotopic=False)
    # assert that the error is smaller than 0.1
    assert np.abs((test_mw_my_implementation-test_mw_bio)/np.mean([test_mw_bio,test_mw_my_implementation])*100) < .1

# check snek.elements.get_protein_formula
def test_get_protein_formula():
    # check sum formula of a simple protein 'ACLI'
    # sum formula calculated with https://spin.niddk.nih.gov/clore/Software/A205.html
    true_sum_formula = 'C18H34N4O5S1'
    protein_sequence = 'ACLI'
    snek_sum_formula = snek.elements.get_protein_formula(protein_sequence)
    assert true_sum_formula == snek_sum_formula

# check snek.elements.element_flux_coefficient
def test_element_flux_coefficient(model):
    c_flux_coefficient = snek.elements.element_flux_coefficient(model,'C','EX_glc__D_e')
    assert c_flux_coefficient == -6

# check snek.elements.element_fluxes
def test_element_fluxes(model):
    C_fluxes = snek.elements.element_fluxes(model,'C')
    assert C_fluxes.loc['EX_glc__D_e','C_flux'] == 60
    H_fluxes = snek.elements.element_fluxes(model,'H')
    assert np.isclose(H_fluxes.loc['EX_h2o_e','H_flux'],-58.351654271131565,atol=model.tolerance,rtol=0)


# check snek.elements.get_unique_elements
def test_get_unique_elements(model):
    assert snek.elements.get_unique_elements(model) == ['C', 'H', 'O', 'P', 'N', 'S']
    assert snek.elements.get_unique_elements(model.reactions.ACKr) == ['C', 'H', 'O', 'N', 'P']
    assert snek.elements.get_unique_elements(model.metabolites.ac_c) == ['C', 'H', 'O']
    assert snek.elements.get_unique_elements("H2O") == ['H', 'O']

print('End of test_elements.py')
