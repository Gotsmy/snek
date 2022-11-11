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

# check snek.elements.element_composition
def test_element_composition():
    composition = snek.elements.element_composition('C16H15N3.99')
    assert composition == {'C':16,'H':15,'N':3.99}

# check snek.elements.element_formula
def test_element_formula():
    formula = snek.elements.element_formula({'C':16,'H':15,'N':3.99})
    assert formula == 'C16H15N3.99'

# check snek.elements.molecular_weight
def test_molecular_weight():
    formula = 'C6H12O6'
    dummy = cobra.Metabolite()
    dummy.formula = formula
    assert np.isclose(dummy.formula_weight,snek.elements.molecular_weight(formula))

# check snek.elements.sum_formula_pDNA
def test_sum_formula_pDNA():
    # test if my two functions work as intended and compare results to the Bio.SeqUtils.molecular_weight function.
    # the big difference: with my functions I also get the sum formula of the pDNA
    sequence = 'ATCG'*100
    test_sum_formula, test_sum_formula_dic = snek.elements.sum_formula_pDNA(sequence,return_dic=True)
    test_mw_my_implementation = snek.elements.molecular_weight(test_sum_formula)
    test_mw_bio = SeqUtils.molecular_weight(sequence,seq_type='DNA',double_stranded=True,circular=True,monoisotopic=False)
    # assert that the error is smaller than 0.1
    assert np.abs((test_mw_my_implementation-test_mw_bio)/np.mean([test_mw_bio,test_mw_my_implementation])*100) < .1

# check snek.elements.sum_formula_protein
def test_sum_formula_protein():
    # check sum formula of a simple protein 'ACLI'
    # sum formula calculated with https://spin.niddk.nih.gov/clore/Software/A205.html
    true_sum_formula = 'C18H34N4O5S1'
    protein_sequence = 'ACLI'
    snek_sum_formula = snek.elements.sum_formula_protein(protein_sequence)
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

print('End of test_elements.py')
