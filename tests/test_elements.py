import pytest
import numpy as np
import cobra
from Bio import SeqUtils
import snek

@pytest.fixture
def model():
    model = cobra.io.load_model("textbook")
    return model

# check snek.elements.count_atom
def test_count_atom():
    nr = snek.elements.count_atom('C16H15N3.99','N')
    assert nr == 3.99

# check snek.elements.element_composition
def test_element_composition():
    all = snek.elements.element_composition('C16H15N3.99')
    assert all == {'C':16,'H':15,'N':3.99}

# check snek.elements.unique_elements
def test_unique_elements():
    elements = snek.elements.unique_elements('C16H15N3.99')
    assert len(elements) == 3

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

# check snek.elements.element_fluxes
def test_element_fluxes(model):
    C_fluxes = snek.elements.element_fluxes(model,'C')
    assert C_fluxes.loc['EX_glc__D_e','C_flux'] == 60
    H_fluxes = snek.elements.element_fluxes(model,'H')
    assert np.isclose(H_fluxes.loc['EX_h2o_e','H_flux'],-58.351654271131565,atol=model.tolerance,rtol=0)

# check snek.elements.element_flux_coefficient
def test_element_flux_coefficient(model):
    c_flux_coefficient = snek.elements.element_flux_coefficient(model,'C','EX_glc__D_e')
    assert c_flux_coefficient == -6

print('End of test_elements.py')
