import cobra
import snek
from Bio import SeqUtils
import numpy as np

# check count_atom
nr = snek.elements.count_atom('C16H15N3.99','N')
assert nr == 3.99

# check element_composition
all = snek.elements.element_composition('C16H15N3.99')
assert all == {'C':16,'H':15,'N':3.99}

# check unique_elements
elements = snek.elements.unique_elements('C16H15N3.99')
assert len(elements) == 3

# check molecular_weight and sum_formula_pDNA
# test if my two functions work as intended and compare results to the Bio.SeqUtils.molecular_weight function.
# the big difference: with my functions I also get the sum formula of the pDNA
sequence = 'ATCG'*100
test_sum_formula, test_sum_formula_dic = snek.elements.sum_formula_pDNA(sequence,return_dic=True)
test_mw_my_implementation = snek.elements.molecular_weight(test_sum_formula)
test_mw_bio = SeqUtils.molecular_weight(sequence,seq_type='DNA',double_stranded=True,circular=True,monoisotopic=False)
# assert that the error is smaller than 0.1
assert np.abs((test_mw_my_implementation-test_mw_bio)/np.mean([test_mw_bio,test_mw_my_implementation])*100) < .1


print('END OF TEST')
