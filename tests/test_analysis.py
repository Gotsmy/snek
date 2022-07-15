import os
import pytest
import numpy as np
import pandas as pd
import cobra
import snek

@pytest.fixture
def model():
    location = os.path.dirname(os.path.abspath(__file__))
    model = cobra.io.read_sbml_model(os.path.join(location,'data','textbook.xml'))
    return model

# check snek.analysis.in_flux
def test_in_flux(model):
    in_flux = snek.analysis.in_flux(model)
    assert in_flux['EX_glc__D_e'] == -10

# check snek.analysis.get_constrained_reactions
def test_get_constrained_reactions(model):
    df = snek.analysis.get_constrained_reactions(model)
    assert df.loc['ATPM','lower_bound'] == 8.39

# check snek.analysis.in_flux
def test_in_flux(model):
    df = snek.analysis.in_flux(model)
    assert df.loc['EX_glc__D_e','flux'] == -10

# check snek.analysis.out_flux
def test_out_flux(model):
    df = snek.analysis.in_flux(model)
    assert isinstance(df,type(pd.DataFrame()))

# check snek.analysis.in_bounds
def test_in_bounds(model):
    df = snek.analysis.in_bounds(model)
    assert df.loc['EX_h2o_e','upper_bound'] == 1000

# check snek.analysis.out_bounds
def test_out_bounds(model):
    df = snek.analysis.out_bounds(model)
    assert df.loc['EX_h2o_e','upper_bound'] == 1000

print('End of test_analysis.py')
