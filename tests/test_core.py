import os
import pytest
import cobra
import snek

@pytest.fixture
def model():
    location = os.path.dirname(os.path.abspath(__file__))
    model = cobra.io.read_sbml_model(os.path.join(location,'data','textbook.xml'))
    return model

# check snek.set_objective
def test_set_objective(model):
    snek.set_objective(model,'Biomass_Ecoli_core','max')
    assert 'ATPM' not in str(model.objective.expression)
    snek.set_objective(model,'ATPM','max')
    assert 'ATPM' in str(model.objective.expression)

# check snek.set_bounds
def test_set_bounds(model):
    snek.set_bounds(model,'EX_glc__D_e',lb=-5,ub=3)
    assert model.reactions.get_by_id('EX_glc__D_e').bounds == (-5,3)

# check snek.sensitive_optimize()
def test_sensitive_optimize(model):
    solution1 = model.optimize()
    solution2 = snek.sensitive_optimize(model)
    assert all(solution1.fluxes == solution2.fluxes)
    # check warning for sensitive_optimize pFBA + glpk solver
    # snek.sensitive_optimize(model,pFBA=True)

# check snek.find_biomass_reaction()
def test_find_biomass_reaction(model):
    df = snek.find_biomass_reaction(model)
    assert df.loc[0,'id'] == 'Biomass_Ecoli_core'

# check find_reaction
def test_find_reaction(model):
    df = snek.find_reaction(model,'su')
    assert df.shape == (5, 2)
    assert list(df['id']) == ['EX_succ_e', 'SUCCt2_2', 'SUCCt3', 'SUCDi', 'SUCOAS']

# check find_metabolite
def test_find_metabolite(model):
    df = snek.find_metabolite(model,'coa')
    assert df.shape == (3, 2)
    assert list(df['id']) == ['accoa_c', 'coa_c', 'succoa_c']

# check snek.get_objective()
def test_get_objetive(model):
    model.objective = 'EX_glc__D_e'
    assert snek.get_objective(model) == 'EX_glc__D_e'

# check snek.get_solver()
def test_get_solver(model):
    model.solver = 'glpk'
    assert snek.get_solver(model) == 'glpk'

print('End of test_core.py')
