import os
import pytest
import cobra
import snek

@pytest.fixture
def model():
    location = os.path.dirname(os.path.abspath(__file__))
    model = cobra.io.read_sbml_model(os.path.join(location,'data','textbook.xml'))
    return model

# check snek.sensitive_optimize()
def test_sensitive_optimize(model):
    solution1 = model.optimize()
    solution2 = snek.sensitive_optimize(model)
    assert all(solution1.fluxes == solution2.fluxes)
    # check warning for sensitive_optimize pFBA + glpk solver
    # snek.sensitive_optimize(model,pFBA=True)

# check snek.get_solver()
def test_get_solver(model):
    model.solver = 'glpk'
    assert snek.get_solver(model) == 'glpk'

# check snek.get_objective()
def test_get_objetive(model):
    model.objective = 'EX_glc__D_e'
    assert snek.get_objective(model) == 'EX_glc__D_e'

print('End of test_core.py')
