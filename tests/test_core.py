import cobra
import snek

def test_sensitive_optimize():
    model = cobra.io.load_model("textbook")
    solution1 = model.optimize()
    solution2 = snek.sensitive_optimize(model)
    assert all(solution1.fluxes == solution2.fluxes)
    # check warning for sensitive_optimize pFBA + glpk solver
    # snek.sensitive_optimize(model,pFBA=True)

def test_get_solver():
    model = cobra.io.load_model("textbook")
    model.solver = 'glpk'
    assert snek.get_solver(model) == 'glpk'

def test_get_objetive():
    model = cobra.io.load_model("textbook")
    model.objective = 'EX_glc__D_e'
    assert snek.get_objective(model) == 'EX_glc__D_e'


# check snek.sensitive_optimize()
test_sensitive_optimize()
# check snek.get_solver()
test_get_solver()
# check snek.get_objective()
test_get_objetive()

print('End of test_core.py')
