import cobra
import snek

model = cobra.io.load_model("textbook")
model.solver = 'glpk'

# check core.sensitive_optimize
solution1 = model.optimize()
solution2 = snek.sensitive_optimize(model)
assert all(solution1.fluxes == solution2.fluxes)
# check warning for sensitive_optimize pFBA + glpk solver
snek.sensitive_optimize(model,pFBA=True)

# check get_solver()
model.solver = 'glpk'
assert snek.get_solver(model) == 'glpk'

print('End of test_core.py')
