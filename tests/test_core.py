import cobra
import snek

model = cobra.io.load_model("textbook")
model.solver = 'glpk'

# check core.snesitive_optimize
solution1 = model.optimize()
solution2 = snek.sensitive_optimize(model)
assert all(solution1.fluxes == solution2.fluxes)

snek.sensitive_optimize(model,pFBA=True)
print('End of test_core.py')
