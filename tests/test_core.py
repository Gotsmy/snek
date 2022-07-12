import cobra
import snek

model = cobra.io.load_model("textbook")

# check core.snesitive_optimize
solution1 = model.optimize()
solution2 = snek.sensitive_optimize(model)
assert all(solution1.fluxes == solution2.fluxes)


print('End of test_core.py')
