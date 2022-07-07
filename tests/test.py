import cobra
import snek

model = cobra.io.load_model("textbook")
solution1 = model.optimize()
solution2 = snek.sensitive_optimize(model)
assert type(solution1) == type(solution2)
assert all(solution1.fluxes == solution2.fluxes)
