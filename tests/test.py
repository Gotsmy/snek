import cobra
import snek

model = cobra.io.load_model("textbook")
solution1 = model.optimize()
solution2 = snek.sensitive_optimize(model)
assert all(solution1.fluxes == solution2.fluxes)

c_fluxes = snek.analysis.carbon_fluxes(model)

nr = snek.analysis.count_atom('C16H15N3.99','N')
assert nr == 3.99
all = snek.analysis.element_composition('C16H15N3.99')
print(all)
assert all == {'C':16,'H':15,'N':3.99}
elements = snek.analysis.unique_elements('C16H15N3.99')
print(elements)
