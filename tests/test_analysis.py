import cobra
import snek

model = cobra.io.load_model("textbook")

# check snek.analysis.carbon_fluxes
c_fluxes = snek.analysis.carbon_fluxes(model)
assert len(c_fluxes) == 3



print('End of test_analysis.py')
