import cobra
import snek

model = cobra.io.load_model("textbook")

# check snek.analysis.carbon_fluxes
c_fluxes = snek.analysis.carbon_fluxes(model)
assert len(c_fluxes) == 3

# check snek.analysis.in_flux
in_flux = snek.analysis.in_flux(model)
assert in_flux['EX_glc__D_e'] == -10

print('End of test_analysis.py')
