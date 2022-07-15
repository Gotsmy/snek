import cobra
import snek

def test_carbon_fluxes():
    model = cobra.io.load_model("textbook")
    c_fluxes = snek.analysis.carbon_fluxes(model)
    assert len(c_fluxes) == 3

def test_in_flux():
    model = cobra.io.load_model("textbook")
    in_flux = snek.analysis.in_flux(model)
    assert in_flux['EX_glc__D_e'] == -10

# check snek.analysis.carbon_fluxes
test_carbon_fluxes()
# check snek.analysis.in_flux
test_in_flux()

print('End of test_analysis.py')
