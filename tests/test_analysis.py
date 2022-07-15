import pytest
import cobra
import snek

@pytest.fixture
def model():
    model = cobra.io.load_model("textbook")
    return model

# check snek.analysis.in_flux
def test_in_flux(model):
    in_flux = snek.analysis.in_flux(model)
    assert in_flux['EX_glc__D_e'] == -10

print('End of test_analysis.py')
