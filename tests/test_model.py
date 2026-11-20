import pytest
import numpy as np
from simulation.model import ModelParameters, InflammatoryOncoContinuum

def test_model_initialization():
    params = ModelParameters()
    assert params.r_v == 10.0
    assert params.r_ox == 200.0
    assert params.lambda_C > 0

    model = InflammatoryOncoContinuum(params)
    assert len(model.r) == params.Nr
    assert model.r[0] == pytest.approx(10.0)
    assert model.r[-1] == pytest.approx(200.0)

def test_bessel_k0_decay():
    params = ModelParameters()
    model = InflammatoryOncoContinuum(params)
    C = model.steady_state_drug_profile(C_plasma=params.C_plasma_on)
    assert len(C) == len(model.r)
    # Radial drug decay should be monotonically decreasing
    assert np.all(np.diff(C) <= 0)
    assert C[0] > C[-1]

def test_sanctuary_volume_fraction():
    model = InflammatoryOncoContinuum(ModelParameters())
    assert 0.0 < model.sanctuary_volume_fraction < 1.0
    assert 0.25 < model.sanctuary_volume_fraction < 0.40

def test_short_simulation_run():
    params = ModelParameters(Nr=30)
    model = InflammatoryOncoContinuum(params)
    res = model.simulate(t_span_months=2.0, t_therapy_stop=1.0, dt_months=0.1)
    assert "time" in res
    assert "N_tot" in res
    assert len(res["time"]) > 0
    assert np.all(res["N_tot"] >= 0.0)
