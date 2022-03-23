import numpy as np

import unifhy


class SubSurfaceComponent(unifhy.component.SubSurfaceComponent):
    """The GR4 ("Génie Rural à 4 paramètres" [in French]) model is a
    bucket-type rainfall-runoff model featuring four parameters.

    It is typically used as a daily model, i.e. GR4J model
    (`Perrin et al., 2003`_). It can also be used at other temporal
    resolutions, e.g. hourly in GR4H model, provided an adjustment in
    its time-dependent parameter and constant values is performed
    (`Ficchì et al., 2016`_). The model has recently been expressed in
    a state-space formulation (`Santos et al., 2018`_).

    This version of the GR4 model is based on its explicit state-space
    formulation and its recommended temporal resolutions are daily or hourly.
    With either of these resolutions, time-dependent parameters *x2*, *x3*,
    *x4*, and constant *nu* are expected for the daily case and they are
    adjusted accordingly if temporal resolution is not daily.

    The subsurface component of the GR4 model comprises the runoff
    generation and the runoff routing using the Nash cascade.

    .. _`Perrin et al., 2003`: https://doi.org/10.1016/s0022-1694(03)00225-7
    .. _`Ficchì et al., 2016`: https://doi.org/10.1016/j.jhydrol.2016.04.016
    .. _`Santos et al., 2018`: https://doi.org/10.5194/gmd-11-1591-2018

    :contributors: Thibault Hallouin [1,2]
    :affiliations:
        1. Department of Meteorology, University of Reading
        2. National Centre for Atmospheric Science
    :licence: GPL-2.0
    :codebase: https://github.com/unifhy-org/unifhycontrib-gr4
    """

    _inwards = {
        'canopy_liquid_throughfall_and_snow_melt_flux',
        'transpiration_flux_from_root_uptake'
    }
    _outwards = {
        'surface_runoff_flux_delivered_to_rivers',
        'soil_water_stress_for_transpiration'
    }
    _parameters_info = {
        'x1': {
            'units': 'kg m-2'
        },
        'x4': {
            'units': 'd'
        }
    }
    _states_info = {
        'production_store': {
            'units': 'kg m-2'
        },
        'nash_cascade_stores': {
            'description': 'stores in Nash cascade',
            'units': 'kg m-2',
            'divisions': 'nres'
        }
    }
    _constants_info = {
        'alpha': {
            'description': 'production precipitation exponent',
            'units': '1',
            'default_value': 2
        },
        'beta': {
            'description': 'percolation exponent',
            'units': '1',
            'default_value': 5
        },
        'nu': {
            'description': 'percolation coefficient',
            'units': '1',
            'default_value': 4/9
        },
        'nres': {
            'description': 'number of stores in Nash cascade',
            'units': '1',
            'default_value': 11
        }
    }

    def initialise(self,
                   # component states
                   production_store, nash_cascade_stores,
                   **kwargs):

        if not self.initialised_states:
            production_store.set_timestep(-1, 0.0)
            nash_cascade_stores.set_timestep(-1, 0.0)

    def run(self,
            # from exchanger
            canopy_liquid_throughfall_and_snow_melt_flux,
            transpiration_flux_from_root_uptake,
            # component inputs
            # component parameters
            x1, x4,
            # component states
            production_store, nash_cascade_stores,
            # component constants
            alpha, beta, nu, nres,
            **kwargs):

        # some name binding to be consistent with GR4J nomenclature
        dt = self.timedelta_in_seconds
        pn = canopy_liquid_throughfall_and_snow_melt_flux * dt
        es = transpiration_flux_from_root_uptake * dt
        s_ = production_store.get_timestep(-1)
        sh_ = nash_cascade_stores.get_timestep(-1)

        # convert time dependent parameters and constants
        # Ficchì et al. (2016) https://doi.org/10.1016/j.jhydrol.2016.04.016
        nu = nu * (86400 / dt) ** 0.25
        x4 = x4 * (86400 / dt)

        # determine where energy limited conditions are
        energy_limited = pn > 0.0

        # --------------------------------------------------------------
        # under energy-limited conditions (i.e. remaining water 'pn')
        # >-------------------------------------------------------------

        s_over_x1 = s_ / x1
        pn_over_x1 = pn / x1
        # limited to 13, as per source code of Coron et al. (2017)
        # https://doi.org/10.1016/j.envsoft.2017.05.002
        pn_over_x1[pn_over_x1 > 13.0] = 13.0
        tanh_pn_over_x1 = np.tanh(pn_over_x1)

        ps = np.where(
            energy_limited,
            (x1 * (1 - s_over_x1 ** alpha) * tanh_pn_over_x1)
            / (1 + s_over_x1 * tanh_pn_over_x1),
            0.0
        )

        pr = pn - ps

        # -------------------------------------------------------------<

        # update production store after infiltration and evaporation
        s = s_ + ps - es

        # percolation from production store
        s *= s > 0
        s_over_x1 = s / x1

        perc = (
            s * (1 - (1 + ((nu * s_over_x1) ** (beta - 1)))
                 ** (1 / (1 - beta)))
        )

        pr += perc

        # update production store after percolation
        s -= perc

        # routing through nash cascade (nres stores in series)
        sh = np.zeros(sh_.shape)
        qsh = np.zeros(sh_.shape)

        outflow_coefficient = 1 - np.exp((1 - nres) / x4)

        sh[..., 0] = sh_[..., 0] + pr
        qsh[..., 0] = sh[..., 0] * outflow_coefficient
        sh[..., 0] -= qsh[..., 0]
        for i in range(1, nres):
            sh[..., i] = sh_[..., i] + qsh[..., i - 1]
            qsh[..., i] = sh[..., i] * outflow_coefficient
            sh[..., i] -= qsh[..., i]

        quh = qsh[..., -1]

        # update component states
        production_store.set_timestep(0, s)
        nash_cascade_stores.set_timestep(0, sh)

        return (
            # to exchanger
            {
                'surface_runoff_flux_delivered_to_rivers':
                    quh / dt,
                'soil_water_stress_for_transpiration':
                    s / x1
            },
            # component outputs
            {}
        )

    def finalise(self, **kwargs):
        pass
