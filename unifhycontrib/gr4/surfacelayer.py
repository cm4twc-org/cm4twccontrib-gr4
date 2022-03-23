import numpy as np

import unifhy


class SurfaceLayerComponent(unifhy.component.SurfaceLayerComponent):
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

    The surface layer component of the GR4 model comprises the
    interception scheme.

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
        'soil_water_stress_for_transpiration'
    }
    _outwards = {
        'canopy_liquid_throughfall_and_snow_melt_flux',
        'transpiration_flux_from_root_uptake'
    }
    _inputs_info = {
        'rainfall_flux': {
            'units': 'kg m-2 s-1',
            'kind': 'dynamic'
        },
        'potential_water_evapotranspiration_flux': {
            'units': 'kg m-2 s-1',
            'kind': 'dynamic'
        }
    }
    _parameters_info = {
        'x1': {
            'units': 'kg m-2'
        }
    }
    _constants_info = {
        'alpha': {
            'description': 'production precipitation exponent',
            'units': '1',
            'default_value': 2
        }
    }
    _outputs_info = {
        'actual_water_evapotranspiration_flux': {
            'units': 'kg m-2 s-1'
        }
    }

    def initialise(self,
                   # component states
                   **kwargs):
        pass

    def run(self,
            # from exchanger
            soil_water_stress_for_transpiration,
            # component inputs
            rainfall_flux, potential_water_evapotranspiration_flux,
            # component parameters
            x1,
            # component states
            # component constants
            alpha,
            **kwargs):

        # some name binding to be consistent with GR4J nomenclature
        dt = self.timedelta_in_seconds
        p = rainfall_flux * dt
        e = potential_water_evapotranspiration_flux * dt
        s_over_x1 = soil_water_stress_for_transpiration

        # interception
        e_minus_p = e - p
        en = np.maximum(e_minus_p, 0.0)
        pn = np.maximum(-e_minus_p, 0.0)

        # determine where water-/energy-limited conditions are
        water_limited = e_minus_p >= 0.0
        energy_limited = ~water_limited
        
        # --------------------------------------------------------------
        # under water-limited conditions (i.e. remaining energy 'en')
        # >-------------------------------------------------------------
        
        en_over_x1 = np.where(water_limited, en / x1, 0.0)
        # limited to 13, as per source code of Coron et al. (2017)
        # https://doi.org/10.1016/j.envsoft.2017.05.002
        en_over_x1[en_over_x1 > 13.0] = 13.0
        tanh_en_over_x1 = np.tanh(en_over_x1)

        s = np.where(water_limited, s_over_x1 * x1, 0.0)
        s_alpha_over_x1 = (s ** alpha) / x1

        es = np.where(
            water_limited,
            ((2 * s - s_alpha_over_x1) * tanh_en_over_x1)
             / (1 + (1 - s_over_x1) * tanh_en_over_x1),
            0.0
        )

        ae = np.where(water_limited, es + p, 0.0)

        # -------------------------------------------------------------<

        # --------------------------------------------------------------
        # under energy-limited conditions (i.e. remaining water 'pn')
        # >-------------------------------------------------------------

        ae[energy_limited] = e[energy_limited]

        # -------------------------------------------------------------<

        return (
            # to exchanger
            {
                'canopy_liquid_throughfall_and_snow_melt_flux':
                    pn / dt,
                'transpiration_flux_from_root_uptake':
                    es / dt
            },
            # component outputs
            {
                'actual_water_evapotranspiration_flux':
                    ae / dt
            }
        )

    def finalise(self, **kwargs):
        pass
