import numpy as np

from cm4twc.components import SubSurfaceComponent
from cm4twc.settings import dtype_float


class GR4(SubSurfaceComponent):
    """
    The GR4 ("Génie Rural à 4 paramètres" [in French]) model is a
    bucket-style rainfall-runoff model featuring four parameters. It is
    typically used as a daily model, i.e. GR4J model (`Perrin et al., 2003`_)
    where J stands for "journalier", meaning daily in French. It can also
    be used at other temporal resolutions, e.g. hourly, provided an adjustment
    in its parameter values is performed (`Ficchì et al., 2016`_). The
    model has recently been expressed in a state-space formulation
    (`Santos et al., 2018`_).

    This version of the GR4 model is based on its explicit state-space
    formulation and it can be used at any temporal resolution provided the
    parameters featuring 'timedelta' in their units are adjusted accordingly
    (see `Ficchì et al., 2016`_).

    The subsurface component of the GR4 model comprises the runoff
    generation and runoff routing processes.

    .. _`Perrin et al., 2003`: https://doi.org/10.1016/s0022-1694(03)00225-7
    .. _`Ficchì et al., 2016`: https://doi.org/10.1016/j.jhydrol.2016.04.016
    .. _`Santos et al., 2018`: https://doi.org/10.5194/gmd-11-1591-2018

    :contributors: Thibault Hallouin [1,2]
    :affiliations:
        1. Department of Meteorology, University of Reading
        2. National Centre for Atmospheric Science
    :licence: GPL-2.0
    """

    _parameters_info = {
        'x1': {
            'units': 'kg m-2'
        },
        'x4': {
            'units': 'timedelta'
        }
    }
    _states_info = {
        'production_store': {
            'units': 'kg m-2'
        },
        'nash_cascade_stores': {
            'description': 'number of stores in Nash cascade (nres)',
            'units': 'kg m-2',
            'divisions': 11
        }
    }
    _constants_info = {
        'alpha': {
            'description': 'production precipitation exponent',
            'units': '1'
        },
        'beta': {
            'description': 'percolation exponent',
            'units': '1'
        },
        'nu': {
            'description': 'percolation coefficient',
            'units': '1'
        }
    }

    def initialise(self,
                   # component states
                   production_store, nash_cascade_stores,
                   **kwargs):
        production_store[-1][:] = 0.0
        nash_cascade_stores[-1][:] = 0.0

    def run(self,
            # from exchanger
            transpiration, evaporation_soil_surface, evaporation_ponded_water,
            throughfall, snowmelt,
            # component inputs
            # component parameters
            x1, x4,
            # component states
            production_store, nash_cascade_stores,
            # component constants
            alpha=2, beta=5, nu=4/9,
            **kwargs):

        # some name binding to be consistent with GR4J nomenclature
        dt = self.timedelta_in_seconds
        pn = (throughfall + snowmelt) * dt
        es = (transpiration + evaporation_soil_surface
              + evaporation_ponded_water) * dt
        s_ = production_store[-1]
        sh_ = nash_cascade_stores[-1]
        nres = sh_.shape[-1]

        # determine where energy limited conditions are
        energy_limited = pn > 0.0

        # energy-limited conditions (i.e. remaining water 'pn')
        s_over_x1 = np.where(energy_limited, s_ / x1, 0.0)
        pn_over_x1 = np.where(energy_limited, pn / x1, 0.0)
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

        # update production store after infiltration and evaporation
        s = s_ + ps - es

        # percolation from production store
        s *= s > 0
        s_over_x1 = s / x1

        # nu is time dependent - see Ficchi et al. (2016)
        # https://doi.org/10.1016/j.jhydrol.2016.04.016
        nu = nu * (86400 / dt) ** 0.25

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

        outflow_coefficient = (1 - np.exp((1 - nres)/x4))

        qsh[...] = sh_[...] * outflow_coefficient
        sh[..., 0] = sh_[..., 0] + pr - qsh[..., 0]
        sh[..., 1:nres] = (
            sh_[..., 1:nres] + qsh[..., 0:nres-1] - qsh[..., 1:nres]
        )

        quh = qsh[..., -1]

        # update component states
        production_store[0][:] = s
        nash_cascade_stores[0][:] = sh

        return (
            # to exchanger
            {
                'surface_runoff': 
                    np.zeros(self.spaceshape, dtype_float()),
                'subsurface_runoff': 
                    quh / dt,
                'soil_water_stress': 
                    s / x1
            },
            # component outputs
            {}
        )

    def finalise(self, **kwargs):
        pass
