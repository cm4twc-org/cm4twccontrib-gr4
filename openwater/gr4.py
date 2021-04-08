import numpy as np

from cm4twc.components import OpenWaterComponent


class GR4(OpenWaterComponent):
    """
    The GR4 ("Génie Rural à 4 paramètres" [in French]) model is a
    bucket-style rainfall-runoff model featuring four parameters. It is
    typically used as a daily model, i.e. GR4J model (`Perrin et al., 2003`_)
    where J stands for "journalier", meaning daily in French. It can also
    be used at other temporal resolutions, e.g. hourly, provided an adjustment
    in its parameter values is performed (`Ficchì et al., 2016`_). The
    model has recently been expressed in a state-space formulation
    (`Santos et al., 2018`_).

    This version of the GR4 model is based on its state-space formulation
    and it can be used at any spatial resolution provided the parameters
    featuring 'timedelta' in their units are adjusted accordingly
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
        'x2': {
            'units': 'kg m-2 timedelta-1'
        },
        'x3': {
            'units': 'kg m-2'
        }
    }
    _states_info = {
        'routing_store': {
            'units': 'kg m-2'
        }
    }
    _constants_info = {
        'gamma': {
            'description': 'routing outflow exponent',
            'units': '1'
        },
        'omega': {
            'description': 'exchange exponent',
            'units': '1'
        },
        'phi': {
            'description': 'partition between routing store and direct flow',
            'units': '1'
        }
    }
    _outputs_info = {
        'outgoing_water_volume_transport_along_river_channel': {
            'units': 'm3 s-1',
            'description': 'streamflow at outlet'
        }
    }

    def initialise(self,
                   # component states
                   routing_store,
                   **kwargs):
        routing_store[-1][:] = 0.0

    def run(self,
            # from exchanger
            surface_runoff, subsurface_runoff, evaporation_openwater,
            # component inputs
            # component parameters
            x2, x3,
            # component states
            routing_store,
            # component constants
            gamma=5, omega=3.5, phi=0.9,
            **kwargs):

        # some name binding to be consistent with GR4J nomenclature
        quh = (surface_runoff + subsurface_runoff) * self.timedelta_in_seconds
        r_ = routing_store[-1]

        # split runoff between direct and routing
        q9 = quh * phi
        q1 = quh * (1 - phi)

        # potential inter-catchment exchange
        with np.errstate(over='ignore'):
            f = np.where(
                r_ > 0.0,
                x2 * (r_ / x3) ** omega,
                0.0
            )

        # runoff from routing store
        r = r_ + q9 + f
        r[r < 0.0] = 0.0
        with np.errstate(divide='ignore', invalid='ignore'):
            qr = np.where(
                r > 0.0,
                r * (1 - (1 + ((r / x3) ** (gamma - 1)) ** (1 / (1 - gamma)))),
                0.0
            )

        # runoff from direct branch
        qd = np.maximum(q1 + f, 0.0)

        # total runoff
        q = qr + qd
        q[q < 0.0] = 0.0

        # update component states
        routing_store[0][:] = r

        return (
            # to exchanger
            {
                'water_level': r
            },
            # component outputs
            {
                'outgoing_water_volume_transport_along_river_channel': q
            }
        )

    def finalise(self, **kwargs):
        pass
