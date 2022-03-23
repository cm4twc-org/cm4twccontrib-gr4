import numpy as np

import unifhy


class OpenWaterComponent(unifhy.component.OpenWaterComponent):
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

    The openwater component of the GR4 model comprises the runoff
    routing using the routing store, and the inter-catchment groundwater
    exchange.

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
        'surface_runoff_flux_delivered_to_rivers'
    }
    _outwards = {}
    _parameters_info = {
        'x2': {
            'units': 'kg m-2 d-1'
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
            'units': '1',
            'default_value': 5
        },
        'omega': {
            'description': 'exchange exponent',
            'units': '1',
            'default_value': 3.5
        },
        'phi': {
            'description': 'partition between routing store and direct flow',
            'units': '1',
            'default_value': 0.9
        }
    }
    _outputs_info = {
        'outgoing_water_volume_transport_along_river_channel': {
            'units': 'kg m-2 s-1',
            'description': 'streamflow at outlet'
        }
    }

    def initialise(self,
                   # component states
                   routing_store,
                   **kwargs):

        if not self.initialised_states:
            routing_store.set_timestep(-1, 0.0)

    def run(self,
            # from exchanger
            surface_runoff_flux_delivered_to_rivers,
            # component inputs
            # component parameters
            x2, x3,
            # component states
            routing_store,
            # component constants
            gamma, omega, phi,
            **kwargs):

        # some name binding to be consistent with GR4J nomenclature
        dt = self.timedelta_in_seconds
        quh = surface_runoff_flux_delivered_to_rivers * dt
        r_ = routing_store.get_timestep(-1)

        # convert time dependent parameters and constants
        # Ficchì et al. (2016) https://doi.org/10.1016/j.jhydrol.2016.04.016
        x2 = x2 * (86400 / dt) ** -0.125
        x3 = x3 * (86400 / dt) ** 0.25

        # split runoff between direct and routing
        q9 = quh * phi
        q1 = quh - q9

        # potential inter-catchment exchange
        r_over_x3 = r_ / x3
        f = x2 * r_over_x3 ** omega

        # runoff from routing store
        r = r_ + q9 + f
        r_over_x3 = r / x3
        qr = r * (1 - (1 + r_over_x3 ** (gamma - 1)) ** (1 / (1 - gamma)))

        r -= qr
        r *= r > 0

        # runoff from direct branch
        qd = np.maximum(q1 + f, 0.0)

        # total runoff
        q = qr + qd
        q *= q > 0

        # update component states
        routing_store.set_timestep(0, r)

        return (
            # to exchanger
            {},
            # component outputs
            {
                'outgoing_water_volume_transport_along_river_channel':
                    q / dt
            }
        )

    def finalise(self, **kwargs):
        pass
