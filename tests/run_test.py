import unittest
from datetime import datetime, timedelta
import cm4twc

from cm4twccontrib.gr4 import (
    SurfaceLayerComponent, SubSurfaceComponent, OpenWaterComponent
)


class TestContribution(unittest.TestCase):

    def test_gr4j(self):
        td = cm4twc.TimeDomain.from_start_end_step(
            start=datetime(2001, 1, 1, 0, 0, 0),
            end=datetime(2002, 1, 1, 0, 0, 0),
            step=timedelta(days=1)
        )

        sd = cm4twc.LatLonGrid.from_extent_and_resolution(
            latitude_extent=(51, 52),
            latitude_resolution=1,
            longitude_extent=(0, 1),
            longitude_resolution=1
        )

        ds = cm4twc.DataSet(['in/rainfall_flux.nc',
                             'in/potential_water_evapotranspiration_flux.nc'])

        area = 1844  # km2

        x1 = (346.9290884, 'kg m-2')
        x2 = (-0.0458, 'kg m-2 d-1')
        x3 = (119.780094, 'kg m-2')
        x4 = (0.9384945, 'd')

        sl = SurfaceLayerComponent(
            saving_directory='out',
            timedomain=td,
            spacedomain=sd,
            dataset=ds,
            parameters={
                'x1': x1
            }
        )

        ss = SubSurfaceComponent(
            saving_directory='out',
            timedomain=td,
            spacedomain=sd,
            dataset=None,
            parameters={
                'x1': x1,
                'x4': x4
            }
        )

        ow = OpenWaterComponent(
            saving_directory='out',
            timedomain=td,
            spacedomain=sd,
            dataset=None,
            parameters={
                'x2': x2,
                'x3': x3
            },
            records={
                'outgoing_water_volume_transport_along_river_channel': {
                    timedelta(days=1): ['mean']
                }
            }
        )

        model = cm4twc.Model(
            identifier='test-gr4j',
            config_directory='out',
            saving_directory='out',
            surfacelayer=sl,
            subsurface=ss,
            openwater=ow
        )

        model.simulate()

        from_file = cm4twc.DataSet(
            'in/outgoing_water_volume_transport_along_river_channel.nc'
        )

        from_model = cm4twc.DataSet(
            'out/test-gr4j_openwater_run_records_daily.nc'
        )

        var_name = 'outgoing_water_volume_transport_along_river_channel'

        self.assertTrue(
            from_file[var_name].field.equals(from_model[var_name].field,
                                             verbose=3)
        )


if __name__ == '__main__':
    test_loader = unittest.TestLoader()
    test_suite = unittest.TestSuite()

    test_suite.addTests(
        test_loader.loadTestsFromTestCase(TestContribution))

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(test_suite)
