A unifhy-compliant version of the rainfall-runoff model GR4
-----------------------------------------------------------

.. image:: https://img.shields.io/pypi/v/unifhycontrib-gr4?style=flat-square&color=00b0f0
   :target: https://pypi.python.org/pypi/unifhycontrib-gr4
   :alt: PyPI Version
.. image:: https://img.shields.io/badge/dynamic/json?url=https://zenodo.org/api/records/5780135&label=doi&query=doi&style=flat-square&color=00b0f0
   :target: https://zenodo.org/badge/latestdoi/355948261
   :alt: DOI
.. image:: https://img.shields.io/github/workflow/status/unifhy-org/unifhycontrib-gr4/Run%20tests?style=flat-square&label=tests
   :target: https://github.com/unifhy-org/unifhycontrib-gr4/actions/workflows/run_tests.yml
   :alt: Tests Status

The GR4 ("Génie Rural à 4 paramètres" [in French]) model is a
bucket-type rainfall-runoff model featuring four parameters.

It is typically used as a daily model, i.e. GR4J model
(`Perrin et al., 2003`_). It can also be used at other temporal resolutions,
e.g. hourly in GR4H model, provided an adjustment in its time-dependent
parameter and constant values is performed (`Ficchì et al., 2016`_). The model
has recently been expressed in a state-space formulation
(`Santos et al., 2018`_).

This version of the GR4 model is based on its explicit state-space
formulation and its recommended temporal resolutions are daily or hourly.
With either of these resolutions, time-dependent parameters *x2*, *x3*,
*x4*, and constant *nu* are expected for the daily case and they are
adjusted accordingly if temporal resolution is not daily.

The surface layer component of the GR4 model comprises the interception scheme.

The subsurface component of the GR4 model comprises the runoff generation
and the runoff routing using the Nash cascade.

The openwater component of the GR4 model comprises the runoff routing
using the routing store, and the inter-catchment groundwater exchange.

.. _`Perrin et al., 2003`: https://doi.org/10.1016/s0022-1694(03)00225-7
.. _`Ficchì et al., 2016`: https://doi.org/10.1016/j.jhydrol.2016.04.016
.. _`Santos et al., 2018`: https://doi.org/10.5194/gmd-11-1591-2018

:contributors: Thibault Hallouin [1,2]
:affiliations:
    1. Department of Meteorology, University of Reading
    2. National Centre for Atmospheric Science
:licence: GPL-2.0
:codebase: https://github.com/unifhy-org/unifhycontrib-gr4


How to install
~~~~~~~~~~~~~~

This package is available on PyPI, so you can simply run:

.. code-block:: bash

   pip install unifhycontrib-gr4
