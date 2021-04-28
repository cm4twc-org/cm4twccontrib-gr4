from setuptools import setup, find_packages


with open("README.rst", 'r') as fh:
    long_desc = fh.read()

with open("gr4/version.py", 'r') as fv:
    exec(fv.read())


setup(
    name='cm4twccontrib.gr4',

    version=__version__,

    description='cm4twc components for GR4',
    long_description=long_desc,
    long_description_content_type="text/x-rst",

    author='Thibault Hallouin',

    project_urls={
        'Source Code': 'https://github.com/hydro-jules/cm4twccontrib-gr4'
    },

    license='GPL-2.0',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Hydrology',
    ],

    packages=find_packages(),

    install_requires=[
        'cm4twc'
    ],

)
