#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2022--, katharoseq development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer

setup(
    name="q2-katharoseq",
    packages=find_packages(),
    version=versioneer.get_version(),
    #version='1.0.0',
    cmdclass=versioneer.get_cmdclass(),
    author="Daniela Perry",
    author_email="dsperry@ucsd.edu",
    description="katharo seq protocol for low biomass samles",
    license='BSD-3',
    entry_points={
        'qiime2.plugins': ['q2-katharoseq=q2_katharoseq.plugin_setup:plugin']
    },
    package_data={
        "q2_katharoseq": ['citations.bib'],
        "example": ['*']
    },
    zip_safe=False,
    install_requires=['scipy',
                      'scikit-learn',
                      'matplotlib',
                      'seaborn',
                      'pandas']
)
