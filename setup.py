#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2022--, katharoseq development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="q2-katharoseq",
    packages=find_packages(),
    version='0.0.1',
    author="Daniela Perry",
    author_email="dsperry@ucsd.edu",
    description="katharo seq protocol for low biomass samles",
    license='BSD-3',
    entry_points={
        'qiime2.plugins': ['q2-katharoseq=q2_katharoseq.plugin_setup:plugin']
    },
    package_data={
        "q2_katharoseq": ['citations.bib'],
    },
    zip_safe=False,
    install_requires=['scipy',
                      'matplotlib',
                      'seaborn',
                      'pandas']
)
