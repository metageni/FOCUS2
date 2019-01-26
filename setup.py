#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

requirements = [
    'numpy >= 1.12.1',
    'scipy >= 0.19.0',
    'pysam'
]

test_requirements = [
    'pytest'
]

setup_requirements = []

setup(
    name='metagenomics-focus2',
    use_scm_version=True,
    description='FOCUS2: Agile and sensitive classification of metagenomics data using a reduced database',
    author='Genivaldo G. Z. Silva',
    author_email='genivaldo.gueiros@gmail.com',
    url='https://github.com/metageni/focus2',
    packages=[
        'focus2_app',
    ],
    package_dir={'focus2_app': 'focus2_app'},
    include_package_data=True,
    install_requires=requirements,
    setup_requires=setup_requirements,
    zip_safe=False,
    keywords='focus2_app',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    entry_points={
        'console_scripts': [
            'focus2 = focus2_app.focus2:main',
        ]
    },
)
