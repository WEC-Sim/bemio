"""A setuptools based setup module for bemio
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

from bemio.__version__ import base

setup(

    name='bemio',

    version=base(),

    description='bemio',

    long_description='Boundary Element Method Input/Output',

    url='https://github.com/WEC-Sim/bemio',

    author='Michael Lawson, Yi-Hsiang Yu, Carlos Michelen',

    author_email='michael.lawson@nrel.gov',

    license='Apache 2.0',

    classifiers=[
        'Development Status :: Release',
        'Intended Audience :: Ocean research community',
        'License ::Apache 2.0',
        'Programming Language :: Python :: 2.7'
    ],

    keywords='bemio',

    packages=find_packages(exclude=['doc', 'tutorials']),

    install_requires=['numpy', 'scipy', 'h5py', 'progressbar2', 'astropy'],

    extras_require={
        'dev': [],
        'test': [],
    },

)
