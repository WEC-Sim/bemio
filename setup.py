"""A setuptools based setup module for bemio
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

setup(

    name='bemio',

    version='1.0',

    description='bemio',

    long_description='Boundary Element Method Input/Output',

    url='https://github.com/WEC-Sim/bemio',

    author='Michael Lawson',

    author_email='michael.lawson@nrel.gov',

    license='Apache 2.0',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Wave Energy Research Community',
        'License ::Apache 2.0',
        'Programming Language :: Python :: 2.7'
    ],

    keywords='bemio',

    packages=find_packages(exclude=['doc', 'tutorials']),

    install_requires=['progressbar2'],

    extras_require={
        'dev': [],
        'test': [],
    },

)