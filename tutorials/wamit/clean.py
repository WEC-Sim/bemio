#! /usr/bin/python

import os

os.chdir('data')
os.system('rm *.mat *.h5 *.p')
os.chdir('../')