#! /usr/bin/python
import shutil
try:
     shutil.rmtree('build')
except:
    pass
try:
    shutil.rmtree('bemio.egg-info')
except:
    pass
try:
    shutil.rmtree('dist')
except:
    pass
