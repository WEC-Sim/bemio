import os

os.chdir('data')
os.system('rm -rf *body*.mat *.h5 *.p')
os.chdir('../')