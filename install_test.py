#!/usr/bin/python
import os

def run_test_case(test_case):
	print '****Running the ' + str(test_case) + ' case****'
	try:
		os.chdir(test_case)
		execfile('run.py')
		os.chdir('..')
		print '****The ' +  str(test_case) + ' test case ran successfully****\n'
	except:
		raise Exception('The ' +  str(test_case) + ' test case failed')

if __name__ == "__main__":

	starting_dir = os.path.abspath(os.curdir)
	os.chdir('tutorials')

	run_test_case('wamit')
	run_test_case('aqwa')
	run_test_case('nemoh')

	try:

		import vtk
		run_mesh = True

	except:

		print 'Not testing the mesh functionality because VTK is not installed'

	if run_mesh is True:

		run_test_case('mesh')

	os.chdir(starting_dir)
