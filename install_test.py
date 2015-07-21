#!/usr/bin/python
import os

def run_test_case(test_case):
	print '****Running the ' + str(test_case) + ' case****'
	try:
		starting_dir = os.path.abspath(os.curdir)
		os.chdir(test_case)
		execfile('run.py')
		os.chdir(starting_dir)
		print '****The ' +  str(test_case) + ' test case ran successfully****\n'
	except:
		raise Exception('The ' +  str(test_case) + ' test case failed')

if __name__ == "__main__":

	run_test_case('tutorials/wamit/COER_hydrodynamic_modeling_comp')
	run_test_case('tutorials/wamit/ecm_ellipsoid')
	run_test_case('tutorials/wamit/oswec')
	run_test_case('tutorials/wamit/rm3')
	run_test_case('tutorials/wamit/sphere')
	run_test_case('tutorials/wamit/wec3')
	run_test_case('tutorials/nemoh')
	run_test_case('tutorials/mesh/scale_and_translate')
	run_test_case('tutorials/mesh/wamit_to_nemoh')

	print '\nAll test cases ran successfully!'
