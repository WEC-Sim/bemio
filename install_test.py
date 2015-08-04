#!/usr/bin/python
import os
from bemio.__version__ import base

def run_test_case(test_case):
	print '****Running the ' + str(test_case) + ' case****'
	try:
		starting_dir = os.path.abspath(os.curdir)
		os.chdir(test_case)
		execfile('run.py')
		os.chdir(starting_dir)
		print '****The ' +  str(test_case) + ' test case ran successfully****\n'
		status = str(test_case) + ': SUCCESS'
	except:
		print 'The ' +  str(test_case) + ' test case failed'
		status = status = str(test_case) + ': *FAILED*'
	return status

if __name__ == "__main__":
	print 'Testing bemio version ' + base() + '\n'
	status = []
	status.append(run_test_case('tutorials/wamit/COER_hydrodynamic_modeling_comp'))
	status.append(run_test_case('tutorials/wamit/ecm_ellipsoid'))
	status.append(run_test_case('tutorials/wamit/oswec'))
	status.append(run_test_case('tutorials/wamit/rm3'))
	status.append(run_test_case('tutorials/wamit/sphere'))
	status.append(run_test_case('tutorials/wamit/wec3'))
	status.append(run_test_case('tutorials/nemoh'))
	status.append(run_test_case('tutorials/mesh/scale_and_translate'))
	status.append(run_test_case('tutorials/mesh/wamit_to_nemoh'))
	status.append(run_test_case('tutorials/utilities/wave_excitation'))

	print 'All test cases run:'
	for i, status_i in enumerate(status):
		print '\t' + status_i
