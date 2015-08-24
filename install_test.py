#!/usr/bin/python
import os
import sys
# from bemio.__version__ import base

def run_test_case(test_case):

	print '\n****Running the ' + str(test_case) + ' case****'
	starting_dir = os.path.abspath(os.curdir)

	try:

		os.chdir(test_case)
		execfile('run.py')
		os.chdir(starting_dir)
		print '****The ' +  str(test_case) + ' test case ran successfully****\n'
		status = str(test_case) + ': SUCCESS'

	except:

		print '****The ' +  str(test_case) + ' test case failed. For the specific information on the error run the ' + str(test_case) + ' test case manually****\n'
		status = str(test_case) + ': FAILED'
		os.chdir(starting_dir)

	return status

if __name__ == "__main__":

	# print 'Testing bemio version ' + base()

	status = []

	if len(sys.argv) == 2 and sys.argv[1] == 'all':

		print 'Running all test cases:\n'
		status.append(run_test_case('tutorials/wamit/ecm_ellipsoid'))
		status.append(run_test_case('tutorials/wamit/sphere'))
		status.append(run_test_case('tutorials/wamit/oswec'))
		status.append(run_test_case('tutorials/wamit/rm3'))
		status.append(run_test_case('tutorials/mesh/scale_and_translate'))
		status.append(run_test_case('tutorials/mesh/wamit_to_nemoh'))
		status.append(run_test_case('tutorials/utilities/wave_excitation'))

	else:

		print 'Running standard test cases:\n'

	status.append(run_test_case('tutorials/wamit/COER_hydrodynamic_modeling_comp'))
	status.append(run_test_case('tutorials/wamit/wec3'))
	status.append(run_test_case('tutorials/nemoh'))
	status.append(run_test_case('tutorials/aqwa'))


	print 'Test cases runs completed:'
	for i, status_i in enumerate(status):
		print '\t' + status_i
