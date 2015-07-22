#!/usr/bin/python

import os

def run_test_case(test):
	print '****Running the ' + str(test) + ' case****'
	try:
		os.chdir(test)
		execfile('run.py')
		os.chdir('..')
		print '****The ' +  str(test) + ' test case ran successfully****\n'
	except:
		raise Exception('The ' +  str(test) + ' test case did not run successfully') 
	

if __name__ == "__main__":

	os.chdir('tutorials')

	run_test_case('wamit')
	run_test_case('aqwa')
	run_test_case('nemoh')
	#run_test_case('mesh') # for code developers only
	
	os.chdir('..')