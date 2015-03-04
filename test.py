import os

def run_test_case(test):
	
	try:
		os.chdir(test)
		execfile('run.py')
		os.chdir('..')
		print '\n****The ' +  str(test) + ' test case ran successfully****\n'
	except:
		raise Exception('The ' +  str(test) + ' test case did not run successfully') 
	

if __name__ == "__main__":

	os.chdir('tutorials')

	run_test_case('wamit')
	run_test_case('aqwa')
	run_test_case('mesh')
	# run_test_case('aqwa') # Uncomment when the aqwaio script is updated and completed	
	
	os.chdir('..')