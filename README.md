## bemio
For more information and code documentation please visit the bemio website http://wec-sim.github.io/bemio for more.

## Easy Install
Download the source code ``git clone https://github.com/WEC-Sim/bemio.git`` and in the resulting bemio directory execute the following command ``python setup.py install --user``. This will copy the bemio directory into python as a new package. NOTE: Using this method, each time bemio is 'pulled' the code will need to be installed again.

## Alternate Install
Download the source code ``git clone https://github.com/WEC-Sim/bemio.git`` and add the bemio directory to the python path. 
- On windows this is done by following these instructions: http://stackoverflow.com/questions/3701646/how-to-add-to-the-pythonpath-in-windows-7
- On a MAC or LINUX, this is done by putting the following line in your ``~/.bash_rc`` or ``~/.bash_profile`` file: 
``export PYTHONPATH=$PYTHONPATH:/Users/your_username/Documents/WEC-SIM/git_global_bemio``
NOTE: If the 'Easy Install' was previously run, the bemio directory must also be removed from your python install 
