.. _Installing VTK with Python Bindings:

Installing VTK with Python Bindings
===================================
In order to use ``bemio.mesh_utilities`` packaged with bemio you must have the Visualization Toolkit (VTK)  and the VTK Python 2.7.x bindings installed on your system.

.. Note::
    Installing the VTK Python bindings is only required if you would like to use bemio's :mod:`bemio.mesh_utilities` features.

.. Note::
	Some Python distributions, such as Anaconda and PythonXY, may include VTK bindings, however, it is recommended that bemio users install VTK using the instructions in this section because the VTK versions distributed with Anaconda and PythonXY are not kept current.

OSX
---
Method 1: Simple
~~~~~~~~~~~~~~~~
Install VTK using `MacPorts <https://www.macports.org/>`_, `Fink <http://www.finkproject.org/>`_, or `HomeBrew <http://brew.sh/>`_. If you would like to install VTK using this method you should also install Python and other bemio dependencies using this method.

Method 2: Complex
~~~~~~~~~~~~~~~~~
Pre-Installation Tasks
.......................
**Step 1:** Insure the following software is installed on your system:

* `Python-2.7.x <https://www.python.org/downloads/>`_
* `Apple Xcode <https://developer.apple.com/xcode/downloads/>`_
* `cmake <http://www.cmake.org/>`_
* `tcl <http://www.tcl.tk/>`_
* `tk <http://www.tcl.tk/>`_

Python is preloaded on OSX systems and XCode can be installed through the `OSX App Store <https://itunes.apple.com/us/app/xcode/id497799835;jsessionid=wnox0jq0k0vj2wbtcl3kohuf?mt=12#>`_ . The other packages can easily be installed using `MacPorts <https://www.macports.org/>`_, `Fink <http://www.finkproject.org/>`_, or `HomeBrew <http://brew.sh/>`_. Here is an example of how to install cmake, tck, and tk using MacPorts:

.. code-block:: shell

    sudo port install cmake
    sudo port install tcl
    sudo port install tk

**Step 2**: Download and unzip the latest version of VTK from here - `http://www.vtk.org/download <http://www.vtk.org/download/>`_.

.. Note::

    The folder that contains the VTK source code will be referred to as ``$VTK_SOURCE``

Install
............
Open a terminal window and use the following procedure to compile VTK and the Python bindings.

.. code-block:: shell
    :emphasize-lines: 17, 35, 39, 46

    # Move to the VTK source directory
    mlawson@mbp:~$ cd $VTK_SOURCE

    # Make a directory for the build and move into the build directory
    mlawson@mbp:VTK-6.2.0$ mkdir build
    mlawson@mbp:VTK-6.2.0$ cd build

    # Run cmake in interactive mode. Insure that the ``PYTHON_VERSION`` is set to ``2``, the
    # ``VTK_WRAP_PYTHON`` is set to ``ON, and ``BUILD_SHARED_LIBS`` is set to ON. ``BUILD_SHARED_LIBS``
    # is not always visable and you may need to enter the "advnacec mode" by pressing ``t`` to confirm this option is
    # correctly set. After the options are properly set, press ``c`` to configure and then ``g``
    # to generate the required make files.
    mlawson@mbp:VTK-6.2.0$ sudo ccmake ..

        BUILD_DOCUMENTATION              OFF
        BUILD_EXAMPLES                   OFF
        BUILD_SHARED_LIBS                ON       # May need to toggle advanced mode to view this option
        BUILD_TESTING                    OFF
        BUILD_USER_DEFINED_LIBS          OFF
        CMAKE_BUILD_TYPE                 Release
        CMAKE_FRAMEWORK_INSTALL_PREFIX   /usr/local/frameworks
        IOS_DEVICE_ARCHITECTURES         arm64;armv7;armv7s
        IOS_SIMULATOR_ARCHITECTURES      i386;x86_64
        OPENGL_ES_VERSION                2.0
        VTK_ANDROID_BUILD                OFF
        VTK_Group_Imaging                OFF
        VTK_Group_MPI                    OFF
        VTK_Group_Qt                     OFF
        VTK_Group_Rendering              ON
        VTK_Group_StandAlone             ON
        VTK_Group_Tk                     OFF
        VTK_Group_Views                  OFF
        VTK_Group_Web                    OFF
        VTK_IOS_BUILD                    OFF
        VTK_PYTHON_VERSION               2
        VTK_SMP_IMPLEMENTATION_TYPE      Sequential
        VTK_USE_LARGE_DATA               OFF
        VTK_WRAP_JAVA                    OFF
        VTK_WRAP_PYTHON                  ON
        VTK_WRAP_TCL                     OFF

        BUILD_DOCUMENTATION: Build the VTK documentation
        Press [enter] to edit option                                                 CMake Version 3.2.2
        Press [c] to configure
        Press [h] for help           Press [q] to quit without generating
        Press [t] to toggle advanced mode (Currently Off)

    # Build VTK using make - this will take a while :(
    mlawson@mbp:VTK-6.2.0$ sudo make

Linux
-----
Install VTK and the VTK Python bindings using your system\'s package manager or adapt the `OSX`_ instructions described above.

Windows
-------
No instructions available at this time.
