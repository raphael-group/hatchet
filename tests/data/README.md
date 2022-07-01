Test Data
=========
This folder contains test data that is useful for running packaged unit tests for Hatchet.

Note that tests and their associated data files are for CI testing, or **NOT** installed when installing the
Hatchet package using setuptools, either in regular or develop mode.

The "fw" (fixed-width) subdirectory contains required files for the scripts "test_solver.py" and "test_steps.py".
The "vw" (variable-width) subdirectory contains files required for the scripts "test_steps_vw.py" and "test_phase.py".
