#!/usr/bin/env python

import os
import sys

from numpy.distutils.core import setup, Extension
import numpy.distutils.system_info as sysinfo


if "interface" in sys.argv:
    # Generate the Fortran signature/interface.
    cmd = "cd lightcurve;"
    cmd += "f2py periodogram.f90 -m _periodogram -h periodogram.pyf"
    cmd += " --overwrite-signature"
    os.system(cmd)
    sys.exit(0)

# Find LAPACK.
lapack_info = sysinfo.get_info("lapack")
assert len(lapack_info) > 0, "You need LAPACK."

# Define the Fortran extension.
f_ext = Extension("lightcurve._periodogram", ["lightcurve/periodogram.pyf",
    "lightcurve/periodogram.f90"], **lapack_info)

setup(
    name="lightcurve",
    author="Dan Foreman-Mackey",
    packages=["lightcurve"],
    ext_modules=[f_ext],
)
