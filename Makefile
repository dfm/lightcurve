
default: interface
	python setup.py build_ext --inplace

interface: lightcurve/periodogram.f90
	python setup.py interface

clean:
	rm -rf lightcurve/*.pyc lightcurve/*.so lightcurve/_periodogram-f2pywrappers.f \
		lightcurve/_periodogrammodule.c lightcurve/periodogram.pyf
