"""
Let's find us some variable stars.

"""

__all__ = ["lc_model", "fit", "chi2", "find_period", "get_model"]

from multiprocessing import Pool

import numpy as np
from numpy.linalg import lstsq
import scipy.optimize as op

import _periodogram

_default_order = 12

_var_m, _var_b = -1.82858601548, 3.23343273764


def lc_model(omega, amplitudes, order):
    a = amplitudes
    return lambda t: a[0] + np.sum([a[2*i+1] * np.sin((i+1)*omega*t)+\
                                    a[2*i+2] * np.cos((i+1)*omega*t)
                                            for i in range(order)],axis=0)


def get_model(period, time, flux, ferr, order=None):
    r = fit(2*np.pi/period, time, flux, ferr, order=order)
    return r


def chi2(model, time, flux, ferr):
    if ferr is None:
        ferr = np.ones_like(flux)
    return np.sum( (flux - model(time))**2 / ferr**2)


class _fit_wrapper(object):
    def __init__(self, time, flux, ferr, order):
        self.time = time
        self.flux = flux
        self.ferr = ferr
        self.order = order

    def __call__(self, omega):
        # Check to make sure that the period isn't _very_ close to a day.
        if np.abs(omega - 2 * np.pi) < 0.05:
            return 1e10

        # Do the fit.
        r = 0.0
        for k in self.time:
            a, chi2, info = _periodogram.get_chi2(omega, self.time[k],
                    self.flux[k], self.ferr[k], self.order)

            if info != 0:
                print "Fit failed"
                return 1e10

            # Calculate the amplitudes and make sure that the 1st order
            # dominates.
            amp = a[1::2] ** 2 + a[2::2] ** 2
            if a[0] < 0 or np.any(amp[0] < amp[1:]):
                return 1e10

            r += chi2

        return r


class _op_wrapper(object):
    def __init__(self, time, flux, ferr, order):
        self.time = time
        self.flux = flux
        self.ferr = ferr
        self.order = order

    def chi(self, w):
        return sum([_periodogram.get_chi(w, self.time[k], self.flux[k],
                self.ferr[k], self.order)[1] for k in self.time])

    def __call__(self, omega):
        res = op.leastsq(self.chi, omega, full_output=True)
        return res[0], np.sum(res[2]["fvec"] ** 2)


def find_period(time, flux, ferr=None, order=None, N=30, Ts=[0.2, 1.3],
        pool=None):
    """
    Find the best fit period of an RR Lyrae lightcurve by doing a grid
    search then a non-linear refinement step.

    """
    if order is None:
        order = _default_order

    # Set up a grid to do a grid search in frequency space.
    domega = 0.2 / min([time[k].max() - time[k].min() for k in time])
    omegas = 2 * np.pi * np.arange(1 / max(Ts), 1 / min(Ts), domega)

    # Do a parallel grid search.
    if pool is None:
        pool = Pool()
    chi2 = map(_fit_wrapper(time, flux, ferr, order), omegas)

    # Sort the results by chi2.
    inds = np.argsort(chi2)

    # Refine the top `N` best fits.
    ref = map(_op_wrapper(time, flux, ferr, order), omegas[inds[:N]])

    # Clean up... otherwise we get the error: _Too many open files_.
    pool.close()
    pool.join()
    del pool

    # Sort the refinements and return the best.
    omega = min(ref, key=lambda x: x[1])
    return 2 * np.pi / float(omega[0])
