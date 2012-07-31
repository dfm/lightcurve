import time
import numpy as np

from lightcurve.lightcurve import find_period


N, w, o = 5, 0.1, 3
t = np.linspace(0, 10, 100)
error = 1.0 * np.ones_like(t)
f = 45. * np.sin(2 * np.pi * t / 0.5 + 3) + error * np.random.randn(len(t))

s = time.time()
print find_period({"g": t}, {"g": f}, ferr={"g": error}, order=o)
print time.time() - s
