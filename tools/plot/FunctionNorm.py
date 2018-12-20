from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy import ma
import matplotlib.cbook as cbook
from matplotlib.colors import Normalize

def ArcsinhNorm(scale=1.,*args,**kwargs):
  fn = lambda d: np.arcsinh(d*scale)
  invfn = lambda d: np.sinh(d)/scale
  return FunctionNorm(fn,invfn,*args,**kwargs)

class FunctionNorm(Normalize):
    """
    Normalize a given value to the ``[0, 1]`` interval with a power-law
    scaling. This will clip any negative data points to 0.
    """
    def __init__(self, fn, invfn, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self, vmin, vmax, clip)
        self.fn = fn
        self.invfn = invfn

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax = self.vmin, self.vmax
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            result.fill(0)
        else:
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)
            resdat = result.data
            resdat = self.fn(resdat)
            resdat -= self._lower
            resdat /= (self._upper - self._lower)

            result = np.ma.array(resdat, mask=result.mask, copy=False)
        if is_scalar:
            result = result[0]
        return result

    def _transform_vmin_vmax(self):
        """
        Calculates vmin and vmax in the transformed system.
        """
        vmin, vmax = self.vmin, self.vmax
        arr = np.array([vmax, vmin]).astype(np.float)
        self._upper, self._lower = self.fn(arr)
    
    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")

        vmin, vmax = self.vmin, self.vmax
        if cbook.iterable(value):
            val = ma.asarray(value)
            val = val * (self._upper - self._lower) + self._lower
            return self.invfn(val)
        else:
            value = value * (self._upper - self._lower) + self._lower
            return self.invfn(value)

    def autoscale(self, A):
        """
        Set *vmin*, *vmax* to min, max of *A*.
        """
        self.vmin = ma.min(A)
        self.vmax = ma.max(A)
        self._transform_vmin_vmax()

    def autoscale_None(self, A):
        ' autoscale only None-valued vmin or vmax'
        if self.vmin is None and np.size(A) > 0:
            self.vmin = ma.min(A)

        if self.vmax is None and np.size(A) > 0:
            self.vmax = ma.max(A)
        self._transform_vmin_vmax()

