#!/usr/bin/env python
from read import read
import numpy as np
import sys
from string import join

if __name__=="__main__":
  if len(sys.argv) > 1:
    f = read(sys.argv[1])
    keys = f.keys()

    # filter out directories
    keys = [ k for k in keys if not k.endswith('/.')]

    # get /config directory
    config = [ k for k in keys if k.startswith('/config/') ]
    prefix = ''
    for k in sorted(config):
      val = f[k]
      key = k[7:]
      p, postfix = key.rsplit('/',1)
      p = p[1:]
      if prefix!=p:
        prefix = p
        print "#  [{0}]".format(prefix)
      print "#{:>23}:{:>14}".format(postfix,val)

    # get all data items (the rest)
    data = [ k for k in keys if not k.startswith('/config/') ]
    outkeys = list()
    for k in data:
      val = f[k]
      if type(val) == np.ndarray and len(val.shape)==2 \
        and val.shape[0]==f['/config/mesh/inum'] and val.shape[1]==f['/config/mesh/jnum']:
        outkeys.append(k)
    outkeys = sorted(outkeys)
    ks = [ i.rsplit('/',1)[1] for i in outkeys ]
    print "# ",join(ks,'\t')
    for i in f:
      vals = [f[k] for k in outkeys]
      for a in np.nditer(vals):
        a = [ "%.4e" % np.float64(i) for i in a ]
        print " ",join(a, '\t')
      print "\n"

