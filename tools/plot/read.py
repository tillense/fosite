#!/usr/bin/env python
from struct import unpack, calcsize
import numpy as np
import os
import glob
import mmap
from collections import OrderedDict
# backport of concurrent.futures in PyPI package "futures"
#from concurrent import futures
import itertools
import functools


"""@package read
This package provides a wrapper for the binary and xdmf output.

In order to load some data.
"""


class read(object):

    def __init__(self, filename, usecache=False, mode='r'):
        if isinstance(filename, (tuple, list)):
            self.filenames = list(filename)
            if len(self.filenames) == 0:
                raise IOError("No file given!")
        else:
            self.filenames = glob.glob(filename)
            if len(self.filenames) == 0:
                raise IOError("No file found, invalid pattern '%s'" % filename)
        self.filenames.sort()
        for file in self.filenames:
            if not os.path.isfile(file):
                raise IOError("Not a file/missing: '%s'" % file)
        self.usecache = usecache
        self.mode = mode

        # may be reduced if a file can not be openend
        self.N = len(self.filenames)
        #self.filetemplate = self._getFiletemplate(filename)
        #self.filename = filename
        # if self.isTemplate:
        #  self.setRange()
        #  self.frame = self.range[0]
        # else:
        #  self.frame = -1
        #self.last = None

        self._ReadHeader()
        if self.usecache:
            self.cache_ahead = 16
            self.cache_growth_rate = 4
            self.cache_size = self.cache_ahead + 2
            self._resetStats()
            self.files = cache(self.cache_size, self.cache_growth_rate)

        self.frame = 0
        self.select(self.frame)

    def _resetStats(self):
        self.stats = {key: 0 for key in list(self.keys())}

    def __iter__(self):
        self.index = 0
        self.select(self.index)
        return self

    def __next__(self):
        if(self.index >= self.N):
            raise StopIteration
        self.select(self.index)
        self.index += 1
        return self

    def _load(self, n):
        if not ((0 <= n) and (n < self.N)):
            return True
        if n not in self.files:
            preload = [k for k, v in list(self.stats.items()) if v > 0]
            #self.files[n] = lambda: self._openFile(n,preload)
            self.files[n] = lambda: memopen(
                self.filenames[n], self.dtype, preload)
            return True
        else:
            return False

    def select(self, n):
        oldframe = self.frame
        if n < 0:
            n += self.N
        self.frame = max(min(n, self.N - 1), 0)
        delta = n - self.frame
        loaded = 0
        if self.usecache:
            for i in itertools.count(n, delta):
                if(self._load(i)):
                    loaded += 1
                if((loaded >= self.cache_growth_rate)
                   or (i >= oldframe + delta * self.cache_ahead)):
                    break
            try:
                self.file = self.files[self.frame]
            except KeyError:
                raise KeyError(
                    "File with number %i not in allowed range [0,%i]." %
                    (self.frame, self.N))
            self._resetStats()
        else:
            self.data = self._openFile(self.frame)

    def getFilename(self):
        return self.filenames[self.frame]

    def getRange(self):
        return [0, self.N - 1]

    def __len__(self):
        return self.N

    def __setitem__(self, key, val):
        raise Exception("Write by getting a reference to the object first!")

    def __contains__(self, k):
        return k in list(self.keys())

    def get(self, k, default=None):
        if self.usecache:
            self.stats[k] += 1
        try:
            return self.data[k].T
        except(KeyError):
            if default is not None:
                return default
            print("Available keys:")
            print(('\n'.join(list(self.keys()))))
            raise Exception("The key '%s' could not be found." % k)

    def __getitem__(self, k):
        return self.get(k)

    def keys(self):
        return sorted(list(self.dtype.keys()))

    def getNumber(self):
        return int(
            self.filenames[
                self.frame].rsplit(
                '_',
                1)[1].split(
                '.',
                1)[0])

    def _ReadHeader(self):
        """ Reads out the header of the binary file.

        The structure of the fixed size binary data is written down
        \link here \endlink.
        """
        self.dtypes = {
            1: lambda d: "%si%i" % (self.endian, self.intsize),
            2: lambda d: "%sf%i" % (self.endian, self.realsize),
            3: lambda d: "S%i" % d,
            4: lambda d: "%si%i" % (self.endian, self.intsize),
            5: lambda d: "%s%sf%i" % ((d,), self.endian, self.realsize),
            6: lambda d: "%s%sf%i" % (d, self.endian, self.realsize),
            7: lambda d: "%s%sf%i" % (d, self.endian, self.realsize),
            8: lambda d: "%s%sf%i" % (d, self.endian, self.realsize),
            9: lambda d: "%s%si%i" % ((d,), self.endian, self.intsize),
            10: lambda d: "%sf%i" % (self.endian, self.realsize),
            11: lambda d: "%si%i" % (self.endian, self.intsize),
            12: lambda d: "%s%sf%i" % (d, self.endian, self.realsize)}

        self.dtype = dict()


        with open(self.filenames[0], 'rb') as m:
            f = mmap.mmap(m.fileno(), 0, mmap.MAP_SHARED, mmap.PROT_READ)
            self.magic, self.endian, self.version, self.realsize, self.intsize \
                = f.read(6), f.read(2), f.read(1), f.read(2), f.read(2)
            if(self.endian == b'II'):
                self.endian = '<'
            else:
                self.endian = '>'
            self.version = unpack('%sB' % self.endian, self.version)[0]
            self.intsize = int(self.intsize)
            self.realsize = int(self.realsize)

            keylen = f.read(self.intsize)
            while len(keylen) == 4:
                keylen, = unpack('%si' % self.endian, keylen)
                key, t, datalen = unpack(
                    "%s%isii" %
                    (self.endian, keylen), f.read(
                        keylen + 2 * self.intsize))
                key = key.decode("utf-8")
                if datalen > 0:
                    dims = datalen
                    if t in [6, 7, 8, 12]:
                        translate = {6: 3, 7: 3, 8: 4, 12: 5}
                        d = translate[t]
                        datalen -= d * self.intsize
                        dims = tuple(
                            unpack(
                                '%s%ii' %
                                (self.endian,
                                 d),
                                f.read(
                                    d *
                                    self.intsize)))[
                            ::-
                            1]
                    offset = f.tell()
                    self.dtype[key] = (self.dtypes[t](dims), np.int64(offset))
                    f.seek(datalen, os.SEEK_CUR)
                keylen = f.read(self.intsize)

    def _openFile(self, i):
        self.file = np.memmap(self.filenames[i], mode=self.mode)
        d = dict()
        for key, value in self.dtype.items():
            dtype, offset = value
            dt = np.dtype(dtype)
            dtscalar = np.dtype(dtype.split(')')[-1])
            d[key] = np.ndarray(
                shape=dt.shape,
                buffer=self.file,
                dtype=dtscalar,
                offset=offset)
        return d

class cache(OrderedDict):

    def __init__(self, maxitems, max_workers=None, *args, **kwargs):
        super(cache, self).__init__(*args, **kwargs)
        #self.pool = futures.ProcessPoolExecutor(max_workers=max_workers)
        self.pool = futures.ThreadPoolExecutor(max_workers=max_workers)
        self.maxitems = maxitems

    def __setitem__(self, key, val):
        while len(self) >= self.maxitems:
            super(cache, self).popitem(False)
        fut = self.pool.submit(val)
        super(cache, self).__setitem__(key, fut)

    def __getitem__(self, key):
        # print "CacheSize: ", len(self)
        return super(cache, self).__getitem__(key).result()

    def __del__(self):
        self.pool.shutdown(False)


def memopen(fname, dtype, preload=None):
    f = np.memmap(fname, dtype=dtype, mode='r')
    #@memoize
    # def get(key):
    #  return f[key]
    if preload is not None:
        #  Ask for request, so the page gets cached.
        #  for k in preload: get(k)
        f[preload]
    # return get
    return f


def memoize(obj):
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer

#d = read('dmr_0000.bin')
# print d['/timedisc/density']
