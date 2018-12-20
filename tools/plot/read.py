#!/usr/bin/env python
from struct import unpack
import numpy as np
import os
import glob
import mmap
from collections import OrderedDict
# backport of concurrent.futures in PyPI package "futures"
from concurrent import futures
import itertools
import functools


class read(object):
  def __init__(self,filename,usecache=False):
    if isinstance(filename,(tuple,list)):
      self.filenames = list(filename)
      if len(self.filenames)==0:
        raise IOError("No file given!")
    else:
      self.filenames = glob.glob(filename)
      if len(self.filenames)==0:
        raise IOError("No file found, invalid pattern '%s'" % filename)
    self.filenames.sort()
    for file in self.filenames:
      if not os.path.isfile(file):
        raise IOError("Not a file/missing: '%s'" % file)
    self.usecache=usecache

    # may be reduced if a file can not be openend
    self.N = len(self.filenames)
    #self.filetemplate = self._getFiletemplate(filename)
    #self.filename = filename
    #if self.isTemplate:
    #  self.setRange()
    #  self.frame = self.range[0]
    #else:
    #  self.frame = -1
    #self.last = None

    self._ReadHeader()
    if self.usecache:
      self.cache_ahead = 16
      self.cache_growth_rate = 4
      self.cache_size = self.cache_ahead + 2
      self._resetStats()
      self.files = cache(self.cache_size,self.cache_growth_rate)
    
    self.frame = 0
    self.select(self.frame)
    
  def _resetStats(self):
    self.stats = {key: 0 for key in self.keys()}

  def __iter__(self):
    self.index = 0
    self.select(self.index)
    return self

  def next(self):
    if(self.index >= self.N):
      raise StopIteration
    self.select(self.index)
    self.index += 1
    return self

  def _load(self,n):
    if not ((0 <= n) and (n < self.N)):
      return True
    if n not in self.files:
      preload = [ k for k, v in self.stats.items() if v > 0 ]
      #self.files[n] = lambda: self._openFile(n,preload)
      self.files[n] = lambda: memopen(self.filenames[n],self.dtype,preload)
      return True
    else:
      return False

  def select(self,n):
    oldframe = self.frame
    if n<0:
      n+=self.N
    self.frame = max(min(n,self.N-1),0)
    delta = n-self.frame
    loaded = 0
    if self.usecache:
      for i in itertools.count(n,delta):
        if(self._load(i)):
          loaded += 1
        if((loaded>=self.cache_growth_rate)\
           or (i>=oldframe+delta*self.cache_ahead)):
          break
      try:
        self.file = self.files[self.frame]
      except KeyError:
        raise KeyError("File with number %i not in allowed range [0,%i]." % (self.frame,self.N))
      self._resetStats()
    else:
      self.file = self._openFile(self.frame)
    #self.file = 
    #self.file = self.files[self.frame]
    #self.file = np.memmap(self._getFilename(),dtype=self.dtype,mode='r')
    #self.file = np.fromfile(self._getFilename(),dtype=self.dtype)

  def getFilename(self):
    return self.filenames[self.frame]

  def getRange(self):
    return [0,self.N-1]

  def __len__(self):
    return self.N

  def __setitem__(self, key, val):
    raise Exception("Writing not allowed")

  def __contains__(self,k):
    return k in self.keys()

  def get(self,k,default=None):
    if self.usecache:
      self.stats[k] += 1
    try:
      #self.last = self.file(name)[0].T
      self.last = self.file[k][0].T
      return self.last
    except(ValueError):
      if default is not None:
        return default
      print "Available keys:"
      print '\n'.join(self.keys())
      raise Exception("The key '%s' could not be found." % k)
    except:
      return self.last

  def __getitem__(self,k):
    return self.get(k)

  def keys(self):
    return sorted(self.dtype.names)

  def getNumber(self):
    return int(self.filenames[self.frame].rsplit('_',1)[1].split('.',1)[0])

  def _ReadHeader(self):
    self.dtypes = { \
      1 : lambda d: "%si%i" % (self.endian,self.intsize),
      2 : lambda d: "%sf%i" % (self.endian,self.realsize),
      3 : lambda d: "S%i" % d,
      4 : lambda d: "%si%i" % (self.endian,self.intsize),
      5 : lambda d: "%s%sf%i" % ((d,),self.endian,self.realsize),
      6 : lambda d: "%s%sf%i" % (d,self.endian,self.realsize),
      7 : lambda d: "%s%sf%i" % (d,self.endian,self.realsize),
      8 : lambda d: "%s%sf%i" % (d,self.endian,self.realsize),
      9 : lambda d: "%s%si%i" % ((d,),self.endian,self.intsize),
      10 : lambda d: "%sf%i" % (self.endian,self.realsize),
      11 : lambda d: "%si%i" % (self.endian,self.intsize), 
      12 : lambda d: "%s%sf%i" % (d,self.endian,self.realsize)}

    self.dtype = dict()

    with open(self.filenames[0],'rb') as m:
      f = mmap.mmap(m.fileno(),0,mmap.MAP_SHARED,mmap.PROT_READ)
      self.magic,self.endian,self.version,self.realsize,self.intsize \
        = f.read(6),f.read(2),f.read(1),f.read(2),f.read(2)
      if(self.endian=='II'): self.endian = '<'
      else: self.endian = '>'
      self.version = unpack('%sB' % self.endian,self.version)[0]
      self.intsize = int(self.intsize)
      self.realsize = int(self.realsize)
      #print self.magic,self.endian,self.version,self.realsize,self.intsize

      keylen = f.read(self.intsize)
      while len(keylen)==4:
        keylen, = unpack('%si' % self.endian,keylen)
        key,t,datalen = unpack("%s%isii" % (self.endian,keylen),f.read(keylen+2*self.intsize))
        if datalen>0:
          dims = datalen
          if t in [6,7,8,12]:
            translate = {6:2, 7:3, 8:4, 12:5}
            d = translate[t]
            datalen -= d*self.intsize
            dims = tuple(unpack('%s%ii' % (self.endian,d),f.read(d*self.intsize)))[::-1]
          offset = f.tell()
          self.dtype[key] = (self.dtypes[t](dims),offset)
          #print "%s: " % key, self.dtype[key]
          f.seek(datalen,os.SEEK_CUR)
        keylen = f.read(self.intsize)

    self.dtype = np.dtype(self.dtype)

  def _openFile(self,i):
    try:
      f = np.memmap(self.filenames[i],dtype=self.dtype,mode='r')
      return f
    except:
      print "Warning: Could not open file '%s'." % self.filenames[i]

  def _openFiles(self):
    self.files = list()
    self.frame = 0
    for filename in self.filenames:
      try:
        file = np.memmap(filename,dtype=self.dtype,mode='r')
      except:
        print "Warning: Could not open file '%s'." % filename
        pass
      self.files.append(file)
      self.frame += 1
    self.N = len(self.files)

class cache(OrderedDict):
  def __init__(self,maxitems,max_workers=None,*args,**kwargs):
    super(cache,self).__init__(*args,**kwargs)
    #self.pool = futures.ProcessPoolExecutor(max_workers=max_workers)
    self.pool = futures.ThreadPoolExecutor(max_workers=max_workers)
    self.maxitems = maxitems

  def __setitem__(self, key, val):
    while len(self) >= self.maxitems:
      super(cache,self).popitem(False)
    fut = self.pool.submit(val)
    super(cache,self).__setitem__(key,fut)

  def __getitem__(self, key):
    #print "CacheSize: ", len(self)
    return super(cache,self).__getitem__(key).result()

  def __del__(self):
    self.pool.shutdown(False)

def memopen(fname,dtype,preload=None):
  f = np.memmap(fname,dtype=dtype,mode='r')
  #@memoize
  #def get(key):
  #  return f[key]
  if preload is not None:
  #  # Ask for request, so the page gets cached.
  #  #for k in preload:
  #  #  get(k)
    f[preload]
  #return get
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
#print d['/timedisc/density']
