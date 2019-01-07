#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib as mpl
#from coplot import *
from coroutine import coroutine
from read import read
import argparse
import itertools as it
from scipy.optimize import curve_fit
import itertools

from matplotlib.colors import LogNorm,SymLogNorm
from FunctionNorm import FunctionNorm,ArcsinhNorm

from centercmap import centercmap
stdcmap = plt.cm.viridis

try:
  from tqdm import tqdm
except ImportError:
  tqdm = lambda a: a


from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=False)
#rc('text.latex', preamble=r"\usepackage{wasysym}")
#rc('text.latex', unicode=True)
# overwrite some standard keymapings, so they don't collide with our
# own.
# See: http://matplotlib.github.com/users/navigation_toolbar.html
rc('keymap',back='c',forward='v',zoom='z',home='h')

# Add Keypress to plot coroutine
def connectKey(fig,key,action):
  def press(event):
    if(event.key == key):
      action()
      fig.canvas.draw()
  fig.canvas.mpl_connect('key_press_event', press)

@coroutine
def once(send,data):
  first = True
  while True:
    d = (yield)
    if first:
      send(data)
    send(d)

@coroutine
def sum(send, axis=None):
  while True:
    d = (yield)
    send(np.sum(d,axis))

@coroutine
def gaussian_filter(send,*args,**kwargs):
  from scipy.ndimage.filters import gaussian_filter
  while True:
    d = (yield)
    a = gaussian_filter(d,*args,**kwargs)
    send(a)

@coroutine
def gaussiankernel_filter(send,stddev,*args,**kwargs):
  from astropy.convolution import Gaussian1DKernel, convolve
  #http://docs.astropy.org/en/stable/convolution/
  g = Gaussian1DKernel(stddev=stddev)
  kwargs.setdefault("boundary","extend")
  while True:
    d = (yield)
    a = convolve(d, g, *args,**kwargs)
    send(a)

@coroutine
def median_filter(send,*args,**kwargs):
  from scipy.signal import medfilt
  while True:
    d = (yield)
    a = medfilt(d,*args,**kwargs)
    send(a)

@coroutine
def uniform_filter(send,*args,**kwargs):
  from scipy.ndimage.filters import uniform_filter
  while True:
    d = (yield)
    a = uniform_filter(d,*args,**kwargs)
    send(a)

@coroutine
def mean(send, axis=None):
  while True:
    d = (yield)
    send(np.mean(d,axis))

@coroutine
def customfilter(send,fn):
  while True:
    d = (yield)
    if fn(d):
      send(d)

@coroutine
def min(send, axis=None):
  while True:
    d = (yield)
    send(np.min(d,axis))

@coroutine
def max(send, axis=None):
  while True:
    d = (yield)
    send(np.max(d,axis))

minmaxmean = lambda d: (np.min(d),np.mean(d),np.max(d))

@coroutine
def minmax(send, axis=None):
  while True:
    d = (yield)
    send((np.min(d,axis),np.max(d,axis)))

@coroutine
def clip(send, min, max):
  while True:
    d = (yield)
    send(np.clip(d, min, max))

@coroutine
def symclip(send):
  while True:
    d = (yield)
    l = np.min([np.abs(np.min(d)),np.max(d)])
    send(np.clip(d, -l, l))

@coroutine
def printer(send,fn=lambda d: d):
  while True:
    d = (yield)
    if 'f' in send:
      f = send['f']
    print fn(d)
    send(d)

@coroutine
def produce(send,d):
  while True:
    (yield)
    send(d)

@coroutine
def generate(send,d=None):
  while True:
    r = (yield)
    if d is None:
      d = r
    for i in d:
      send(i)

def getval(p,*args):
  @coroutine
  def exitandraise(send):
    d = (yield)
    raise StopIteration(d)
  try:
    (p | exitandraise()).send(*args)
  except StopIteration as inst:
    if(len(inst.args)==1):
      return inst.args[0]
    else:
      return inst.args


@coroutine
def pipe(send,p=None):
  if p is None:
    while True:
      send((yield))
  else:
    @coroutine
    def helper(send2):
      while True:
        d2 = (yield)
        send(d2)
    subpipe = p | helper()
    while True:
      d = (yield)
      subpipe.send(d)

def inject(send,p):
  @coroutine
  def helper(send2):
    while True:
      d2 = (yield)
      send(d2)
  subpipe = p | helper()
  return subpipe

@coroutine
def call(send,f,*args,**kwargs):
  while True:
    d = (yield)
    f(*args,**kwargs)
    send(d)

@coroutine
def iterate(send):
  while True:
    d = (yield)
    for i in d:
      send(i)

@coroutine
def select(send,frame=None):
  f = send['f']
  while True:
    d = (yield)
    if frame is None:
      if isinstance(d,str):
        i = int(d.rstrip(".bin").split("_")[-1])
      else:
        i = d
    else:
      i = frame
    if isinstance(f,tuple) or isinstance(f,list):
      for file in f:
        file.select(i)
    else:
      f.select(i)
    send(d)

@coroutine
def imshow(send,aspect=1,cmap=stdcmap,cax=plt,clim=None,cscale=None,norm=None,xscale=None,yscale=None,xlabel='x',ylabel='y',clabel=None,*args,**kwargs):
  from matplotlib.colors import LogNorm
  kwargs.setdefault("rasterized",True)
  p = None
  fig,sub = send['fig'],send['sub']
  if clim==None:
    rescale = lambda p: p.set_clim(np.min(p.get_array()),np.max(p.get_array()))
  else:
    rescale = lambda p: p.set_clim(clim[0],clim[1])
  ex = kwargs.get('extent')
  if ex is not None:
    print ex
    f = send['f']
    try:
      kwargs['extent'] = ex(f)
    except:
      pass
    
  while True:
    d = (yield)
    if p==None:
      if cscale=='log':
        norm=LogNorm()
      elif cscale=='symlog':
        norm=SymLogNorm(0.01)
      elif cscale=='arcsinh':
        norm=ArcsinhNorm()
      p = sub.imshow(d.T,
        interpolation='none',
        cmap=cmap,
        clim=clim,
        norm=norm,
        origin='lower',
        *args,
        **kwargs)
      ax = p.axes
      if aspect is not None:
        ax.set_aspect(aspect)
      ax.set_xlabel(xlabel)
      ax.set_ylabel(ylabel)
      if xscale is not None:
        ax.set_xscale(xscale)
      if yscale is not None:
        ax.set_yscale(yscale)
      cbar = cax.colorbar(p)
      if clabel!=None and cbar is not None:
        cbar.set_label(clabel)
      if clim==None:
        rescale2 = lambda: p.set_clim(np.min(p.get_array()),np.max(p.get_array()))
      else:
        rescale2 = lambda: p.set_clim(clim[0],clim[1])
      connectKey(fig,'backspace',rescale2)
      rescale2()
    else:
      p.set_array(d.T)
      rescale(p)
    send(d)

@coroutine
def surface(send,f,fig,sub,aspect=1,cmap=stdcmap,cax=plt,clim=None,cscale=None,norm=None,xscale=None,yscale=None,*args,**kwargs):
  from matplotlib.colors import LogNorm
  p = None
  while True:
    d = (yield)
    if p==None:
      if cscale=='log':
        norm=LogNorm()
      bccart = f['/mesh/bary_centers']
      X = bccart[:,:,0]
      Y = bccart[:,:,1]
      #p = sub.plot_surface(
      p = sub.plot_trisurf(
        X.flatten(), Y.flatten(),
        d.flatten(),
        linewidth=0.,
#        rstride = 1,
#        cstride = 1,
        norm=norm,
        cmap=cmap,
        *args,
        **kwargs)
      ax = p.axes
      ax.set_aspect(aspect)
      if xscale is not None:
        ax.set_xscale(xscale)
      if yscale is not None:
        ax.set_yscale(yscale)
      cax.colorbar(p)
      if clim==None:
        rescale = lambda: p.set_clim(np.min(p.get_array()),np.max(p.get_array()))
      else:
        rescale = lambda: p.set_clim(clim[0],clim[1])
      connectKey(fig,'backspace',rescale)
      rescale()
    else:
      #p.set_array(d.T)
      p.set_array(d.flatten())
    send(d)

@coroutine
def plot(send,xdata=None,xlabel='x',ylabel='y',fmt='',xscale=None,yscale=None,xlim=None,ylim=None,aspect=None,text = lambda f:'n = %i\nt = %.2e' % (f.frame, f['/timedisc/time']),**kwargs):
  f = send.get('f')
  fig,sub = send['fig'],send['sub']
  p, = sub.plot([],[],fmt,**kwargs)
  ax = p.axes
  ax.grid(True)
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  if aspect is not None:
    ax.set_aspect(aspect)
  title = ax.text(0.95,0.9,'',transform=ax.transAxes,ha='right')
  def rescale():
    minmax = lambda data: (np.min(data), np.max(data))
    x, y = p.get_data()
    if(xlim==None):
      sub.set_xlim(*minmax(x))
    else:
      sub.set_xlim(xlim)
    if(ylim==None):
      sub.set_ylim(*minmax(y))
    else:
      sub.set_ylim(ylim)
  connectKey(fig,"backspace",rescale)
  first = True

  while True:
    d = (yield)
    try:
      p.set_data(*d)
    except:
      if xdata is None:
        p.set_data(np.arange(len(d)),d)
      else:
        p.set_data(xdata,d)
    if f is not None:
      title.set_text(text(f))
    if first:
      if xscale is not None:
        sub.set_xscale(xscale)
      if yscale is not None:
        sub.set_yscale(yscale)
      rescale()
      first = not first
    send(d)

@coroutine
def multiplot(send,xlabel='x',ylabel='y',xlim=(None,None),ylim=(None,None),fmts=(),labels=(),repeat=0,xscale=None,yscale=None,loc='best',ncol=1,ltitle=None,rtitle=None,color='rgbcmykw',marker=None,linestyle='-',**kwargs):
  #marker='xo+ps*'
  fig,sub = send['fig'],send['sub']
  minmax = lambda d: (np.min(d),np.max(d))
  ax = None
  p = []
  def rescale():
    xl = minmax(np.concatenate([l.get_xdata() for l in p]))
    yl = minmax(np.concatenate([l.get_ydata() for l in p]))
    ax.set_xlim(xl)
    ax.set_ylim(yl)
#    sub.set_xlim(auto=True)
#    sub.set_ylim(auto=True)
  def press(event):
    print 'press', event.key
    if(event.key in keys):
      keys[event.key]()
      fig.canvas.draw()
  keys = dict(backspace = lambda: rescale())
  fig.canvas.mpl_connect('key_press_event', press)
  for i,fmt,label in it.izip_longest(it.count(0),fmts,labels,fillvalue=None):
    x,y = (yield)
    if repeat==0 or i<repeat:
      c = color[len(p) % len(color)]
      ls = linestyle[len(p) % len(linestyle)]
      if marker is not None:
        m = marker[len(p) % len(marker)]
        #print marker,len(p),len(marker),len(p) % len(marker),m
        #kwargs.setdefault('marker',m)
      else:
        m= None
      kwargs.setdefault('markevery',int(len(y)/10.))
      if fmt != None:
        p.append(sub.plot(x,y,fmt,label=label,color=c,marker=m,linestyle=ls,**kwargs)[0])
      else:
        p.append(sub.plot(x,y,label=label,color=c,marker=m,linestyle=ls,**kwargs)[0])
      if len(p)==1:
        ax = p[0].axes
        ax.grid(True)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xscale is not None:
          ax.set_xscale(xscale)
        if yscale is not None:
          ax.set_yscale(yscale)
      ax.legend(loc=loc,ncol=ncol)
      rescale()
      if ltitle is not None:
        ax.set_title(ltitle(),loc='left')
      if rtitle is not None:
        ax.set_title(rtitle(),loc='right')
      ax.set_xlim(*xlim)
      ax.set_ylim(*ylim)
    else:
      p[i%repeat].set_ydata(y)
    send((x,y))

@coroutine
def makelist(send,n=lambda f: len(f),l=None):
  f = send.get('f')
  l = None
  i = 0
  nolist=False
  while True:
    d = (yield)
    if not isinstance(d, (list,tuple)):
      d = (d,)
      nolist=True
    i += 1
    if l==None:
      l = tuple([e] for e in d)
    else:
      for a,e in zip(l,d):
        a.append(e)
    if(i%n(f)==0):
      if nolist:
        l, = l
      send(l)
      l=None

@coroutine
def cut(send,s):
  sel = tuple()
  for i in s:
    if i==None:
      sel += (slice(None,None),)
    else:
      sel += (i,)
  while True:
    d = (yield)
    send(d[sel])

#with interpolation
@coroutine
def streamplot2(send,res=50,clabel=None,xlim=None,ylim=None,scale=1.,cmap=stdcmap,color=lambda u,v: np.sqrt(u**2+v**2),lw=lambda c: 4*c/np.nanmax(c),*args,**kwargs):
  import streamplot as sp
  from scipy.interpolate import griddata
  f,fig,sub = send['f'],send['fig'],send['sub']
  minmax = lambda d: [np.min(d),np.max(d)]
  bc = f['/mesh/bary_centers']/scale
  bx,by = bc[:,:,0].flatten(),bc[:,:,1].flatten()
  if xlim is None:
    xlim = minmax(bx)
  if ylim is None:
    ylim = minmax(by)
  x,y = np.linspace(xlim[0],xlim[1],res),np.linspace(ylim[0],ylim[1],res)
  xi,yi = np.meshgrid(x,y)

  #mask for inner hole
  minr2 = np.min(bx**2+by**2)
  mask = xi**2+yi**2 < minr2

  minmax = lambda d: [np.nanmin(d),np.nanmax(d)]

  interpolate = lambda d: griddata(zip(bx,by),d.flatten(),(xi,yi))
  firsttime = True
  p = None
  while True:
    u,v = (yield)
    speed = color(u,v,f)
    gu,gv=interpolate(u),interpolate(v)
    gu[mask] = np.nan
    gspeed=interpolate(speed)
    if p!=None:
      p.lines.remove()
      ##p.arrows.remove() # does not work: not implemented for this collection!
      #workaŕound from http://stackoverflow.com/questions/19621521/deleting-streamplots-matplotlib-without-clearing-the-graph
      keep = lambda x: not isinstance(x, mpl.patches.FancyArrowPatch)
      ax = sub.axes
      ax.patches = [patch for patch in ax.patches if keep(patch)]
      #sub.cla() # only available method at the moment
      #firsttime=True
    #p = sub.streamplot(x,y,gu,gv,color=gspeed,linewidth=lw,cmap=cmap,*args,**kwargs)
    #try:
    #  p = sp.streamplot(sub,x,y,gu,gv,color=gspeed,linewidth=lw(gspeed),cmap=cmap,*args,**kwargs)
    #except sp.InvalidIndexError:
    if True:
      print "Could not start at index point. Disabling start points for now"
      kw = kwargs
      kw['start_points'] = None
      kw['density'] = 0.5
      kw['maxlength'] = 4.
      kw['integration_direction'] = 'both'
      p = sp.streamplot(sub,x,y,gu,gv,color=gspeed,linewidth=lw(gspeed),cmap=cmap,*args,**kw)
      
    ax = sub.axes
    if p is not None and firsttime:
      #from mpl_toolkits.axes_grid1 import make_axes_locatable
      ## create an axes on the right side of ax. The width of cax will be 5%
      ## of ax and the padding between cax and ax will be fixed at 0.05 inch.
      #divider = make_axes_locatable(ax)
      #cax = divider.append_axes("right", size="5%", pad=0.05)
      #cbar = plt.colorbar(p.lines,cax=cax)
      cbar = plt.colorbar(p.lines,use_gridspec=True)
      if clabel is not None:
        cbar.set_label(clabel)
      firsttime=False
      ax.set_aspect('equal')
      ax.set_xlim(xlim)
      ax.set_ylim(ylim)
    #else:
      #cbar = plt.colorbar(p.lines)
      #cbar.update_normal(p.lines)
      #cbar.remove()
      #cbar = plt.colorbar(p.lines)
    send((u,v))


@coroutine
def streamplot(send,scale=1.,*args,**kwargs):
  f,fig,sub = send['f'],send['fig'],send['sub']
  bc = f['/mesh/bary_centers']/scale
  x,y = bc[:,0,0],bc[0,:,1]
  p = None
  while True:
    u,v = (yield)
    if p!=None:
#      p.lines.remove()
#      p.arrows.remove() # does not work: not implemented for this collection!
      sub.cla() # only available method at the moment
    p = sub.streamplot(x,y,u.T,v.T,*args,**kwargs)
    send((u,v))

@coroutine
def contour(send,*args,**kwargs):
  from matplotlib.colors import LogNorm
  f,sub = send['f'],send['sub']
  bary = f['/mesh/bary_centers']
  x = bary[:,:,0].T
  y = bary[:,:,1].T
  p = None
  cscale=kwargs.get('cscale')
  if cscale=='log':
    norm=LogNorm()
  else:
    norm=None
  while True:
    d = (yield)
    if p==None:
      p = sub.contour(x,y,d.T,*args,norm=norm,**kwargs)
      #plt.clabel(p,inline=1)
    else:
      p.set_array(d.T)
    send(d)

@coroutine
def nextsub(send):
  subs = send['subs']
  send['sub'] = subs.next()
  while True:
    d = (yield)
    send(d)

@coroutine
def selectsub(send,i):
  grid = send['grid']
  send['sub'] = grid[i]
  while True:
    d = (yield)
    send(d)
  

@coroutine
def pcolormesh(send,cmap=stdcmap,clim=None,xlabel='x',ylabel='y',aspect='equal',scale=1.,clabel=None,xlim=None,ylim=None,text=None,X=None,Y=None,xscale=None,yscale=None,cscale=None,xticks=None,yticks=None,norm=None,ltitle=None,rtitle=None,cax=None,autoscale=False,edgecolors='None',linewidth=0,tcolor='k',zoom=None,nbins=None,*args,**kwargs):
  kwargs.setdefault("rasterized",True)
  p = None
  d = None
  fig,sub = send['fig'],send['sub']
  minmax = lambda d: (np.min(d),np.max(d))
  rescale = lambda: p.set_clim(*minmax(d))
  def press(event):
    #print 'press', event.key
    if(event.key in keys):
      keys[event.key]()
      fig.canvas.draw()
  keys = dict(backspace = lambda: rescale())
  fig.canvas.mpl_connect('key_press_event', press)
  f=send.get('files')
  if f is not None:
    f = f[0]
  while True:
    d = (yield)
    if p==None:
      if cscale=='log':
        norm=LogNorm()
      elif cscale=='symlog':
        norm=SymLogNorm(0.01)
      elif cscale=='arcsinh':
        norm=ArcsinhNorm()
      if X is None: X=send.get('X')
      if Y is None: Y=send.get('Y')
      if X is None:
        X=f['/mesh/grid_x']/scale
      else:
        X=X(f)
      if Y is None: 
        Y=f['/mesh/grid_y']/scale
      else:
        Y=Y(f)
      # Squeeze any additional dimensions with extent=1 from the grid,
      # e.g. in case of 3D->2D data:
      X = np.squeeze(X)
      Y = np.squeeze(Y)
      p = sub.pcolormesh(X,Y,
             d,cmap=cmap,
             edgecolors=edgecolors,
             linewidth=linewidth,
             norm=norm,
             *args,**kwargs)#,
             #picker=True)
      plist = send.get('p')
      if plist is None:
        plist = [p]
      else:
        plist.append(p)
      send['p'] = plist
      ax = p.axes
      ax.set_aspect(aspect)
      ax.set_xlabel(xlabel)
      ax.set_ylabel(ylabel)
      plt.autoscale(tight=True)
      if xscale is not None:
        ax.set_xscale(xscale)
      if yscale is not None:
        ax.set_yscale(yscale)
      if xlim!=None:
        ax.set_xlim(*xlim)
      if ylim!=None:
        ax.set_ylim(*ylim)
      if yticks is not None:
        ax.set_yticks(yticks[0])
        ax.set_yticklabels(yticks[1])
      if zoom is not None:
        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        from mpl_toolkits.axes_grid1.inset_locator import mark_inset
        axins = zoomed_inset_axes(sub, 4, loc=2) # zoom = 6
        pz = axins.pcolormesh(X,Y,
             d,cmap=cmap,
             edgecolors=edgecolors,
             linewidth=linewidth,
             norm=norm,
             *args,**kwargs)
        pz.set_clim(*clim)
        axins.set_xlim(*zoom[0])
        axins.set_ylim(*zoom[1])
        plt.xticks(visible=False)
        plt.yticks(visible=False)
        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(sub, axins, loc1=3, loc2=1, fc="none", ec="0.5")

        #fix cax
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(sub)
        cax = divider.append_axes("right", size="5%", pad=0.05)

      if nbins is not None:
        #sub.locator_params(nbins=nbins)
        #from matplotlib.ticker import MaxNLocator
        ax.locator_params(nbins=nbins)
        #ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='both'))
        #ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='both'))
      textbox = ax.text(0.95,0.90,'',color=tcolor,transform=ax.transAxes,ha='right')
      if isinstance(cax,(int,long)):
        grid = send['grid']
        cax = grid.cbar_axes[cax]
      if cax is None:
        #from mpl_toolkits.axes_grid1 import make_axes_locatable
        ## create an axes on the right side of ax. The width of cax will be 5%
        ## of ax and the padding between cax and ax will be fixed at 0.05 inch.
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        #cbar = fig.colorbar(p,cax=cax)
        cbar = fig.colorbar(p,use_gridspec=True)
      elif cax!='' and cax!=None:
        cbar = plt.colorbar(p,cax=cax)
        #try:
        #  #cbar = cax.colorbar(p,rasterized=True)
        #  cbar = plt.colorbar(p,cax=cax,rasterized=True)
        #except:
        #  cbar = None
        #  pass
      else:
        cbar = None
      if cbar is not None:
        cbar.update_ticks()
        # to use scaling factor
        if norm is None:
          cbar.formatter.set_scientific(True)
          cbar.formatter.set_powerlimits((0,0))
        cbar.update_ticks()
      if clabel!=None and cbar is not None:
        cbar.set_label(clabel)
      if clim is not None:
        p.set_clim(*clim)
        p.set_clim(*clim)
    else:
      p.set_array(d.ravel())
      if autoscale:
        if clim is None:
          # Yes, two times is correct here.
          # There seems to be a bug? feature? in mpl,
          # and therefore the colorbar is only updated
          # after the second call
          # Its the same, when calling "rescale" per
          # hand, but it doesn't matter there too much.
          p.set_clim(*minmax(d))
          p.set_clim(*minmax(d))
        fig.canvas.draw()
    if ltitle is not None:
      ax.set_title(ltitle(),loc='left')
    if rtitle is not None:
      ax.set_title(rtitle(),loc='right')
    if text is not None:
      textbox.set_text(text(f))
    else:
      time = send['files'][0]['/timedisc/time']
      textbox.set_text('n = %i\nt = %.2e' % (f.getNumber(), time))
    send(d)

@coroutine
def text(send,x,y,label,*args,**kwargs):
  first = True
  while True:
    d = (yield)
    if first:
      sub,f = send['sub'],send.get('f')
      kwargs.setdefault('transform',sub.transAxes)
      kwargs.setdefault('ha','right')
      sub.text(x,y,label(f),*args,**kwargs)
      first = False
    send(d)

@coroutine
def annotate(send,label,xy,xytext,*args,**kwargs):
  first = True
  while True:
    d = (yield)
    if first:
      sub,f = send['sub'],send.get('f')
      kwargs.setdefault('xycoords','axes fraction')
      sub.annotate(label(f),xy=xy,xytext=xytext,*args,**kwargs)
      first = False
    send(d)

@coroutine
def colorbar(send,i=0,label="",*args,**kwargs):
  fig = send.get('fig')
  p = send.get('p')
  grid = send['grid']
  try:
    cax = grid.cbar_axes[i]
  except:
    import matplotlib as mpl
    cax,kw = mpl.colorbar.make_axes([ax for ax in grid.flat])
  #cbar = cax.colorbar(p[i],*args,**kwargs)
  cbar = fig.colorbar(p[i],cax=cax)
  try:
    cax.toggle_label(True)
    cax.axis[cax.orientation].set_label(label)
  except:
    cbar.set_label(label)
    
  #print cbar
  #cbar.update_ticks()
  # to use scaling factor
  #if norm is None:
  #cbar.formatter.set_scientific(True)
  #cbar.formatter.set_powerlimits((0,0))
  while True:
    d = (yield)
    send(d)

@coroutine
def addtitle(send,i=0,label="",*args,**kwargs):
  grid = send['grid']
  grid[i].set_title(label,*args,**kwargs)
  while True:
    d = (yield)
    send(d)

def pcolorflat(*args,**kwargs):
  scale=1.
  if 'scale' in kwargs:
    scale = kwargs['scale']
  #get = lambda f,i: f['/mesh/bary_centers'][:,:,i]
  def get(f,i): 
    if i==0:
      return f['/mesh/grid_x']
    elif i==1:
      return f['/mesh/grid_y']
    else:
      raise "not allowed"
  X = lambda f: np.sqrt(get(f,0)**2+get(f,1)**2)/scale
  def Y(f):
    d = np.arctan2(get(f,1),get(f,0))
    if d[0,-1]<0.:
      d[:,-1] += 2.*np.pi
    d /= np.pi
    return d
  kwargs.update({'X':X,'Y':Y})
  return pcolormesh(*args,**kwargs)

@coroutine
def timelistgrid(send,scale=1.,tscale=1.):
  first=True
  def get(f,i): 
    if i==0:
      return f['/mesh/grid_x']
    elif i==1:
      return f['/mesh/grid_y']
    else:
      raise "not allowed"
  while True:
    time,radius,data = (yield)
    if first:
      #t = []
      #for f in tqdm(file):
      #  t.append(file['/timedisc/time'])
      t = np.zeros(len(time)+1)
      t[0:-1] = time
      t[-1] = (time[-1] + (time[-1]-time[-2]))
      send['X'] = lambda d: np.array(t) / tscale
      send['Y'] = lambda d: radius/scale
      #send['Y'] = lambda d: np.sqrt(get(f,0)**2+get(f,1)**2)[:,0] / scale
      first = not first
    send(data)

@coroutine
def moviesave(send,filename,fps=60,ffmpeg_params=['-preset','slow','-crf','4'],*args,**kwargs):
  from moviepy.editor import VideoClip
  from moviepy.video.io.bindings import mplfig_to_npimage
  fig = send['fig']
  f = send['f']
  #duration = int(len(f)/fps)
  duration = len(f)/fps
  while True:
    d = (yield)
    def make_frame(t):
      i = int(t*fps)
      f.select(i)
      send(d)
      # hack! probably not always right..
      #fig.canvas.draw()
      return mplfig_to_npimage(fig)
    anim = VideoClip(make_frame, duration=duration)
    #ani = animation.FuncAnimation(fig, anifn, frames=frames)#, blit=True, repeat=False)
    #ani.save(filename,*args,**kwargs)
    anim.write_videofile(filename,fps,ffmpeg_params=ffmpeg_params,*args,**kwargs)

@coroutine
def anisave(send,filename,frames,*args,**kwargs):
  fig = send['fig']
  f = send['f']
  while True:
    d = (yield)
    def anifn(i):
      f.select(i)
      send(d)
      # hack! probably not always right..
      #fig.canvas.draw()
      return fig.axes[0],
    ani = animation.FuncAnimation(fig, anifn, frames=frames)#, blit=True, repeat=False)
    ani.save(filename,*args,**kwargs)

@coroutine
def animate(send,*args,**kwargs):
  while True:
    names = (yield)
    def anifn(a):
      send(names)
      # hack! probably not always right..
      fig.canvas.draw()
      return fig.axes[0],
    ani = animation.FuncAnimation(fig, anifn, f, interval=interval, blit=True, repeat=False)


@coroutine
def timeseries(send):
  f = send['f']
#  it = range(len(f))
  while True:
    names = (yield)
    for i in tqdm(f): #it:
      #f.select(i)
      send(names)

@coroutine
def series(send,it):
  f = send['f']
  while True:
    names = (yield)
    try:
      for i in tqdm(it(f)):
        f.select(i)
        send(names)
    except:
      for i in tqdm(it):
        f.select(i)
        send(names)

def getData(f, names):
  if(isinstance(names,list)):
    return list([f[name] for name in names])
  elif(isinstance(names,tuple)):
    return tuple((f[name] for name in names))
  else:
    return f[names]


@coroutine
def get(send):
  f = send['f']
  while True:
    names = (yield)
    send(getData(f,names))

@coroutine
def savetxt(send,name,*args,**kwargs):
  while True:
    d = (yield)
    np.savetxt(name,np.transpose(d),*args,**kwargs)
    send(d)

@coroutine
def loadtxt(send,*args,**kwargs):
  kwargs.setdefault('unpack',True)
  while True:
    d = (yield)
    try:
      d = np.loadtxt(d,*args,**kwargs)
    except:
      d = np.loadtxt(d[0],*args,**kwargs)
    send(d)

@coroutine
def savefig(send,name,*args,**kwargs):
  fig = send['fig']
  kwargs.setdefault("bbox_inches","tight")
  while True:
    d = (yield)
    if isinstance(name,list):
      n = name.pop(0)
    elif(hasattr(name, '__call__')):
      n = name(send['f'])
    else:
      n = name
    #print "Save: %s" % n
    fig.savefig(n,*args,**kwargs)
    send(d)

def savepng(name,*args,**kwargs):
  name = name.rstrip(".png")
  return savefig(name + ".png",*args,**kwargs)

def savepdf(name,*args,**kwargs):
  name = name.rstrip(".pdf")
  return savefig(name + ".png",*args,**kwargs) \
         | savefig(name + ".pdf",*args,**kwargs) 

@coroutine
def axvline(send,x,*args,**kwargs):
  sub = send['sub']
  sub.axvline(x,*args,**kwargs)
  while True:
    d = (yield)
    send(d)

@coroutine
def saveimg(send,name,*args,**kwargs):
  fig = send['fig']
  kwargs.setdefault("bbox_inches","tight")
  while True:
    d = (yield)
    if isinstance(name,list):
      n = name.pop(0)
    elif(hasattr(name, '__call__')):
      n = name(send['f'])
    else:
      n = name
    sp = n.rsplit(".",1)
    if len(sp) != 2:
        raise Exception("No suffix/filetype has been given in '%s'." % n)
    ext = sp[1]
    n = sp[0]
    if ext=='pgf':
      try:
        fig.savefig(n + ".png",*args,**kwargs)
      except:
        pass # dont fail on latex errors?!
    fig.savefig(n + "." + ext,*args,**kwargs)
    send(d)

@coroutine
def savepgf(send,name,*args,**kwargs):
  fig = send['fig']
  kwargs.setdefault("bbox_inches","tight")
  while True:
    d = (yield)
    if isinstance(name,list):
      n = name.pop(0)
    elif(hasattr(name, '__call__')):
      n = name(send['f'])
    else:
      n = name
    #print "Save: %s" % n
    n = n.rstrip(".pgf")
    try:
      fig.savefig(n + ".png",*args,**kwargs)
    except:
      pass # dont fail on latex errors?!
    fig.savefig(n + ".pgf",*args,**kwargs)
    send(d)

@coroutine
def save(send,name,*args,**kwargs):
  while True:
    d = (yield)
    np.save(name,d,*args,**kwargs)
    send(d)

@coroutine
def load(send,*args,**kwargs):
  while True:
    name = (yield)
    d = np.load(name,*args,**kwargs)
    send(d)

@coroutine
def psave(send,name,*args,**kwargs):
  import cPickle as pickle
  while True:
    d = (yield)
    pickle.dump(d,open(name,"wb"),*args,**kwargs)
    send(d)

@coroutine
def pload(send,*args,**kwargs):
  import cPickle as pickle
  while True:
    name = (yield)
    if isinstance(name, (tuple,list)):
      name = name[0]
    d = pickle.load(open(name,"rb"),*args,**kwargs)
    send(d)

@coroutine
def rfft(send):
  from numpy.fft import rfft,rfft2
  while True:
    d = (yield)
    rank = len(d.shape)
    if(rank==1):
      send(rfft(d))
    elif(rank==2):
      send(rfft2(d))
    else:
      throw('Only 1D or 2D ffts available')

@coroutine
def irfft(send):
  from numpy.fft import irfft,irfft2
  while True:
    d = (yield)
    rank = len(d.shape)
    if(rank==1):
      send(irfft(d))
    elif(rank==2):
      send(irfft2(d))
    else:
      throw('Only 1D or 2D ffts available')

@coroutine
def abs(send):
  while True:
    send(np.abs((yield)))

@coroutine
def log10(send):
  while True:
    send(np.log10((yield)))

@coroutine
def log(send):
  while True:
    send(np.log((yield)))

@coroutine
def only(send,i,p):
  d = None
  @coroutine
  def helper(send2):
    while True:
      di = (yield)
      t = tuple()
      for j,e in enumerate(d):
        if j==i:
          t += (di,)
        else:
          t += (e,)
      send(t)
  subpipe = p | helper()
  subpipe.update(send)
  while True:
    d = (yield)
    subpipe.send(d[i])

@coroutine
def waiter(send, n):
  from itertools import count
  for i in count(0):
    d = (yield)
    if  i % n == 0:
      send(d)

@coroutine
def onlyif(send, fn):
  f = send['f']
  while True:
    d = (yield)
    if fn(f):
      send(d)

@coroutine
def each(send, p):
  serial = reflist([])
  @coroutine
  def join(sendj,serial,l):
    while True:
      res = (yield)
      serial.append(res)
      if(len(serial)==l):
        send(tuple(serial.get()))
        serial.set([])
  pipes = []
  while True:
    d = (yield)
    if len(pipes)<len(d):
      for p in pipes:
        p.close()
      pipes = []
      for e in d:
        pipes.append(p() | join(serial,len(d)))
    for e,subpipe in zip(d,pipes):
      send.ref(subpipe)
      subpipe.send(e)
    

# This defines a list reference as new datatype
class reflist:
  def __init__(self,l): self.l = l
  def get(self): return self.l
  def set(self,l): self.l = l
  def append(self,l): self.l.append(l)
  def __len__(self): return len(self.l)

@coroutine
def split(send,*args):
  serial = reflist([])
  @coroutine
  def join(sendj,serial):
    while True:
      res = (yield)
      serial.append(res)
      if(len(serial)==len(args)):
        send(tuple(serial.get()))
        serial.set([])

  while True:
    d = (yield)
    for e,p in zip(d,args):
      subpipe = (p | join(serial))
      send.ref(subpipe)
      subpipe.send(e)

@coroutine
def repeat(send,no):
  while True:
    d = (yield)
    d = (d,)*no
    send(d)

@coroutine
def broadcast(send,*args):
  while True:
    d = (yield)
    for p in args:
      p.update(send)
      p.send(d)
    send(d)

@coroutine
def normalize(send,i=None,norm=lambda d,d0: d/d0):
  if i==None:
    i=0
  while True:
    d = (yield)
    if isinstance(i, int):
      res = [norm(e,d[i]) for e in d]
    else:
      res = [norm(e,i) for e in d]
    send(res)

@coroutine
def normalize2(send, norm=lambda x,y: x/y):
  d0 = None
  while True:
    d = (yield)
    if d0 is None:
      d0 = d
    send(norm(d,d0))

@coroutine
def numexpr(send,expr):
  import numexpr as ne
  while True:
    d = (yield)
#    send(ne.evaluate(expr,local_dict=send))
    if 'f' in send:
      f = send['f']
      vars = dict(d=d)
      ex = ""
      for i,l in enumerate(expr.split("f['")):
        if i==0:
          ex = l
        else:
          l=l.split("']",1)
          label = l[0].rsplit('/',1)[-1]
          ex += label + l[1]
          vars[label] = f[l[0]]
      send(ne.evaluate(ex,vars))
    else:
      send(ne.evaluate(expr))

@coroutine
def apply(send,fn):
  while True:
    d = (yield)
    #try:
    #  res = fn(*d)
    #except:
    #  res = fn(d)
    res = fn(d)
    send(res)

@coroutine
def movingaverage(send,order=3):
  from collections import deque
  last = deque(maxlen=order)
  while True:
    d = (yield)
    last.append(d)
    res = np.sum(last)/np.float64(len(last))
    #print d, "-->", res
    send(res)
    

@coroutine
def forall(send,fn):
  while True:
    d = (yield)
    d = tuple(fn(dd) for dd in d)
    send(d)

@coroutine
def diff_backward(send,order=1):
  from collections import deque
  from scipy.special import binom
  last = deque(maxlen=order+1)
  t0 = None
  while True:
    time, d = (yield)
    last.appendleft(d)
    n = len(last)
    if n > 1:
      ord = n - 1
      dt = time - t0
      res = np.sum([(-1.)**i*binom(ord,i)*last[i] for i in range(n)]) / dt**ord
      send((time,res))
    t0 = time


@coroutine
def timediff(send,reset=None):
  from itertools import count
  if reset is None:
    t0, d0 = None, None
  for i in count(0):
    if reset is not None:
      if i % reset == 0:
        t0, d0 = None, None
    time, d = (yield)
    if not t0==None:
      send((0.5*(time+t0),(d-d0)/(time-t0)))
    t0, d0 = time, d

@coroutine
def timediff2(send):
  while True:
    t, d = (yield)
    tdiff = 0.5*(t[1:]+t[:-1])
    ddiff = (d[1:]-d[:-1])/(t[1:]-t[:-1])
    send((tdiff,ddiff))

@coroutine
def vecabs(send):
  while True:
    d = (yield)
    if isinstance(d,tuple):
      x,y = d
      send(np.sqrt(x**2+y**2))
    else:
      send(np.sqrt(d[:,:,0]**2+d[:,:,1]**2))
      

def diff(f, d, dir):
  dlx = f['/mesh/dlx']
  dly = f['/mesh/dly']
  bhx = f['/mesh/bhx']
  bhy = f['/mesh/bhy']
  dlx = dlx/bhx
  dly = dly/bhy
  if dir==0:
    gx = np.empty(d.shape)
    gx[1:-1,:] = 0.5*(d[2:,:]-d[:-2,:])/dlx[1:-1,:]
    gx[0,:] = (d[1,:] - d[0,:])/dlx[0,:]
    gx[-1,:] = (d[-1,:] - d[-2,:])/dlx[-1,:]
    return gx
  elif dir==1:
    gy = np.empty(d.shape)
    gy[:,1:-1] = 0.5*(d[:,2:]-d[:,:-2])/dly[:,1:-1]
    gy[:,0] = (d[:,1] - d[:,0])/dly[:,0]
    gy[:,-1] = (d[:,-1] - d[:,-2])/dly[:,-1]
    return gy

@coroutine
def grad(send):
  f = send['f']
  while True:
    d = (yield)
    send((diff(f, d,0)/f['/mesh/bhx'],diff(f, d,1)/f['/mesh/bhy']))

@coroutine
def delR(send):
  f = send['f']
  while True:
    d = (yield)
    send(diff(f, d,0)/f['/mesh/bhx'])

def absgrad(f):
  return grad(f) | vecabs()

@coroutine
def curl(send):
  f = send['files'][0]
  bhx = f['/mesh/bhx']
  bhy = f['/mesh/bhy']
  invsqrtg = 1./(bhx*bhy)
  while True:
    vx,vy = (yield)
    # the curl of a plane 2D vector field is a scalar
    send(invsqrtg*(diff(f,bhy*vy,0)-diff(f,bhx*vx,1)))

@coroutine
def div(send):
  f = send['f']
  bhx = f['/mesh/bhx']
  bhy = f['/mesh/bhy']
  invsqrtg = 1./(bhx*bhy)
  while True:
    vx,vy = (yield)
    send(invsqrtg*(diff(f,bhy*vx,0)+diff(f,bhx*bhy,1)))

@coroutine
def cross(send,f,a):
  a1, a2 = a
  while True:
    b1,b2 = (yield)
    send(a1*b2-a2*b1)

@coroutine
def Convert2Cartesian(send):
  f = send['f']
  rot = f['/mesh/rotation']
  while True:
    vxi, veta = (yield)
    vx = np.cos(rot)*vxi + np.sin(rot)*veta
    vy = -np.sin(rot)*vxi + np.cos(rot)*veta
    send((vx,vy))

@coroutine
def inner(send):
  #component sum
  def compsum(l):
    res = l.pop()
    while len(l) > 0:
      res += l.pop()
    return res

  while True:
    d = (yield)
    send(compsum([np.multiply(i,j) for i,j in d]))

@coroutine
def readdirs(send):
  from read import read
  from os.path import join
  while True:
    fnames = (yield)
    if not isinstance(fnames,(tuple,list)): fnames = [fnames]
    fnames = [ join(name,"*.bin") for name in fnames]
    f = [read(fname) for fname in fnames]
    send['files'] = f
    if len(f)==1: f=f[0]
    send['f'] = f
    send(f)


@coroutine
def readbin(send):
  from read import read
  while True:
    fnames = (yield)
    if not isinstance(fnames,(tuple,list)): fnames = [fnames]
    if not any('*' in x or '?' in x for x in fnames):
      fnames.sort()
      fnames = [list(v) for k,v in itertools.groupby(fnames,key=lambda x:x.rsplit('_',1)[0])]
    f = [read(fname) for fname in fnames]
    send['files'] = f
    if len(f)==1: f=f[0]
    send['f'] = f
    send(f)

@coroutine
def pltshow(send):
  plt.show()
  while True:
    send((yield))

@coroutine
def mkplt(send,*args,**kwargs):
  plt.clf()
  show = kwargs.pop('show',False)
  fig = plt.figure(*args,**kwargs)
  sub = fig.add_subplot(111)
  send['fig'],send['sub'] = fig,sub
  try:
    while True:
      send((yield))
      if show:
        plt.show()
  except GeneratorExit:
    plt.close(fig)

@coroutine
def mksubplots(send,grid=(1,2),*args,**kwargs):
  plt.clf()
  show = kwargs.pop('show',False)
  kwargs.setdefault('sharex',True)
  kwargs.setdefault('sharey',True)
  fig,grid = plt.subplots(grid[0],grid[1],*args,**kwargs)
  subs = itertools.cycle(grid.flat)
  send.update({'fig' : fig, 'subs' : subs, 'grid' : grid, 'sub' : subs.next()})
  try:
    while True:
      send((yield))
  except GeneratorExit:
    plt.close(fig)

@coroutine
def mkgrid(send,grid,figsize=None,**kwargs):
  import itertools
  plt.clf()
  c = {"nrows_ncols" : grid,
       "direction" : "column",
       "axes_pad" : 0.1,
       "add_all" : True,
       "label_mode" : "L",
       "share_all" : True,
       "cbar_pad" : 0.1}
       #"cbar_mode" : "single",
       #"cbar_location" : "right",
       #"cbar_size" : "2%",
       #"cbar_pad" : 0.2}
  c.update(kwargs)
  rc = c['nrows_ncols']
  fig = plt.figure(1,figsize)
  from mpl_toolkits.axes_grid1 import ImageGrid
  grid = ImageGrid(fig, 111, **c)
  subs = itertools.cycle(grid)
  send.update({'fig' : fig, 'subs' : subs, 'grid' : grid, 'sub' : subs.next()})
  try:
    while True:
      send((yield))
  except GeneratorExit:
    plt.close(fig)



@coroutine
def control(send,i=None):
  import sys
  d = None
  show = True
  fig = plt.figure()
  sub = fig.add_subplot(111)
  send['fig'],send['sub'] = fig,sub
  f = send.get('f')
  def select(i):
    if isinstance(f,tuple) or isinstance(f,list):
      for file in f:
        file.select(i)
    else:
      f.select(i)
  def mvframe(i):
    if isinstance(f,tuple) or isinstance(f,list):
      for file in f:
        file.select(file.frame+i)
    else:
      select(f.frame+i)
  def setframe(i):
    select(i)
  def press(event):
    #print 'press', event.key
    if(event.key in keys):
      keys[event.key]()
      try:
        send(d)
      except StopIteration:
        pass
      fig.canvas.draw()
  fig.canvas.mpl_connect('key_press_event', press)
  keys = dict(\
    right   = lambda: mvframe(1),
    left  = lambda: mvframe(-1),
    up    = lambda: mvframe(10),
    down  = lambda: mvframe(-10),
    home  = lambda: setframe(0),
    end   = lambda: setframe(-1),
#    a   = lambda: toggleAnimate(),
#    backspace = lambda: plotfn.rescale(),
#    pageup =
#    pagedown = 
    q   = lambda: sys.exit(0),
  )
  if i!=None:
    select(i)
  try:
    while True: 
      d = (yield)
      send(d)
      if show:
        plt.show()
  except GeneratorExit:
    plt.close(fig)

@coroutine
def append(send, l):
  while True:
    d = (yield)
    l.append(d)
    send(d)

@coroutine
def togglemax(send):
  C = True
  A = False
  B = False
  a = 0.0
  while True:
      d = (yield)
      d = d[500,:]       # this value needs to be adjusted
#      print d.shape
      # toggles between amax and amin
      if C:
        amax = np.amax(d)
        if amax > a:
          A = True
        if amax < a:
          B = True
        if (A and B) and (amax > a):
          A = B = False
          C = False
          amax = np.amin(d)
        a = amax
      else:
        amin = np.amin(d)
        if amin < a:
          A = True
        if amin > a:
          B = True
        if (A and B) and (amin < a):
          A = B = False
          C = True
          amin = np.amax(d)
        a = amin
      print a
      send(a);

@coroutine
def printmeta(send):
  while True:
    d = (yield)
    print send.data
    send(d)

@coroutine
def setmeta(send,d):
  send.update(d)
  while True:
    send((yield))

@coroutine
def consumemeta(send,name):
  while True:
    d = (yield)
    send[name] = d
    send(d)

# relevant physical quantities:

def vorticity(norm=lambda x,y: x):
  #return produce(("/timedisc/xvelocity","/timedisc/yvelocity")) | get(f) | only(1, numexpr("d+%f*f['/mesh/bhy']" % omega,f)) | curl(f)
  #return produce(("/timedisc/xvelocity","/timedisc/yvelocity")) | get(f) | only(1, numexpr("d-f['/mesh/bhy']**(-0.5)+%f*f['/mesh/bhy']" % omega,f)) | curl(f)
  return apply(lambda f: (f["/timedisc/xvelocity"],f["/timedisc/yvelocity"]+f.get('/mesh/radius',0.)*f['/config/mesh/omega'])) |  each(lambda: normalize2(norm=norm)) | curl()
  #return produce(("/timedisc/xvelocity","/timedisc/yvelocity")) | get(f) | only(1, numexpr("d-f['/mesh/bhy']**(-0.5)+%f*f['/mesh/bhy']" % omega,f)) | curl(f)
  #return apply(lambda f: (f['/timedisc/xvelocity'],f['/timedisc/yvelocity']-f['/mesh/bhy']**(-0.5)+f['/config/mesh/omega']*f['/mesh/bhy'])) | curl()

def kepvelocity(f):
  units = f['/config/physics/units']
  if units == 1:
    GN = 6.6742E-11
  elif units == 3:
    GN = 1.
  else:
    raise "unknown units"
  mass = f["/sources/grav/pmass/mass"]
  return np.sqrt(GN*mass/f['/mesh/radius'])

def kepvorticity():
  return apply(lambda f: (f["/timedisc/xvelocity"],f["/timedisc/yvelocity"]+f.get('/mesh/radius',0.)*f['/config/mesh/omega'] - kepvelocity(f)))  | curl()

def potvorticity(f,omega=0.):
    return apply(lambda d: (f['/timedisc/xvelocity'],f['/timedisc/yvelocity'])) | curl(f) | apply(lambda d: d + 2*omega) | apply(lambda d: d/f['/timedisc/density'])

def vortensity(f,omega=0.):
  return vorticity(f,omega) | numexpr("d/f['/timedisc/density']",f)

def specifictorque(f,omega=0.):
  #return numexpr("f['/mesh/bhy'] * (f['/timedisc/yvelocity'] + f['/mesh/bhy'] * %f)" % omega, f)
  return produce(("/timedisc/xvelocity","/timedisc/yvelocity")) | get(f) | Convert2Cartesian(f) | cross(f,(f['/mesh/bary_centers_x'],f['/mesh/bary_centers_y'],))

def mass(f):
  return numexpr("f['/timedisc/density'] * f['/mesh/volume']",f)

# e.g. control(f,fig) | specifictorqueflux(f,0.1) | imshow(fig,sub,aspect=4)
def specifictorqueflux(f,csiso=None,omega=0.):
  if csiso==None:
    return produce((1,2)) | split(\
      specifictorque(f,omega) | grad(f), \
      produce(('/timedisc/xvelocity','/timedisc/yvelocity')) | get(f)\
      ) | inner() | numexpr("d + f['/mesh/bhy']*f['/timedisc/pressure']", f)
  else:
    return produce((1,2)) | split(\
      produce((1,2)) | split(\
        specifictorque(f,omega) | grad(f), \
        produce(('/timedisc/xvelocity','/timedisc/yvelocity')) | get(f)\
        ) | inner(),\
      produce((0.,2)) | only(1,numexpr("f['/mesh/bhy']*f['/timedisc/density']*%f**2" % csiso, f)) | div(f) \
      ) | apply(lambda a,b: a + b)

# only for cartesian coords. for now
def angularmomentum(f):
  return produce(("/timedisc/xvelocity","/timedisc/yvelocity")) | get(f) | Convert2Cartesian(f) | cross(f,(f['/mesh/bary_centers_x'],f['/mesh/bary_centers_y'],))  \
    | apply(lambda d: d * f['/mesh/volume'] * f['/timedisc/density'])

# Ready to use plots:

line = lambda v: produce(('/mesh/bary_centers_x',v)) | control(f,fig) | get(f) | plot(f,fig,sub)

def im(v = "/timedisc/density",normalize=False,**kwargs):
  if normalize:
    return produce(v) | control(f,fig) | get(f) | normalize2() | imshow(fig,sub,**kwargs)
  else:
    return produce(v) | control(f,fig) | get(f) | imshow(fig,sub,**kwargs)

def angularmomentumconservation():
  return control(f,fig) | timeseries(f) | angularmomentum(f) | sum() | makelist(len(f)) | apply(lambda d: d,) | normalize() | numexpr('abs(d-1)') | apply(lambda d: d.clip(1.E-16,1.)) | printer() | log10() | plot(f,fig,sub)

def cellaspectratio():
  c = f['/mesh/corners']
  return control(f,fig) | apply(lambda d: (c[:,0,1,0]-c[:,0,0,0])/(c[:,0,3,1]-c[:,0,2,1])) | plot(f,fig,sub)

def gettemp(f):
  MU = 6.02E-4 #kg/mol, 75/25 Massenmischung H/He, vollständig ionisiert
  RG = 8.3144598 #J/(mol*K) allgemeine gaskonstante
  mH = 1.67E-27 #kg, mass of H-atom
  kB = 1.381E-23 #J/K, boltzmann constant
  fac = MU/RG
  #fac = 0.5 * mH/kB
  if '/timedisc/pressure' in f.keys():
    return f['/timedisc/pressure']/f['/timedisc/density']*fac
  else:
    return fac*f['/physics/bccsound']**2

def temperature():
  return apply(lambda f: gettemp(f))

def massspectrum(logscale=True,omega=0.):
  @coroutine
  def sort(send,f,logscale):
    while True:
      m, l = (yield)
      m, l = np.ravel(m),np.ravel(l)
      s = np.argsort(l)
      l = l[s]
      m = np.cumsum(m[s])
      if logscale:
        l = np.log10(l)
        m = np.log10(m)
      send((l,m))
  return produce((1,2)) | control(f,fig) | split(specifictorque(f,omega),mass(f)) | sort(f,logscale) | plot(f,fig,sub)

omega = lambda f: f['/timedisc/yvelocity']/f['/mesh/radius'] + f['/config/mesh/omega']

epicyclicfrequency = lambda f: np.sqrt(2.*omega(f)/f['/mesh/radius']*diff(f,f['/mesh/radius']**2*omega(f),0)/f['/mesh/bhx'])

#toomre stable for >1
def toomre(kappa=omega):
  return apply(lambda f: f['/physics/bccsound']*np.abs(kappa(f))/(3.14159*6.67384E-11*f['/timedisc/density']))

#def toomre(omega=1.,gamma=2.):
#  gn = 6.67408e-11
#  return apply(lambda f: (np.sqrt(gamma*f['/timedisc/pressure'])*omega)/(np.pi*gn*(f['/timedisc/density'])**(1.5)))
#
#def toomre2(f,omega=1.,gamma=2.):
#  gn = 6.67408e-11
#  cs = gn*10*np.pi
#  return apply(lambda d: (cs*omega)/(np.pi*gn*f['/timedisc/density']))

def kinetic_energy(f,omega=1.):
  gn = 6.67408e-11
  return apply(lambda d: (0.5*f['/timedisc/density']*(f['/timedisc/xvelocity']**2+(f['/timedisc/yvelocity']+1.5*omega*f['/mesh/bary_centers'][:,:,0])**2)))

def gravitational_energy(f):
  return apply(lambda d: f['/timedisc/density']*f['/sources/grav/self/phi'])

def thermal_energy(f):
  return apply(lambda d: 1.5*f['/timedisc/pressure'])

def cartvelocity():
  return apply(lambda f: (f['/timedisc/xvelocity'],f['/timedisc/yvelocity']+f['/mesh/radius']*f['/config/mesh/omega'])) | Convert2Cartesian()

def mach():
  return apply(lambda f: f['/timedisc/xvelocity']/f['/physics/bccsound'],(f['/timedisc/yvelocity']+f['/mesh/radius']*f['/config/mesh/omega'])/f['/physics/omega'])

def overR(p,i=0):
  return apply(lambda f: (f['/mesh/radius'][:,i],p(f)[:,i]))



@coroutine
def getvars(send,l):
  f = send.get('f')
  # Remove not available keys
  l = [n for n in l if n in f]
  while True:
    (yield)
    send(np.array([f[i].T for i in l]).T)

pvars = lambda: getvars(['/timedisc/%s' % s for s in ['density','xvelocity','yvelocity','pressure']])
cvars = lambda: getvars(['/timedisc/%s' % s for s in ['density','xmomentum','ymomentum','energy']])

@coroutine
def L1(send):
  f = send.get('f')
  inum,jnum = f['/config/mesh/inum'],f['/config/mesh/jnum']
  d0 = None
  while True:
    d = (yield)
    if d0 is None:
      d0 = d
    dU = np.abs(d-d0)
    sumdU = np.sum(dU,(0,1))
    L1 = np.sqrt(np.sum(sumdU**2)) / (inum*jnum)
    send(L1)
    

if __name__=="__main__":
  parser = argparse.ArgumentParser()
  group = parser.add_mutually_exclusive_group()
  group.add_argument('-p','--plot',help='enable 2d plotting',action="store_true")
  group.add_argument('-P','--plot3d',help='enable 3d plotting',action="store_true")
  parser.add_argument('-i','--import',nargs='?',help='import extra coroutines')
  parser.add_argument('pipe',metavar='pipe',help='pipe')
  parser.add_argument('filename',nargs='+',metavar='filename',help='filename')
  args = parser.parse_args()
  #f = [read(fname) for fname in args.filename]
  #if len(f) == 1: f = f[0]
  d = vars(args)
  if d['import']: execfile(d['import'])
  #if args.plot:
    #fig = plt.figure()
    #sub = fig.add_subplot(111)
  elif args.plot3d:
    fig = plt.figure()
    sub = fig.add_subplot(111, projection='3d')
  p = eval(args.pipe)
  #print p,type(p),args.pipe
  p.send(args.filename)
  #if args.plot:
  #  plt.show()

  # example:
  # ~/fosite/tools/plot/coplot.py -p "line('/timedisc/density')" pringle.xmf
