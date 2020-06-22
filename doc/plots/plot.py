#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as spi

def kNyquist(N):
    return math.floor((N-1)/2)

def low_pass_filter(x,kmax):
    den=math.sin(math.pi*x)
    num=math.sin(math.pi*x*(2*kmax+1))
    eps=0.001
    if math.fabs(den)<eps:
        return 2*kmax+1
    return num/den 

class sample_func:
    def __init__(self,f):
        self.name = f.__name__
        self.W = f
    def __call__(self,x,N):
        return self.W(x,N)

@sample_func
def gaussian(x,N):
    sigma=1/N
    expf= math.exp(-(x*x)/(2*sigma*sigma))
    #if expf<0.005:
    #    expf = 0
    return expf/math.sqrt(2*sigma*sigma*math.pi)

class gaussian_cut:
    def __init__(self,ncut,name):
        self.ncut = ncut
        self.name=name
    def __call__(self,x,N):
        sigma=1/N
        expf= math.exp(-(x*x)/(2*sigma*sigma))
        if x/sigma > self.ncut:
            expf = 0
        return expf/math.sqrt(2*sigma*sigma*math.pi)


@sample_func
def NGP(x,N):
    l=1/N
    xmod = math.fabs(x)/l
    if xmod < .5:
        return 1/l
    return 0

## sin filter
@sample_func
def SIN(x,N):
    l=1/N
    T = 10*l
    xmod = math.fabs(x)/T
    if xmod < 0.25:
        return math.pi*math.cos(2*math.pi*xmod)/T
    return 0

@sample_func
def CIC(x,N):
    l=1/N
    xmod = math.fabs(x)/l
    if xmod < 1.:
        return (1-xmod)/l
    return 0

@sample_func
def TSC(x,N):
    l=1/N
    xmod = math.fabs(x)/l
    if xmod < .5:
        return (.75-xmod*xmod)/l
    if xmod < 1.5:
        return (0.5*(1.5-xmod)**2)/l
    return 0

@sample_func
def PCS(x,N):
    l=1/N
    xmod = math.fabs(x)/l
    if xmod < 1.:
        return (1/6 * (4-6*xmod*xmod + 3*xmod**3))/l
    if xmod < 2:
        return (1/6 * (2-xmod)**3)/l
    return 0

def W2(x,N,f):
    kcut = kNyquist(N)
    return spi.quad(lambda xp : low_pass_filter(x-xp,kcut)*f(xp,N) ,-0.5,0.5)[0]

@sample_func
def Wgauss(x,N):
    return W2(x,N,gaussian)


def show_delta(N):
    kcut=kNyquist(N)
    fig=plt.figure()
    ax=fig.add_subplot()
    X = np.linspace(-0.5,0.5,1000)
    Y = [ low_pass_filter(x,kcut) for x in X ]
    ax.plot(X,Y)
    ax.set_xlabel('x')
    ax.set_ylabel('amplitude')
    plt.savefig('delta-W.pdf')
    

def show(N,shape):
    fig=plt.figure()
    ax=fig.add_subplot()
    X = np.linspace(-0.5,0.5,1000)
    Y = [ W2(x,N,shape) for x in X ]
    ax.plot(X,Y,label="low-pass convolution")
    Y = [ shape(x,N) for x in X ]
    ax.plot(X,Y,label=shape.name)
    ax.legend()
    ax.set_xlabel('x')
    ax.set_ylabel('amplitude')
    plt.savefig('%s-W.pdf' % shape.name)

def pk(N,kcut,shape):
    X = list(range(1,kcut))
    Y1 = [ spi.quad( lambda x: shape(x,N)*math.cos(-2*math.pi*x*k) ,-0.5,0.5)[0] for k in X  ]
    Y2 = [ spi.quad( lambda x: shape(x,N)*math.sin(-2*math.pi*x*k) ,-0.5,0.5)[0] for k in X  ]
    Pk = [ y[0]**2 + y[1]**2 for y in zip(Y1,Y2)  ]
    return X,Pk

def plotpk(N):
    fig=plt.figure()
    ax=fig.add_subplot()
    for f in [NGP,CIC,TSC,PCS,gaussian]:
        X,Pk = pk(N,N,f)
        ax.semilogy(X,Pk,'o-',label=f.name)
        #ax.plot(X,Pk,'o-',label=f.name)
    knyquist = kNyquist(N)
    ax.vlines([knyquist],0,1,transform=ax.get_xaxis_transform(),color='black',label='k_nyquist')
    ax.set_title("Power spectrum")
    ax.legend()
    ax.set_xlabel('k (modes)')
    ax.set_ylabel('Pk')
    plt.savefig('pk.pdf')

def showall(N):
    for f in [NGP,CIC,TSC,PCS,gaussian]:
        show(N,f) 
    

def plotpk_single(N,fun,fname):
    fig=plt.figure()
    ax=fig.add_subplot()
    X,Pk = pk(N,N,fun)
    ax.semilogy(X,Pk,'o-',label=fun.name)
    knyquist = kNyquist(N)
    ax.vlines([knyquist],0,1,transform=ax.get_xaxis_transform(),color='black',label='k_nyquist')
    ax.set_title("Power spectrum")
    ax.set_xlabel('k (modes)')
    ax.set_ylabel('Pk')
    plt.savefig(fname)

def show_single(N,shape,fname):
    fig=plt.figure()
    ax=fig.add_subplot()
    X = np.linspace(-0.5,0.5,1000)
    Y = [ W2(x,N,shape) for x in X ]
    ax.plot(X,Y,label="low-pass convolution")
    Y = [ shape(x,N) for x in X ]
    ax.plot(X,Y,label=shape.name)
    ax.legend()
    ax.set_xlabel('x')
    ax.set_ylabel('amplitude')
    plt.savefig(fname)

if __name__ == "__main__":
    show_delta(64)
    showall(64)
    plotpk(64)
    plotpk_single(64,gaussian_cut(4,'gaussian w. cut-off'),'pk-gauss-cut.pdf')
    show_single(64,gaussian_cut(4,'gaussian w. cut-off'),'wlow-gauss-cut.pdf')
