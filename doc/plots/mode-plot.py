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
    eps=1e-8
    if math.fabs(den)<eps:
        return 2*kmax+1
    return num/den 

class sample_func:
    def __init__(self,f):
        self.name = f.__name__
        self.W = f
    def __call__(self,x,N):
        return self.W(x,N)
        
def pbc(x):
    while x>0.5:
        x-=1
    while x<-0.5:
        x+=1
    return x
    
@sample_func
def gaussian(x,N):
    x=pbc(x)
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
def low_pass(x,N):
    x=pbc(x)
    kmax = kNyquist(N)
    den=math.sin(math.pi*x)
    num=math.sin(math.pi*x*(2*kmax+1))
    eps=1e-8
    if math.fabs(den)<eps:
        return 2*kmax+1
    return num/den 

@sample_func
def Sinc(x,N):
    x=pbc(x)
    den=math.sin(math.pi*x)
    num=math.sin(math.pi*x*N)
    eps=1e-8
    if math.fabs(den)<eps:
        return N
    return num/den 

@sample_func
def NGP(x,N):
    x=pbc(x)
    l=1/N
    xmod = math.fabs(x)/l
    eps = 1e-14
    if xmod < .5:
        return 1/l
    if xmod < .5+eps:
        return .5/l
    return 0

## sin filter
@sample_func
def SIN(x,N):
    x=pbc(x)
    l=1/N
    T = 10*l
    xmod = math.fabs(x)/T
    if xmod < 0.25:
        return math.pi*math.cos(2*math.pi*xmod)/T
    return 0

@sample_func
def CIC(x,N):
    x=pbc(x)
    l=1/N
    xmod = math.fabs(x)/l
    if xmod < 1.:
        return (1-xmod)/l
    return 0

@sample_func
def TSC(x,N):
    x=pbc(x)
    l=1/N
    xmod = math.fabs(x)/l
    if xmod < .5:
        return (.75-xmod*xmod)/l
    if xmod < 1.5:
        return (0.5*(1.5-xmod)**2)/l
    return 0

@sample_func
def PCS(x,N):
    x=pbc(x)
    l=1/N
    xmod = math.fabs(x)/l
    if xmod < 1.:
        return (1/6 * (4-6*xmod*xmod + 3*xmod**3))/l
    if xmod < 2:
        return (1/6 * (2-xmod)**3)/l
    return 0

def W2(x,N,f):
    kcut = kNyquist(N)
    return spi.quad(lambda xp : low_pass_filter(x-xp,kcut)*f(xp,N) ,0.0,1.0)[0]

@sample_func
def Wgauss(x,N):
    return W2(x,N,gaussian)


def show_delta(N):
    kcut=kNyquist(N)
    fig=plt.figure()
    ax=fig.add_subplot()
    X = np.linspace(0.,1.,1000)
    Y = [ low_pass_filter(x,kcut) for x in X ]
    ax.plot(X,Y)
    ax.set_xlabel('x')
    ax.set_ylabel('amplitude')
    plt.savefig('delta-W.pdf')
    

def show_W(N,W):
    fig=plt.figure()
    ax=fig.add_subplot()
    X = np.linspace(-0.5,0.5,1000)
    Y = [ W(x,N) for x in X ]
    ax.plot(X,Y,label=W.name)
    ax.set_xlabel('x')
    ax.set_ylabel('amplitude')
    ax.set_title('%s (N=%d)' % (W.name,N))
    #plt.show()
    plt.savefig('%s-real-%d.pdf' % (W.name,N))

def base_re(k,x):
    return math.cos(2*math.pi * x * k)
def base_im(k,x):
    return math.sin(2*math.pi * x * k)

def prod(k,W,N):
    v_re,_  = spi.quad(lambda xp : base_re(k,-xp)*W(xp,N) ,-0.5,0.5)
    v_im,_  = spi.quad(lambda xp : base_im(k,-xp)*W(xp,N) ,-0.5,0.5)
    return math.sqrt(v_re**2 + v_im**2)

def DFT(W,N):
    a=0.
    b=1.
    L = b-a
    vec = np.array([ W(a + L*i/N,N)  for i in range(N) ])
    vec_t = np.fft.fft(vec)
    vec_t = np.absolute(vec_t)/N
    kmax = (N-1)//2
    res = np.zeros(kmax+1)
    res[0]=vec_t[0]
    for i in range(1,kmax+1):
        res[i] = (vec_t[i]+vec_t[N-i-1])/2 
    return res
    
def show_innerP(N,W):
    fig=plt.figure()
    ax=fig.add_subplot()
    K = list(range(0,N))
    Y = [ (prod(k,W,N)+prod(-k,W,N))/2 for k in K ]
    ax.plot(K,Y,'-o',label='Fourier modes')
    Yd = DFT(W,N)
    ax.plot(list(range(0,len(Yd))),Yd,'-o',label='DFT modes')
    ax.legend()
    ax.set_title('%s (N=%d)' % (W.name,N))
    ax.set_xlabel('k')
    ax.set_ylabel('amplitude')
    #plt.show()
    plt.savefig('%s-modes-%d.pdf' % (W.name,N))

def show_interpolate(N,W,f):
    a=0
    b=1
    L = b-a
    pts = [ a+i*L/N for i in range(N)  ]
    vec = [ f(x) for x in pts ]
    
    def f2(x):
        s=0
        for i in range(N):
            s+=vec[i] * W(x-i*L/N-a ,N )
        return L*s/N
    
    fig=plt.figure()
    ax=fig.add_subplot()
    X = np.linspace(a,b,10000)
    Y = [ f(x) for x in X ]
    Y2 = [ f2(x) for x in X  ]
    ax.plot(X,Y,'-',label='true function')
    ax.plot(X,Y2,'-',label='interpolated')
    ax.plot(pts,vec,'o',label='grid points')
    ax.legend()
    ax.set_title('%s (N=%d)' % (W.name,N))
    ax.set_xlabel('x')
    ax.set_ylabel('amplitude')
    #plt.show()
    plt.savefig('%s-interpolation-%d.pdf' % (W.name,N))

if __name__ == "__main__":
    n = 20
    for w in [NGP,CIC,TSC,PCS,gaussian,Sinc,low_pass]:
        show_innerP(n,w)
        #show_W(n,w)
    n = 21
    for w in [NGP,CIC,TSC,PCS,gaussian,Sinc,low_pass]:
        show_innerP(n,w)
        #show_W(n,w)
    def f(x):
        return math.sin(2*math.pi*9*x) + math.sin(2*math.pi*5*x)
    n = 20
    for w in [Sinc,low_pass]:
        show_interpolate(n,w,f)
    n = 21
    for w in [Sinc,low_pass]:
        show_interpolate(n,w,f)
