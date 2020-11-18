#!/usr/bin/env python3

import matplotlib.pyplot as plt

def plot(ax,shape,fname):
    ax.set_title(fname)
    x=[]
    y=[]
    with open(shape+'_'+fname+'.dat','r') as fd:
        for l in fd:
            a,b = map(lambda x: float(x),l.split())
            #print(a,b)
            x.append(a)
            y.append(b)
    ax.plot(x,y,label=fname)
    
def plot_pt(ax,shape,fname):
    x=[]
    y=[]
    with open(shape+'_'+fname+'.dat','r') as fd:
        for l in fd:
            a,b = map(lambda x: float(x),l.split())
            #print(a,b)
            x.append(a)
            y.append(b)
    ax.plot(x,y,'o',label=fname)
    

def plot_interpolation(shape,funlist,subname):
    fig=plt.figure()
    ax=fig.add_subplot()
    plot_pt(ax,shape,'points')
    for fun in funlist:
        plot(ax,shape,fun)
    ax.legend()
    #plt.show()
    plt.savefig(shape+'_'+subname+'.png')

if __name__=="__main__":
    plot_interpolation('triangle',['exact','NGP'],'ngp')
    plot_interpolation('triangle',['exact','CIC'],'cic')
    plot_interpolation('triangle',['exact','TSC'],'tsc')
    plot_interpolation('triangle',['exact','PCS',],'pcs')
    plot_interpolation('triangle',['exact','Gaussian'],'gaussian')
    plot_interpolation('sine',['exact','NGP'],'ngp')
    plot_interpolation('sine',['exact','CIC'],'cic')
    plot_interpolation('sine',['exact','TSC'],'tsc')
    plot_interpolation('sine',['exact','PCS',],'pcs')
    plot_interpolation('sine',['exact','Gaussian'],'gaussian')
