#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys,os

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
    os.chdir(sys.argv[1])
    plot_interpolation('triangle',['exact','ngp'],'ngp')
    plot_interpolation('triangle',['exact','cic'],'cic')
    plot_interpolation('triangle',['exact','tsc'],'tsc')
    plot_interpolation('triangle',['exact','pcs',],'pcs')
    plot_interpolation('triangle',['exact','gaussian'],'gaussian')
    plot_interpolation('sine',['exact','ngp'],'ngp')
    plot_interpolation('sine',['exact','cic'],'cic')
    plot_interpolation('sine',['exact','tsc'],'tsc')
    plot_interpolation('sine',['exact','pcs',],'pcs')
    plot_interpolation('sine',['exact','gaussian'],'gaussian')
