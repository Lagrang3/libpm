#!/usr/bin/env python3

import matplotlib.pyplot as plt

def plot(ax,iname,fname):
    ax.set_title(fname)
    x=[]
    y=[]
    with open(iname+'_'+fname+'.dat','r') as fd:
        for l in fd:
            a,b = map(lambda x: float(x),l.split())
            #print(a,b)
            x.append(a)
            y.append(b)
    ax.plot(x,y,label=iname)
    
def plot_pt(ax,iname,fname):
    ax.set_title(fname)
    x=[]
    y=[]
    with open(iname+'_'+fname+'.dat','r') as fd:
        for l in fd:
            a,b = map(lambda x: float(x),l.split())
            #print(a,b)
            x.append(a)
            y.append(b)
    ax.plot(x,y,'o',label=iname)
    

def plot_interpolation(name,funlist,subname):
    fig=plt.figure()
    ax=fig.add_subplot()
    plot_pt(ax,'pts',name)
    for f in funlist:
        plot(ax,f,name)
    ax.legend()
    plt.savefig(name+'_'+subname+'.png')

if __name__=="__main__":
    plot_interpolation('triangle',['exact','ngp','cic','tsc','pcs','gauss','low_pass'],'all')
    plot_interpolation('triangle',['exact','ngp'],'ngp')
    plot_interpolation('triangle',['exact','cic'],'cic')
    plot_interpolation('triangle',['exact','tsc'],'tsc')
    plot_interpolation('triangle',['exact','pcs',],'pcs')
    plot_interpolation('triangle',['exact','gauss'],'gauss')
    plot_interpolation('triangle',['exact','low_pass'],'low_pass')
    plot_interpolation('sine_5_6',['exact','ngp'],'ngp')
    plot_interpolation('sine_5_6',['exact','cic'],'cic')
    plot_interpolation('sine_5_6',['exact','tsc'],'tsc')
    plot_interpolation('sine_5_6',['exact','pcs',],'pcs')
    plot_interpolation('sine_5_6',['exact','gauss'],'gauss')
    plot_interpolation('sine_5_6',['exact','low_pass'],'low_pass')
