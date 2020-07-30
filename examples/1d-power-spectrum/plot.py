#!/usr/bin/env python3

import matplotlib.pyplot as plt

def plot(ax,fname,iname):
    x=[]
    y=[]
    with open(iname+'_'+fname+'.dat','r') as fd:
        l=fd.readline()
        y=list(map(float,l.split()))
        x=list(range(len(y)))
    ax.plot(x[1:],y[1:],label=fname)

def plot_pw(name,funlist,subname):
    fig=plt.figure()
    ax=fig.add_subplot()
    for f in funlist:
        plot(ax,f,name)
    ax.legend()
    plt.show()
    #plt.savefig(name+'_'+subname+'.png')

if __name__=="__main__":
    plot_pw('pw-random',['exact','ngp','cic','tsc','pcs','gauss'],'all')

