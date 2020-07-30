#!/usr/bin/env python3

import matplotlib.pyplot as plt
import math

def plot(ax,fname,iname):
    x=[]
    y=[]
    with open(iname+'_'+fname+'.dat','r') as fd:
        l=fd.readline()
        y=list(map(float,l.split()))
        x=list(range(len(y)))
    ax.loglog(x[1:],y[1:],label=fname)
    
def plot_poweri4(ax,fname,label):
    y=[0]
    with open(fname,'r') as fd:
        l=fd.readline()
        for l in fd:
            parse = l.split()
            y.append(float(parse[2])/(320**3))
    x=list(range(len(y)))
    ax.loglog(x[1:],y[1:],label=label)
    

def plot_pw(name,funlist,subname):
    fig=plt.figure()
    ax=fig.add_subplot()
    for f in funlist:
        plot(ax,f,name)
    plot_poweri4(ax,"newton_power006",'Poweri4')
    ax.legend()
    plt.xlabel("mode number")
    plt.ylabel("P_k")
    plt.show()
    #plt.savefig(name+'_'+subname+'.png')



if __name__=="__main__":
    plot_pw('pw-gevolution',['ngp','cic','tsc','pcs','gauss'],'all')

