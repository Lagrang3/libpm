#!/usr/bin/env python3

import matplotlib.pyplot as plt
import math,numpy

def plot(ax,fname,iname,base):
    x=[]
    y=[]
    with open(iname+'_'+fname+'.dat','r') as fd:
        l=fd.readline()
        y=list(map(float,l.split()))
        x=list(range(len(y)))
    y = numpy.array(y)
    y = (y-base)/base
    ax.plot(x[1:],y[1:],label=fname)
    
def get_poweri4(ax,fname):
    y=[]
    with open(fname,'r') as fd:
        l=fd.readline()
        for l in fd:
            parse = l.split()
            y.append(float(parse[2])/(320**3))
    return numpy.array(y)
    

def plot_pw(name,funlist,subname):
    fig=plt.figure()
    ax=fig.add_subplot()
    base = get_poweri4(ax,"newton_power006")
    for f in funlist:
        plot(ax,f,name,base)
    
    ax.legend()
    plt.xlabel("mode number")
    plt.ylabel("P_k")
    plt.show()
    #plt.savefig(name+'_'+subname+'.png')



if __name__=="__main__":
    plot_pw('pw-gevolution',['ngp','cic','tsc','pcs','gauss'],'all')

