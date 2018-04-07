#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Jan 28, 2018

@author: drews
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import random as rd
import os.path

def readDataOrder(t,Yn,sigY,order,plot=False):
    rc('text',usetex=True)
    rc('font',family='sanserif')
    ts = []
#    print("order: {}".format(order))
    for o in range(order):
        ts.append([x**(o+1) for x in t])
    F = np.vstack((np.ones([1,len(t)]),ts))
#    print("T: {}".format(t))
#    print("F: {}".format(F))
    gamma = np.divide(F,sigY)
#    print("F/sigy: {}".format(gamma))
    alpha = np.dot(gamma,np.transpose(gamma))
    delta = np.divide(Yn,sigY)
#    delta[np.isnan(delta)] = np.nanmean(delta)
    beta = np.dot(delta,np.transpose(gamma))
    epsilon = np.linalg.inv(alpha)
    Afit = np.dot(beta,epsilon)
    Yfit = np.dot(Afit,F)

    
    if plot:
        plt.figure()
        plt.errorbar(t,Yn,sigY,fmt='b.')
        plt.plot(t,Yfit,'r-')
        plt.show()
        print("order: {}".format(order))
        print("Afit: {}".format(Afit))        
    return Afit,Yfit
    
def readData(*args):
    rc('text',usetex=True)
    rc('font',family='sanserif')
    
    if len(args) <= 1:
        if len(args) == 1:
            filename = args[0]
        else:
            filename = "/Volumes/AS0216$/Library/Comp Methods/LLSData.txt"
        print("Filename is: {}".format(filename))
        if not os.path.isfile(filename):
            print("No file found at {}".format(filename))
        else:
            with open(filename,'r') as F:
                data = F.read().splitlines()
                t,Yn,sigY = [],[],[]
                for line in data[1:]:
                    line = line.split()
                    t.append(float(line[0]))
                    Yn.append(float(line[1]))
                    sigY.append(float(line[2]))
                t = np.array(t)
    else:
        t,Yn,sigY = args[:]
    order = 1
    while True:
        ts = []
        for o in range(order):
            ts.append([x**(o+1) for x in t])
        
        F = np.vstack((np.ones([1,len(t)]),ts))
        gamma = np.divide(F,sigY)
        alpha = np.dot(gamma,np.transpose(gamma))
        delta = np.divide(Yn,sigY)
        beta = np.dot(delta,np.transpose(gamma))
        epsilon = np.linalg.inv(alpha)
        Afit = np.dot(beta,epsilon)
        sAfit = np.sqrt(np.diagonal(epsilon))
        Yfit = np.dot(Afit,F)
        
        residual = [abs(Yfit[i] - Yn[i]) for i in range(len(t))]
        avgDev = np.sqrt(np.mean([x**2 for x in residual]))

        if avgDev < np.mean(sigY) or order > 9:
            break
        else:
            order += 1
    
    
    plt.figure()
    plt.errorbar(t,Yn,sigY,fmt='b.')
    plt.plot(t,Yfit,'r-')
    plt.show()
    print("order: {}".format(order))
    print("deviance: {}".format(avgDev))
    print("Afit: {}".format(Afit))        
    print("sAfit: {}".format(sAfit))
    print("sAfitM: {}".format(np.mean(sAfit)))
    return Afit
'''
I couldn't finish this as I didn't know how to compute reasonable 
noise for a given function.
'''    
def createPoly():
    numElements = 100
    
    A,ts = [],[]
    roots = []
    order = rd.randint(2,5)
#    print("order: {}".format(order))
#    print("step: {}".format(step))
    for o in range(order):
        A.append(10*rd.random())
        roots.append(rd.randint(-10,10))
    roots.append(100*rd.random() - 50)
    A.append(10*rd.random())

    print("roots: " + str(roots))
#    A = np.polynomial.polynomial.polyfromroots(roots)
    print("coef: {}".format(A))
    start = -10
    end = 10
#    start,end = max(a,b),min(a,b)
    step = (end-start)/numElements
    print("order: {}".format(order))
    print("start={} end={} step={}".format(start,end,step))
    t = np.arange(int(start),int(end),step)
    
    for o in range(order):
        ts.append([x**(o+1) for x in t])
        
    F = np.vstack((np.ones([1,len(t)]),ts))
    Y = np.dot(A,F)
    sigy = 0.75*np.ones(len(Y))*step/4
    noise = rd.random()*np.ones(len(Y))*10
    Yn = Y + noise
    readData(False,t,Yn,sigy)
    

if __name__ == '__main__':
    # readData() => runs the example
    # readData(filename) => runs the given file
    # readData(t,Yn,sigY) => runs the given data
    readData()
