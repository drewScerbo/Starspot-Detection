#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 12:27:38 2018

@author: as0216
"""
from scipy.stats import kstest
from transit import t2z
import numpy as np
import kplr
from LLSDataReader import readDataOrder
import matplotlib.pyplot as plt

class Model:
        
    def getKS(self,array):
        return(kstest(array,'norm'))
    
    def getData(self):
        client = kplr.API()
        #return [client.planet('2b'),
        return [client.koi(340.01)]
    
    def make_LC_model(self,koi):
        # use Mandel and Agol 2002 equations
        # transit.occultquad()
        pass
    
    def convolve(self,lc,model):
        return np.convolve(lc,model)
    
    def findTransits(self,Time,Flux,num_transits,time0,period):
        transits = []
        
        ### this finds transits times ###
        i, c = 0, 0
        up = True
        while c < num_transits:
            if up:
                t = time0 + i*period
                if t < max(Time):
                    transits.append(t)
                    i += 1
                    c += 1
                else:
                    up = False
                    i = 1
            else:
                t = time0 - i*period
                if t > min(Time):
                    transits.append(t)
                    c += 1
                    i += 1
                else:
                    break
        return transits
    
    def normalize(self, koiLC,time0,duration,period,numTransits,plot=False):
        """
        normalize():
            find peak flux of lightcurve
            fit out transit with 2nd degree poly
            subtract data by this poly to get diffLightcurve
            add peak flux value to diffLightcurve
            divide diffLightcuve by peak flux value
            return diffLightcurve
        """
        time, flux = [], []
        with koiLC.open() as f:
            # The lightcurve data are in the first FITS HDU.
            hdu_data = f[1].data
            time.append(hdu_data["time"])
            flux.append(hdu_data["sap_flux"])
        
        Flux,Time = [],[]
        for i in range(len(time)):
            Flux.extend(flux[i])
            Time.extend(time[i])
            
        # remove NaN's
        b = np.array([i for i in range(len(Flux)) if str(Flux[i]) == 'nan'])
        Time = [Time[i] for i in range(len(Time)) if i not in b]
        Flux = [Flux[i] for i in range(len(Flux)) if i not in b]
        
        transits = self.findTransits(Time,Flux,numTransits,time0,period)
        if plot:
            plt.figure()
            plt.plot(Time,Flux,'g.')
            Ymin, Ymax = min(Flux), max(Flux)
            for t in transits:
                plt.plot([t,t],[Ymin,Ymax],'r-')
                idx = (np.abs(np.array(Time)-t)).argmin()
#                t1 = (np.abs(np.array(Time)-(t-duration))).argmin()
#                t2 = (np.abs(np.array(Time)-(t+duration))).argmin()
#                print(t)
#                print(t+duration)
#                print(t-duration)
#                print(duration)
                plt.plot([Time[idx-int(duration)],Time[idx+int(duration)]],\
                         [Flux[idx],Flux[idx]],'r-')
            plt.show()
        print(transits)
        for t in transits:
            print("at: " + str(t))
            idx = (np.abs(np.array(Time)-t)).argmin()
            end = idx+int(duration)
            start = idx-int(duration)
            
            endOut = idx+3*int(duration)
            startOut = idx-3*int(duration)
            outTransitT = Time[startOut:start] + Time[end+1:endOut+1]
            outTransitF = Flux[startOut:start] + Flux[end+1:endOut+1]
    
            sigY = [np.nanstd(outTransitF) for i in range(len(outTransitF))]
            A,Y = readDataOrder(outTransitT,outTransitF,sigY,2)
            
            ts = []
            for o in range(2):
                ts.append([x**(o+1) for x in Time[startOut:endOut+1]])
            F = np.vstack((np.ones([1,len(Time[startOut:endOut+1])]),ts))
            peakFlux = np.percentile(Flux,95)
            Yfit = np.dot(A,F)
            diffLC = Flux[startOut:endOut+1] - Yfit
            diffLC += peakFlux
            diffLC /= peakFlux
            if plot:
                plt.figure()
                intransitT = Time[start:end+1]
                intransitF = Flux[start:end+1]
                plt.plot(intransitT,intransitF,'go')
                plt.plot(outTransitT,outTransitF,'b.')
                plt.plot(outTransitT,Y,'k-')
                plt.show()
                
                plt.figure()
                plt.plot(Time[startOut:endOut+1],diffLC,'b.')
                plt.show()
            
        return [diffLC,Time[startOut:endOut+1]]
    
    def run_occultquad(self,z,p0,gamma):
        """
        INPUTS:
            z -- sequence of positional offset values
    
            p0 -- planet/star radius ratio
    
            gamma -- two-sequence.
               quadratic limb darkening coefficients.  (c1=c3=0; c2 =
               gamma[0] + 2*gamma[1], c4 = -gamma[1]).  If only a single
               gamma is used, then you're assuming linear limb-darkening.
           
        OUTPUTS:
            the function based at times z
        """
        return transit.occultquad(z,p0,gamma)
    
    t2z