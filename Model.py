#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 12:27:38 2018

@author: as0216
"""
from scipy.stats import kstest
from transit import t2z,occultquad
import numpy as np
import kplr
from LLSDataReader import readDataOrder
import matplotlib.pyplot as plt

class Model:
        
    def getKS(self,array):
        return(kstest(array,'norm'))
    
    def getData(self):
        client = kplr.API()
        return [client.planet('2b'),client.koi(340.01)]
#        return [client.koi(340.01)]
    
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
        while c < num_transits:
            t = time0 + i*period
            if t < max(Time) and t > min(Time):
                transits.append(t)
            t = time0 - i*period
            if t < max(Time) and t > min(Time):
                transits.append(t)
            i += 1
            c += 1
                
#        for i in range(len(transits) - 1):
#            # remove if f(t) == f(t+1)
#            idx = (np.abs(np.array(Time)-transits[i])).argmin()
#            idx2 = (np.abs(np.array(Time)-transits[i+1])).argmin()
#            if Flux[idx] == Flux[idx2]:
#                transits.remove(transits[i])
#            else:
#                break
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
#        time, flux = [], []
        Flux,Time,Error = [],[],[]
        with koiLC.open() as f:
            # The lightcurve data are in the first FITS HDU.
            hdu_data = f[1].data
            Time.extend(hdu_data["time"])
            Flux.extend(hdu_data["sap_flux"])
            Error.extend(hdu_data["sap_flux_err"])
        
        # remove NaN's
        b = np.array([i for i in range(len(Flux)) if str(Flux[i]) == 'nan'])
        Time = [Time[i] for i in range(len(Time)) if i not in b]
        Error = [Error[i] for i in range(len(Error)) if i not in b]
        Flux = [Flux[i] for i in range(len(Flux)) if i not in b]
        
        transits = self.findTransits(Time,Flux,numTransits,time0,period)
        if plot:
            plt.figure()
            plt.plot(Time,Flux,'g.')
            Ymin, Ymax = min(Flux), max(Flux)
            for t in transits:
                plt.plot([t,t],[Ymin,Ymax],'r-')
                idx = (np.abs(np.array(Time)-t)).argmin()
                plt.plot([Time[idx-int(duration)],Time[idx+int(duration)]],\
                         [Flux[idx],Flux[idx]],'r-')
            plt.show()
        
        normedLCs, normedTs = [],[]
        for i in range(len(transits)):
            t = transits[i]
            idx = (np.abs(np.array(Time)-t)).argmin()
            end = idx+int(duration)
            start = idx-int(duration)
            endOut = idx+3*int(duration)
            startOut = idx-3*int(duration)
            
            outTransitT = Time[startOut:start] + Time[end+1:endOut+1]
            outTransitT = np.subtract(outTransitT,t)
            outTransitF = Flux[startOut:start] + Flux[end+1:endOut+1]
            outTransitE = Error[startOut:start] + Error[end+1:endOut+1]
            
            A,Y = readDataOrder(outTransitT,outTransitF,outTransitE,2)
            times = np.subtract(Time[startOut:endOut+1],t)
#            ts = [times,[x**(2) for x in times]]
#            for o in range(2):
#                ts.append([x**(o+1) for x in times])
#            F = np.vstack((np.ones([1,len(times)]),ts))
            F = np.vstack((np.ones([1,len(times)]),times,[x**(2) for x in times]))
            peakFlux = np.percentile(Flux,95)
            Yfit = np.dot(A,F)
            diffLC = Flux[startOut:endOut+1] - Yfit
            diffLC += peakFlux
            diffLC /= peakFlux
            normedLCs.append(diffLC)
            normedTs.append(times)
            if plot:
                plt.figure()
                intransitT = np.subtract(Time[start:end+1],t)
                intransitF = Flux[start:end+1]
                plt.plot(intransitT,intransitF,'go')
                plt.plot(outTransitT,outTransitF,'b.')
                plt.plot(outTransitT,Y,'k-')
                plt.show()
                
                plt.figure()
                plt.plot(times,diffLC,'b.--')
                plt.show()
                
                

        return normedLCs,normedTs
    
    def makeModel(self,ror,ldm_coeff1,ldm_coeff2):
        """
        INPUTS:
            z -- sequence of positional offset values
    
            p0 -- planet/star radius ratio
    
            gamma -- two-sequence.
               quadratic limb darkening coefficients.  (c1=c3=0; c2 =
               gamma[0] + 2*gamma[1], c4 = -gamma[1]).  If only a single
               gamma is used, then you're assuming linear limb-darkening.
            koi_ldm_coeff1
            koi_ldm_coeff2
           
        OUTPUTS:
            the function based at times z
        """
        
        # use t2z to make z
        z = []
        return occultquad(z,ror,gamma)
    
    t2z