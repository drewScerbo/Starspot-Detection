
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 12:27:38 2018

@author: as0216
"""
from transit import t2z,occultquad
import numpy as np
import kplr
from LLSDataReader import readDataOrder
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from scipy.signal import fftconvolve
import random

class Model:
    
    def getKS(self,arr1,arr2,plot=False):
        if plot:
            y1 = np.cumsum(arr1)
            y2 = np.cumsum(arr2)
            print(len(y1),len(y2))
            print(len(arr1),len(arr2))
            plt.figure()
            plt.plot(y1,'g--.')
            plt.plot(y2,'b--.')
            plt.title("CDF's of in and out transits")
            plt.show()
        with np.errstate(divide='ignore'):
            ks = ks_2samp(arr1,arr2)
        return ks
    
    def getData(self):
        client = kplr.API()
        # ror^2 >= 0.01 to start, try >= 0.008 next
#        planetKOIs = client.planets(koi_ror=">=0.01")
#        kois = client.kois(koi_ror=">=0.01")
#        print(len(kois))
        planetKOIs = [client.planet('2b')]
#        kois = [client.koi(340.01)]
        kois = []
        kois.extend([client.koi(s.koi_number) for s in planetKOIs])
        return kois
    
    def convolve(self,model):
        """
        
        how long is the shutter open for?
        
        Q:
            There are two kinds of cadence time for Kepler, 5 minutes and 30 
            minutes.
            I just know there is no shutter and the CCDs of Kepler can be 
            locked after exposure.
            
            But how to realize that? Remember we have to not affect the photometry.
            
            The other question is so far how many unique sources has Kepler 
            observed? 160,000?
        
        A:
            They don't stop. The Kepler CCDs read out the signals collected during
            the accumulated time of 6.02 seconds. The fixed read out time is 0.52 
            seconds. So each CCD gets one frame every (6.02 + 0.52) = 6.54 seconds.
            Then Kepler sums up every 9 frames (short cadence) and 270 frames 
            (long cadence). The time between two short cadences is 
            (6.54 x 9) = 58.9 seconds, and (6.54 x 270) = 1766 seconds between two 
            long cadences, but the exposure time is (6.02 x 9) = 54.2 seconds for a
            short cadence, and (6.02 x 270) = 1625 seconds for a long cadence.
        
        """
        k = []
        for i in range(1766*2):
            if i % 13== 0:
                k.append(1)
            else: 
                k.append(0)
        return fftconvolve(model,k,'same')
    
    def findTransits(self,koi,num_transits,time0,period,q=0):
        Flux,Time,Error = [],[],[]
        c = 0
        for lc in koi.get_light_curves(short_cadence=False):
            # lc => data for 1 quarter
#            if c > 0:
#                break
#            if q == c:
            with lc.open() as f:
                
                # The lightcurve data are in the first FITS HDU.
                hdu_data = f[1].data
                Time.extend(hdu_data["time"])
                Flux.extend(hdu_data["sap_flux"])
                Error.extend(hdu_data["sap_flux_err"])
        
        c +=1
        # remove NaN's
        b = np.array([i for i in range(len(Flux)) if str(Flux[i]) == 'nan'])
        Time = [Time[i] for i in range(len(Time)) if i not in b]
        Error = [Error[i] for i in range(len(Error)) if i not in b]
        Flux = [Flux[i] for i in range(len(Flux)) if i not in b]
        
        transits = []
        
        ### this finds transits times ###
        c = 0
        while c < num_transits:
            t = time0 + c*period
            if t < max(Time) and t > min(Time):
                transits.append(t)
            else:
                t = time0 - c*period
                if t < max(Time) and t > min(Time):
                    transits.append(t)
            c += 1
        return Time,Flux,Error,transits
    
    def make_fakes(self,Flux,Time,model,transits,duration):
        """
        num_fakes = number of transits (actually number of sets of number of transits)
        randomly pick until we have that many
        """
        num = int(len(model)//2)
        idx,i = random.randint(num,len(Time) - 1 - num),0
        LCs,times = [],[]
        while i < len(transits):
            if min(abs(Time[int(idx)] - transits)) <= duration/20:
                idx = random.randint(num,len(Time) - 1 - num)
                continue
            LCs.append(Flux[int(idx-num):int(idx+num)+1])
            times.append(Time[int(idx)])
            idx = random.randint(num,len(Time) - 1 - num)
            i += 1
        return LCs,times
        
    
    def normalize(self, Time,Flux,Error,t,duration,period,quarter=0,plot=False):
        """
        normalize():
            find peak flux of lightcurve
            fit out transit with 2nd degree poly
            subtract data by this poly to get diffLightcurve
            add peak flux value to diffLightcurve
            divide diffLightcuve by peak flux value
            return diffLightcurve
        """
        
        
        if plot:
            plt.figure()
            plt.plot(Time,Flux,'g.')
            plt.xlabel("Time [days]")
            plt.ylabel("Flux")
            plt.title("Quarter " + str(quarter))
            Ymin, Ymax = min(Flux), max(Flux)
            plt.plot([t,t],[Ymin,Ymax],'r-')
            idx = (np.abs(np.array(Time)-t)).argmin()
            plt.plot([Time[idx]-duration/48,Time[idx]+duration/48],\
                     [Flux[idx],Flux[idx]],'r-')
            plt.show()
        
        idx = (np.abs(np.array(Time)-t)).argmin()
        end = (np.abs(np.array(Time)-t-(duration/48))).argmin()
        start = (np.abs(np.array(Time)-t+(duration/48))).argmin()
        endOut = (np.abs(np.array(Time)-t-(duration/24))).argmin()
        startOut = (np.abs(np.array(Time)-t+(duration/24))).argmin()

        outTransitT = Time[startOut-1:start] + Time[end:endOut+1]
        outTransitT = np.subtract(outTransitT,t)
        outTransitF = Flux[startOut-1:start] + Flux[end:endOut+1]
        outTransitE = Error[startOut-1:start] + Error[end:endOut+1]

        if end - start < 1:
            return [],[]
        
        A,Y = readDataOrder(outTransitT,outTransitF,outTransitE,2)
        times = np.subtract(Time[startOut:endOut+1],t)
        F = np.vstack((np.ones([1,len(times)]),times,[x**(2) for x in times]))
        peakFlux = np.percentile(Flux,95)
        Yfit = np.dot(A,F)
        diffLC = Flux[startOut:endOut+1] - Yfit
        diffLC += peakFlux
        diffLC /= peakFlux
        
        intransitT = np.subtract(Time[start:end],t)

        if plot:
            print("transit time: {}".format(t))
            plt.figure()
            intransitF = Flux[start:end]
            plt.plot(intransitT,intransitF,'go')
            plt.plot(outTransitT,outTransitF,'b.')
            plt.plot(outTransitT,Y,'k-')
            plt.legend(["In-transit Flux","Out-transit Flux","LLS fit of Out-transit"])
            plt.xlabel("Normalized Time [days]")
            plt.ylabel("Flux")
            plt.show()
            
            plt.figure()
            plt.plot(times,diffLC,'b.')
            plt.xlabel("Normalized Time [days]")
            plt.ylabel("Normalized Flux")
            plt.show()

        return diffLC,times
    
    def makeModel(self,tt,period,incl,t,dor,ror,ldm_coeff1,ldm_coeff2):
        """
        -t2z
        I have modified this function so that the midpoint of every transit is 0
        
        find lmdk from the star 
        """
        # use t2z to make z
        z = t2z(tt,period,incl,t,dor)
        return occultquad(z,ror,[ldm_coeff1,ldm_coeff2])
    
    def getINOUT(self,lc,times,duration,plot=False):
        start = (np.abs(np.array(times)+(duration/48))).argmin()-1
        end = (np.abs(np.array(times)-(duration/48))).argmin()
        lcOUT = np.array(lc[:start])
        lcOUT = np.append(lcOUT,lc[end+1:])
        timeOUT = times[:start]
        timeOUT = np.append(timeOUT,times[end+1:])
        if plot:
            plt.figure()
            plt.title("In-transit vs Out-transit")
            plt.plot(timeOUT,lcOUT,'b.')
            plt.plot(times[start:end],lc[start:end],'g.')
            plt.show()
        return lc[start:end],lcOUT,times[start:end],timeOUT
    
    def applyModel(self,lc,times,model,duration,plot=False):
        """
        returns INresiduals, OUTresiduals, INtime, OUTtime
        """
        if plot:
            plt.figure()
            plt.plot(times,model,'k-')
            plt.plot(times,lc,'b.')
            plt.show()
    
        resINs, resOUTs = [],[] 
        timeINs, timeOUTs = [],[]
        start = (np.abs(np.array(times)+(duration/48))).argmin()-1
        end = (np.abs(np.array(times)-(duration/48))).argmin()
        residual = lc - model
        residualOUT = np.array(residual[:start])
        residualOUT = np.append(residualOUT,residual[end+1:])
        timeOUT = times[:start]
        timeOUT = np.append(timeOUT,times[end+1:])
        if plot:
            plt.figure()
            plt.title("In-transit vs Out-transit")
            plt.plot(timeOUT,residualOUT,'b.')
            plt.plot(times[start:end],residual[start:end],'g.')
            plt.show()
        resINs.append(residual[start:end])
        resOUTs.append(residualOUT)
        timeINs.append(times[start:end])
        timeOUTs.append(timeOUT)
        return residual[start:end],residualOUT,times[start:end],timeOUT


