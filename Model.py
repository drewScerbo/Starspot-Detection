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

class Model:
    
    def getKS(self,array):
        from scipy.stats import kstest
        return(kstest(array,'norm'))
    
    def getData(self):
        client = kplr.API()
#        return [client.koi(340.01),client.planet('2b')]
        return [client.koi(340.01)]
    
    def convolve(self,lc,model):
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
        return transits
    
    def normalize(self, koiLC,time0,duration,period,numTransits,quarter,plot=False):
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
            plt.xlabel("Time [days]")
            plt.ylabel("Flux")
            plt.title("Quarter " + str(quarter))
            Ymin, Ymax = min(Flux), max(Flux)
            for t in transits:
                plt.plot([t,t],[Ymin,Ymax],'r-')
                idx = (np.abs(np.array(Time)-t)).argmin()
                plt.plot([Time[idx-int(duration)],Time[idx+int(duration)]],\
                         [Flux[idx],Flux[idx]],'r-')
            plt.show()
        
        normedLCs, normedTs, transitT = [],[],[]
        transits = [transits[0]]
        for i in range(len(transits)):
            t = transits[i]
            transitT.append(t)
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
                plt.legend(["In-transit Flux","Out-transit Flux","LLS fit of Out-transit"])
                plt.xlabel("Normalized Time [days]")
                plt.ylabel("Flux")
                plt.show()
                
                plt.figure()
                plt.plot(times,diffLC,'b.')
                plt.xlabel("Normalized Time [days]")
                plt.ylabel("Normalized Flux")
                plt.show()
                
                

        return normedLCs,normedTs,transitT
    
    def makeModel(self,tt,period,incl,t,dor,ror,ldm_coeff1,ldm_coeff2):
        """
        -occultquad
        INPUTS:
            z -- sequence of positional offset values
            p0 -- planet/star radius ratio
            gamma -- two-sequence.
            koi_ldm_coeff1
            koi_ldm_coeff2
           
        OUTPUTS:
            the function based at times z
        
        -t2z
        I have modified this function so that the midpoint of every transit is 0
        :INPUTS:
        tt --  scalar. transit ephemeris
        per --  scalar. planetary orbital period (in days)
        inc -- scalar. orbital inclination (in degrees)
        hjd -- scalar or array of times, typically heliocentric or
               barycentric julian date.
        ars -- scalar.  ratio a/Rs,  orbital semimajor axis over stellar radius
        koi_dor
        OUPUTS:
            
        """
        # use t2z to make z
        z = t2z(tt,period,incl,t,dor)
        return occultquad(z,ror,[ldm_coeff1,ldm_coeff2])
    
    def applyModel(self,LCs,times,transits):
        """
        take in a lightcurves and times for each transit
        take in planet period
        
        find places in
        
        """
        pass