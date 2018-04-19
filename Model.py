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
        plan = client.planet('2b')
2        
        print(vars(plan))
        
        # get star of the planet
#        return []
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
                plt.plot([Time[idx]-duration/48,Time[idx]+duration/24],\
                         [Flux[idx],Flux[idx]],'r-')
            plt.show()
        
        normedLCs,normedTs,transitT = [],[],[]

        for i in range(len(transits)):
            t = transits[i]
            transitT.append(t)
            
            idx = (np.abs(np.array(Time)-t)).argmin()
            end = (np.abs(np.array(Time)-t-(duration/48))).argmin()
            start = (np.abs(np.array(Time)-t+(duration/48))).argmin()
            endOut = (np.abs(np.array(Time)-t-(duration/24))).argmin()
            startOut = (np.abs(np.array(Time)-t+(duration/24))).argmin()

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
            intransitT = np.subtract(Time[start:end+1],t)
            normedTs.append(times)
            if plot:
                print("transit time: {}".format(t))
                plt.figure()
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
        -t2z
        I have modified this function so that the midpoint of every transit is 0
        
        
        find lmdk from the star 
        """
        
        # use t2z to make z
        z = t2z(tt,period,incl,t,dor)
        return occultquad(z,ror,[ldm_coeff1,ldm_coeff2])
    
    def applyModel(self,LCs,times,model,transitTimes,duration,quarter,plot=False):
        """
        -take in a lightcurves and times for each transit
        -take in planet period
        
        -find places in lightcurve that doesn't have transits or with in 5 
        durations of a transit and apply a transit
        
        
        PAPER:
            residuals = observed data - expected
            don't use ingress/egress ???
            seperated into in in-residuals and out-residuals
            
        """
#        N = len(transitTimes)
#        if N != len(times) and N != len(LCs):
#            print("length error")
#            # throw error
#            pass
        
        for i in range(len(transitTimes)):
            t = transitTimes[i]
            lc = LCs[i]
            tt = times[i]
            start = (np.abs(np.array(tt)+(duration/48))).argmin()
            end = (np.abs(np.array(tt)-(duration/48))).argmin()
            residual = lc - model
            residualOUT = np.array(residual[:start])
            residualOUT = np.append(residualOUT,residual[end+1:])
            timeOUT = tt[:start]
            timeOUT = np.append(timeOUT,tt[end+1:])
            if plot:
#                plt.figure()
#                plt.plot(tt,residual,'r.')
#                plt.plot([tt[start],tt[start]],[-.002,.003],'g-')
#                plt.plot([tt[end],tt[end]],[-.002,.003],'k-')
#                plt.title("Residual of Quarter {} at Time {}".format(quarter,t))
                
                plt.figure()
                plt.title("In-transit vs Out-transit of Quarter " + str(quarter))
#                plt.plot(tt[:start],residual[:start],'b.')
#                plt.plot(tt[end+1:],residual[end+1:],'b.')
                plt.plot(timeOUT,residualOUT,'b.')
                plt.plot(tt[start:end+1],residual[start:end+1],'g.')
                plt.show()
            return residual[start:end+1],residualOUT,tt[start:end+1],timeOUT


