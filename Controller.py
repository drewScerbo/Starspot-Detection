#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 12:27:13 2018

@author: as0216
"""
from Model import Model
#from View import View
import matplotlib.pyplot as plt

"""
Q's:
    normalize the fakes?
    number of fakes?
    look at getData()
"""
    
model = Model()
#    view = View()
kois = model.getData() # get kois

# for each koi in data, get LC and normalize the first LC
avgKS,avgP = 0.0,0.0
for koi in kois:
    c = 0
#        Flux,Time = [],[]
    
    print("number of quarters: {}".format(len(\
          koi.get_light_curves(short_cadence=False))))
    dur = koi.koi_duration
#    for lc in koi.get_light_curves(short_cadence=False):
#        # lc => data for 1 quarter
#        time, flux = [], []
#        with lc.open() as f:
#            # The lightcurve data are in the first FITS HDU.
#            hdu_data = f[1].data
#            time.append(hdu_data["time"])
#            flux.append(hdu_data["sap_flux"])
#        if c > 0:
#            break
#        c +=1
        
        
        
    Time,Flux,Error,transits = model.findTransits(koi,koi.koi_num_transits,\
                                        koi.koi_time0bk,koi.koi_period,0)
#    print(len(Time))
    IN,OUT=[],[]
    for i in range(len(transits)):
        diffLC,t = model.normalize(Time,Flux,Error,transits[i],dur, \
                           koi.koi_period,c,False)
        if len(diffLC) < 10:
            continue
        
        # tt,period,incl,hjd,dor,ror,ldm_coeff1,ldm_coeff2
        m = model.makeModel(transits[i],koi.koi_period,koi.koi_incl,t,\
                            koi.koi_dor,koi.koi_ror,koi.koi_ldm_coeff1,\
                            koi.koi_ldm_coeff2)

#            m = model.convolve(m)

        
#            lc,times,model,transitTimes,duration,quarter
        INres,OUTres,INt,OUTt = model.applyModel(diffLC,t,m,\
                                        dur,False)
        IN.extend(INres)
        OUT.extend(OUTres)
    ks, p = model.getKS(IN,OUT,False)
    print("real ks: {:.4}, p: {:.4}".format(ks,p))
    
    avgKS,avgP,c = 0,0,0
    num_fake_trials = 20
#    plt.figure()
#    plt.plot(Time,Flux,'g.')
#    for j in range(len(transits)):
#        plt.plot([transits[j],transits[j]],[min(Flux),max(Flux)],'r-')
    for j in range(num_fake_trials):
        fakeLCs,fakeTs = model.make_fakes(Flux,Time,m,transits,dur)
        c,avgKS,avgP = 0,0,0
        for k in range(len(fakeLCs)):
            fakeLCs[k],times = model.normalize(Time,Flux,Error,fakeTs[k],dur, \
                           koi.koi_period,c,False)
            if len(fakeLCs[k]) < 7:
                continue
#        print("lenLC: {}".format(len(fakeLCs[j])))
#        print("lenT: {}".format(len(list(times))))
#        print("lenM: {}".format(len(list(m))))
#        print()
#        plt.plot([fakeTs[j],fakeTs[j]],[min(Flux),max(Flux)],'b--')

            INfakes,OUTfakes,INt,OUTt = model.getINOUT(fakeLCs[k],times,dur,False)
            ks, p = model.getKS(INfakes,OUTfakes,False)
            if str(ks) != 'nan':
                avgKS += ks
                avgP += p
                c += 1
        print("({})fake ks: {:.4}, p: {:.4}".format(j,avgKS/c,avgP,c))
#    plt.show()
    print("fake ks: {:.4}, p: {:.4}".format(avgKS/c,avgP/c))
                
#print("Average KS: {:.4} +- {:.4}".format(avgKS/c,avgP/c))
                