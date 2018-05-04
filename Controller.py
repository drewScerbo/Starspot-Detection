#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 12:27:13 2018

@author: as0216
"""
from Model import Model
#from View import View
import matplotlib.pyplot as plt
import pandas as pd
import os

"""
Q's:
    normalize the fakes?
    number of fakes?
    look at getData()
"""
    
model = Model()
#    view = View()
kois = model.getData() # get kois
print("number of KOIs to run: {}".format(len(kois)))

cols =  ['KOI_name','KS','avg_fake_KS','p','avg_fake_p','num_transits','star_spot']
df = pd.DataFrame(columns=cols)
filename = 'koi_tests.csv'

realKS,realP = 0,0
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
    print("number of transits: {}".format(len(transits)))
#    print(len(Time))
    IN,OUT=[],[]
    for i in range(len(transits)):
        diffLC,t = model.normalize(Time,Flux,Error,transits[i],dur, \
                           koi.koi_period,c,False)
        if len(diffLC) < 1:
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
    realKS, realP = model.getKS(IN,OUT,False)
    print("real ks: {:.4}, p: {:.4}".format(realKS,realP))
    
    avgKS,avgP,c = 0,0,0
    num_fake_trials = 20
#    plt.figure()
#    plt.plot(Time,Flux,'g.')
#    for j in range(len(transits)):
#        plt.plot([transits[j],transits[j]],[min(Flux),max(Flux)],'r-')
    KS,P = 0,0
    for j in range(num_fake_trials):
        fakeLCs,fakeTs = model.make_fakes(Flux,Time,m,transits,dur)
        c,avgKS,avgP = 0,0,0
        for k in range(len(fakeLCs)):
            fakeLCs[k],times = model.normalize(Time,Flux,Error,fakeTs[k],dur, \
                           koi.koi_period,c,False)
            if len(fakeLCs[k]) < 1:
                continue

            INfakes,OUTfakes,INt,OUTt = model.getINOUT(fakeLCs[k],times,dur,False)
            ks, p = model.getKS(INfakes,OUTfakes,False)
            if str(ks) != 'nan':
                avgKS += ks
                avgP += p
                c += 1
#        print("Running fake transit #{}".format(j))
        KS += avgKS/c
        P += avgP/c
#    plt.show()
    KS /= num_fake_trials
    P /= num_fake_trials
    print("fake ks: {:.4}, p: {:.4}".format(KS,P))

    print("diff: {}".format(abs(realKS - KS)))
    # save to pandas
    # 'KOI_name','KS','avg_fake_KS','p','avg_fake_p','star_spot'
    if abs(realKS - KS) >= P:
        print("KOI {} => YES".format(koi.kepoi_name))
        df.loc[len(df.index)] = [koi.kepoi_name,realKS,KS,realP,P,len(transits),'YES']
    else:
        print("KOI {} => NO".format(koi.kepoi_name))
        df.loc[len(df.index)] = [koi.kepoi_name,realKS,KS,realP,P,len(transits),'NO']
#    os.remove(filename)
    df.to_csv(filename,index=False)
    print('wrote to ' + filename)
    print()