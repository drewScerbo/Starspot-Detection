#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 12:27:13 2018

@author: as0216
"""
from Model import Model
#from View import View
#import matplotlib.pyplot as plt
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
filename = 'koi_tests.csv'
if os.path.exists(filename):
    df = pd.read_csv(filename)
    print("loaded in {}".format(filename))
else:
    print("created {}".format(filename))
    cols =  ['KOI_name','KS','avg_fake_KS','p','avg_fake_p','num_transits','star_spot']
    df = pd.DataFrame(columns=cols)

realKS,realP = 0,0
avgKS,avgP = 0.0,0.0
count = 0
for koi in kois:
    count += 1
    error = False
    print()
    print("Looking at {} of {}: {}".format(count,len(kois),koi.kepoi_name))
    if koi.kepoi_name in df.values:
        print("already calculated")
        continue
    c = 0
    Flux,Time = [],[]
    
    print("number of quarters: {}".format(len(\
          koi.get_light_curves(short_cadence=False))))
    dur = koi.koi_duration
        
        
    if koi.koi_num_transits == None:
        print("num_transits = None")
        error = True
#        continue
    if not error:
        Time,Flux,Error,transits = model.findTransits(koi,koi.koi_num_transits,\
                                        koi.koi_time0bk,koi.koi_period,0)
    print("number of transits: {}".format(len(transits)))
#    print(len(Time))
    IN,OUT=[],[]
    for i in range(len(transits)):
        if error:
            break
        diffLC,t = model.normalize(Time,Flux,Error,transits[i],dur, \
                           koi.koi_period,c,False)
        if len(diffLC) < 1:
            continue
        
        # tt,period,incl,hjd,dor,ror,ldm_coeff1,ldm_coeff2
        m = model.makeModel(transits[i],koi.koi_period,koi.koi_incl,t,\
                            koi.koi_dor,koi.koi_ror,koi.koi_ldm_coeff1,\
                            koi.koi_ldm_coeff2)

#            m = model.convolve(m)

        
#       lc,times,model,transitTimes,duration,quarter
        INres,OUTres,INt,OUTt = model.applyModel(diffLC,t,m,\
                                        dur,False)
        IN.extend(INres)
        OUT.extend(OUTres)
    if len(IN) < 1:
        print("len(IN) = 0")
        error = True
#        continue
    if not error:
        realKS, realP = model.getKS(IN,OUT,False)
        print("real ks: {:.4}, p: {:.4}".format(realKS,realP))
    
    avgKS,avgP,c = 0,0,0
    if error:
        num_fake_trials = 0
    else:
        num_fake_trials = 20
    KS,P = 0,0
    for j in range(num_fake_trials):
        if error:
            break
        if j == num_fake_trials/2:
            print('halfway through fake trials')
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
            else:
                k -= 1
        if c == 0:
            j -= 1
            continue
        KS += avgKS/c
        P += avgP/c
    if not error:
        KS /= num_fake_trials
        P /= num_fake_trials
        print("fake ks: {:.4}, p: {:.4}".format(KS,P))
    
        print("diff: {}".format(abs(realKS - KS)))
        # save to pandas
        # 'KOI_name','KS','avg_fake_KS','p','avg_fake_p','num_transits','star_spot'
        if abs(realKS - KS) >= P + realP:
            print("KOI {} => YES".format(koi.kepoi_name))
            df.loc[len(df.index)] = [koi.kepoi_name,realKS,KS,realP,P,len(transits),'YES']
        else:
            print("KOI {} => NO".format(koi.kepoi_name))
            df.loc[len(df.index)] = [koi.kepoi_name,realKS,KS,realP,P,len(transits),'NO']
#    os.remove(filename)
    else:
        df.loc[len(df.index)] = [koi.kepoi_name,-1,-1,-1,-1,-1,'ERROR']
    df.to_csv(filename,index=False)
    print('wrote to ' + filename)
