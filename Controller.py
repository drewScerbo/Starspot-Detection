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
choice = ""
while not choice.endswith('.csv') and choice != 'd':
    choice = input("Enter filename ending in '.csv' or 'd' for default: ")
if choice == 'd':
    filename = 'star_data.csv'
else:
    filename = choice
replace = False
if os.path.exists(filename):
    choice = input("Replace existing file? [y/n] ")
    if choice == 'n':
        df = pd.read_csv(filename)
        replace = input("Replace error files? [y/n] ")
        if replace == 'y':
            replace = True
        else:
            replace = False
        print("loaded in {}".format(filename))
    else:
        print("created {}".format(filename))
        cols =  ['KOI_name','KS','fake_KSs','p','fake_Ps','num_trials',\
                 'fakeAvg','fakeMed','fakeStd','fakeMad','realAvgRatio','realMedRatio']
        df = pd.DataFrame(columns=cols)
else:
    print("created {}".format(filename))
    cols =  ['KOI_name','KS','fake_KSs','p','fake_Ps','num_trials',\
             'fakeAvg','fakeMed','fakeStd','fakeMad','realAvgRatio','realMedRatio']
    df = pd.DataFrame(columns=cols)

while True:
    n = input("Enter integer for number of fake trial tests: ")
    if not isinstance(n,int):
        num_trials = int(n)
        break

print('Loading Kepler data...')
kois = model.getData() # get kois
print("number of KOIs to run: {}".format(len(kois)))
realKS,realP = 0,0
avgKS,avgP = 0.0,0.0
count = -1
for koi in kois:
    count += 1
    error = False
    print()
    print("Looking at {} of {}: {}".format(count+1,len(kois),koi.kepoi_name))
    if koi.kepoi_name in df.KOI_name.values:
        print("already calculated")
        if replace and -1 in df.loc[count].values:
            print("replacing error file")
        else:
            continue
            
    c = 0
    Flux,Time = [],[]
    
    print("number of quarters: {}".format(len(\
          koi.get_light_curves(short_cadence=False))))
    dur = koi.koi_duration
        
    if model.isNone(koi):
        error = True
        transits = []
    else:
        Time,Flux,Error,transits = model.findTransits(koi,koi.koi_num_transits,\
                                        koi.koi_time0bk,koi.koi_period)
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

#       lc,times,model,transitTimes,duration,quarter
        INres,OUTres,INt,OUTt = model.applyModel(diffLC,t,m,\
                                        dur,False)
        IN.extend(INres)
        OUT.extend(OUTres)
    if len(IN) < 1:
        error = True
    if not error:
        realKS, realP = model.getKS(IN,OUT,koi.kepoi_name,False)
        print("real ks: {:.4}, p: {:.4}".format(realKS,realP))
    
    avgKS,avgP,c = 0,0,0
    if error:
        num_fake_trials = 0
    else:
        num_fake_trials = num_trials
    KS,P = '',''
    for j in range(num_fake_trials):
        if error:
            break
        fakeLCs,fakeTs = model.make_fakes(Flux,Time,m,transits,dur)
        
        c,avgKS,avgP = 0,0,0
        IN, OUT = [],[]
        for k in range(len(fakeLCs)):
            fakeLCs[k],times = model.normalize(Time,Flux,Error,fakeTs[k],dur, \
                           koi.koi_period,c,False)
            if len(fakeLCs[k]) < 1:
                continue

            INfakes,OUTfakes,INt,OUTt = model.getINOUT(fakeLCs[k],times,dur,False)
            IN.extend(INfakes)
            OUT.extend(OUTfakes)
        ks, p = model.getKS(INfakes,OUTfakes,koi.kepoi_name,False)
        if str(ks) == 'nan':
            continue
        if j == 0:
            KS += str(ks)
            P += str(p)
        KS += ';' + str(ks)
        P += ';' + str(p)
    if not error and len(KS) > 0:
        print("fake ks: {:.4}, p: {:.4}".format(KS,P))
        
        values = model.getValues(realKS,KS)
        
        # save to file
        if replace and -1 in df.loc[count].values:
            df.loc[count] = [koi.kepoi_name,realKS,KS,realP,P,num_fake_trials,\
               values[0],values[1],values[2],values[3],values[4],values[5]]
        else:
            df.loc[len(df.index)] = [koi.kepoi_name,realKS,KS,realP,P,num_fake_trials,\
               values[0],values[1],values[2],values[3],values[4],values[5]]
        
    else:
        print('Error found')
        df.loc[len(df.index)] = [koi.kepoi_name,-1,'ERROR',-1,'ERROR',num_fake_trials,\
               -1,-1,-1,-1,-1,-1]
    df.to_csv(filename,index=False)
    print('wrote to ' + filename)

