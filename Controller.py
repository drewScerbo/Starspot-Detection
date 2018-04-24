#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 12:27:13 2018

@author: as0216
"""
from Model import Model
from View import View
import matplotlib.pyplot as plt
import numpy as np

class Controller:
    pass
    

if __name__ == '__main__':
    cont = Controller()
    model = Model()
    view = View()
    kois = model.getData() # get kois
    
    # for each koi in data, get LC and normalize the first LC
    for koi in kois:
        c = 0
#        Flux,Time = [],[]
        
        print("number of quarters: {}".format(len(\
              koi.get_light_curves(short_cadence=False))))
        for lc in koi.get_light_curves(short_cadence=False):
            # lc => data for 1 quarter
            time, flux = [], []
            with lc.open() as f:
                # The lightcurve data are in the first FITS HDU.
                hdu_data = f[1].data
                time.append(hdu_data["time"])
                flux.append(hdu_data["sap_flux"])
            if c > 0:
                break
            c +=1
            
            diffLC,t,transits = model.normalize(lc,koi.koi_time0bk,koi.koi_duration, \
                                   koi.koi_period,koi.koi_num_transits,c,True)
            for i in range(len(transits)):
                # tt,period,incl,hjd,dor,ror,ldm_coeff1,ldm_coeff2
                m = model.makeModel(transits[i],koi.koi_period,koi.koi_incl,t[i],\
                                    koi.koi_dor,koi.koi_ror,koi.koi_ldm_coeff1,\
                                    koi.koi_ldm_coeff2)
                plt.figure()
                plt.plot(t[i],diffLC[i],'b.')
                plt.plot(t[i],m,'g-')
                plt.show()
                INres,OUTres,INt,OUTt = model.applyModel(diffLC[i],t[i],m,transits,koi.koi_duration,c,True)
                ks, p = model.getKS(INres,OUTres)
                print("ks: {:.4}, p: {:.4}".format(ks,p))
                