#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 12:27:13 2018

@author: as0216
"""
from Model import Model
from View import View
import matplotlib.pyplot as plt

class Controller:
    pass
    

if __name__ == '__main__':
    cont = Controller()
    model = Model()
    view = View()
    kois = model.getData() # get kois
    
    # for each koi in data, get LC and normalize the first LC
    for koi in kois:
        c = -1
#        Flux,Time = [],[]
        
#        print(len(koi.get_light_curves(short_cadence=False)))
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
            
            diffLC,t = model.normalize(lc,koi.koi_time0bk,koi.koi_duration, \
                                   koi.koi_period,koi.koi_num_transits,True)
#            plt.figure()
#            plt.plot(time,flux,'b.')
#            plt.show()
            
    
    
#    z = linspace(0, 1.2, 100)
#    gammavals = [[0., 0.], [1., 0.], [2., -1.]]
#    figure()
#    for gammas in gammavals:
#        f = occultquad(z, 0.1, gammas)
#        plot(z, f)