#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
import numpy as np
from pytc2 import filtros_digitales as test_module



fpass = 0.4
fstop = 0.6
    
fs = 2.0
ripple = 0.5  # dB
attenuation = 40  # dB
lgrid = 16

@pytest.mark.parametrize(
    "func_ptr, order, band_edges, desired, par_dict, true_val",
    [
       # FIR tipo 1 LP FIRPM
       ( test_module.fir_design_pm,
         22,
         [0, 0.4, 0.6, 1.0],
         [1, 1, 0, 0],
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'lowpass' },
         np.array([-7.19767604e-03, -9.18676455e-17,  1.36083987e-02,  3.73850272e-17,
                   -2.62924775e-02, -1.25558142e-16,  4.86004002e-02,  8.08736877e-17,
                   -9.64741264e-02, -1.81134067e-16,  3.14996419e-01,  5.00000000e-01,
                    3.14996419e-01, -1.81134067e-16, -9.64741264e-02,  8.08736877e-17,
                    4.86004002e-02, -1.25558142e-16, -2.62924775e-02,  3.73850272e-17,
                    1.36083987e-02, -9.18676455e-17, -7.19767604e-03]),         
       ),
       # FIR tipo 1 LP FIRLS
       ( test_module.fir_design_ls,
         22,
         [0, 0.4, 0.6, 1.0],
         [1, 1, 0, 0],
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'lowpass' },
         np.array([ -4.69791235e-03, -1.05904868e-15,  1.19150454e-02,  1.65145675e-15,
                    -2.45363416e-02, -2.49800181e-15,  4.70951653e-02,  3.08086889e-15,
                    -9.54323190e-02, -3.46944695e-15,  3.14641342e-01,  5.00000000e-01,
                     3.14641342e-01, -3.46944695e-15, -9.54323190e-02,  3.08086889e-15,
                     4.70951653e-02, -2.49800181e-15, -2.45363416e-02,  1.65145675e-15,
                     1.19150454e-02, -1.05904868e-15, -4.69791235e-03]),         
       ),
       # FIR tipo 1 BP FIRPM
       ( test_module.fir_design_pm,
         70,
         [0, .2, .3, .5, .6, 1.0],
         [0, 0, 1, 1, 0, 0],
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'multiband' },
         np.array([ -6.02677592e-04, -4.91931801e-05, -7.91129330e-05, -6.71427916e-04,
                     5.70959923e-04,  2.58909882e-03,  9.33335202e-04, -1.93598260e-03,
                    -7.60746961e-04, -4.15004609e-04, -4.36641271e-03, -2.21145698e-03,
                     7.20910028e-03,  7.02173246e-03, -1.76766911e-03, -3.99628646e-04,
                     2.14272655e-03, -1.16929931e-02, -1.65895581e-02,  7.02267493e-03,
                     1.95269372e-02,  3.13295510e-03,  3.97275511e-03,  1.86108426e-02,
                    -1.24296119e-02, -5.23266070e-02, -1.67401424e-02,  3.38624210e-02,
                     1.10469411e-02,  9.24141343e-03,  8.52817347e-02,  4.53310789e-02,
                    -1.66705930e-01, -2.07101144e-01,  8.93428693e-02,  3.00768456e-01,
                     8.93428693e-02, -2.07101144e-01, -1.66705930e-01,  4.53310789e-02,
                     8.52817347e-02,  9.24141343e-03,  1.10469411e-02,  3.38624210e-02,
                    -1.67401424e-02, -5.23266070e-02, -1.24296119e-02,  1.86108426e-02,
                     3.97275511e-03,  3.13295510e-03,  1.95269372e-02,  7.02267493e-03,
                    -1.65895581e-02, -1.16929931e-02,  2.14272655e-03, -3.99628646e-04,
                    -1.76766911e-03,  7.02173246e-03,  7.20910028e-03, -2.21145698e-03,
                    -4.36641271e-03, -4.15004609e-04, -7.60746961e-04, -1.93598260e-03,
                     9.33335202e-04,  2.58909882e-03,  5.70959923e-04, -6.71427916e-04,
                    -7.91129330e-05, -4.91931801e-05, -6.02677592e-04]),         
       ),
       # FIR tipo 1 BP FIRLS
       ( test_module.fir_design_ls,
         70,
         [0, .2, .3, .5, .6, 1.0],
         [0, 0, 1, 1, 0, 0],
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'multiband' },
         np.array([ -3.33382042e-04, -1.38105643e-04, -6.69749731e-05, -4.90105573e-04,
                     3.93135051e-04,  2.04355186e-03,  8.48504401e-04, -1.56546889e-03,
                    -7.23864165e-04, -4.50581191e-04, -3.91667785e-03, -2.11689430e-03,
                     6.42636415e-03,  6.55154816e-03, -1.45437770e-03, -3.57365400e-04,
                     2.00411785e-03, -1.10641449e-02, -1.60234658e-02,  6.42693297e-03,
                     1.88215360e-02,  3.22363166e-03,  3.93995227e-03,  1.83095710e-02,
                    -1.19657025e-02, -5.15946577e-02, -1.69381269e-02,  3.33195224e-02,
                     1.09967245e-02,  9.13388300e-03,  8.50153953e-02,  4.56186888e-02,
                    -1.66233940e-01, -2.07198881e-01,  8.90896872e-02,  3.00693177e-01,
                     8.90896872e-02, -2.07198881e-01, -1.66233940e-01,  4.56186888e-02,
                     8.50153953e-02,  9.13388300e-03,  1.09967245e-02,  3.33195224e-02,
                    -1.69381269e-02, -5.15946577e-02, -1.19657025e-02,  1.83095710e-02,
                     3.93995227e-03,  3.22363166e-03,  1.88215360e-02,  6.42693297e-03,
                    -1.60234658e-02, -1.10641449e-02,  2.00411785e-03, -3.57365400e-04,
                    -1.45437770e-03,  6.55154816e-03,  6.42636415e-03, -2.11689430e-03,
                    -3.91667785e-03, -4.50581191e-04, -7.23864165e-04, -1.56546889e-03,
                     8.48504401e-04,  2.04355186e-03,  3.93135051e-04, -4.90105573e-04,
                    -6.69749731e-05, -1.38105643e-04, -3.33382042e-04]),         
       ),
       # FIR tipo 2 BP FIRPM
       ( test_module.fir_design_pm,
         71,
         [0, .2, .3, .5, .6, 1.0],
         [0, 0, 1, 1, 0, 0],
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'multiband' },
         np.array([ -3.54969240e-04, -5.23347076e-04,  1.51332495e-04, -5.32773747e-04,
                    -3.45740496e-04,  1.81273295e-03,  2.28377249e-03, -9.20211592e-04,
                    -1.85823348e-03,  2.82972112e-05, -2.36614583e-03, -4.72686880e-03,
                     2.65115000e-03,  9.12201443e-03,  2.28858943e-03, -2.80903625e-03,
                     2.44671533e-03, -3.40441036e-03, -1.78854769e-02, -6.75033420e-03,
                     1.78420027e-02,  1.25777553e-02, -9.98777368e-04,  1.37399098e-02,
                     9.77576384e-03, -3.84339694e-02, -4.42504661e-02,  1.55821754e-02,
                     2.98146948e-02, -1.94182354e-03,  4.56945888e-02,  9.34272685e-02,
                    -5.49126096e-02, -2.32762347e-01, -8.48638288e-02,  2.41007249e-01,
                     2.41007249e-01, -8.48638288e-02, -2.32762347e-01, -5.49126096e-02,
                     9.34272685e-02,  4.56945888e-02, -1.94182354e-03,  2.98146948e-02,
                     1.55821754e-02, -4.42504661e-02, -3.84339694e-02,  9.77576384e-03,
                     1.37399098e-02, -9.98777368e-04,  1.25777553e-02,  1.78420027e-02,
                    -6.75033420e-03, -1.78854769e-02, -3.40441036e-03,  2.44671533e-03,
                    -2.80903625e-03,  2.28858943e-03,  9.12201443e-03,  2.65115000e-03,
                    -4.72686880e-03, -2.36614583e-03,  2.82972112e-05, -1.85823348e-03,
                    -9.20211592e-04,  2.28377249e-03,  1.81273295e-03, -3.45740496e-04,
                    -5.32773747e-04,  1.51332495e-04, -5.23347076e-04, -3.54969240e-04]),         
       ),
       # FIR tipo 2 BP FIRLS
       ( test_module.fir_design_ls,
         71,
         [0, .2, .3, .5, .6, 1.0],
         [0, 0, 1, 1, 0, 0],
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'multiband' },
         np.array([ -2.75824911e-04, -2.81438633e-04,  6.28022974e-05, -2.98459216e-04,
                    -2.39471918e-04,  1.50084486e-03,  1.83236691e-03, -8.45390014e-04,
                    -1.62156045e-03,  1.26013513e-05, -2.12466977e-03, -4.17311052e-03,
                     2.56262670e-03,  8.40933561e-03,  2.02051724e-03, -2.69724680e-03,
                     2.23434799e-03, -3.36393798e-03, -1.70656660e-02, -6.23598296e-03,
                     1.74073579e-02,  1.21698414e-02, -1.02047998e-03,  1.33669329e-02,
                     9.33521565e-03, -3.80475141e-02, -4.35239943e-02,  1.57254198e-02,
                     2.96881632e-02, -1.86754999e-03,  4.55638701e-02,  9.29261077e-02,
                    -5.52185769e-02, -2.32616083e-01, -8.46027952e-02,  2.41148499e-01,
                     2.41148499e-01, -8.46027952e-02, -2.32616083e-01, -5.52185769e-02,
                     9.29261077e-02,  4.55638701e-02, -1.86754999e-03,  2.96881632e-02,
                     1.57254198e-02, -4.35239943e-02, -3.80475141e-02,  9.33521565e-03,
                     1.33669329e-02, -1.02047998e-03,  1.21698414e-02,  1.74073579e-02,
                    -6.23598296e-03, -1.70656660e-02, -3.36393798e-03,  2.23434799e-03,
                    -2.69724680e-03,  2.02051724e-03,  8.40933561e-03,  2.56262670e-03,
                    -4.17311052e-03, -2.12466977e-03,  1.26013513e-05, -1.62156045e-03,
                    -8.45390014e-04,  1.83236691e-03,  1.50084486e-03, -2.39471918e-04,
                    -2.98459216e-04,  6.28022974e-05, -2.81438633e-04, -2.75824911e-04]),         
       ),
       # FIR tipo 1 BE FIRPM
       ( test_module.fir_design_pm,
         70,
         [0, .2, .3, .5, .6, 1.0],
         [1, 1, 0, 0, 1, 1],
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'multiband' },
         np.array([  6.02677592e-04,  4.91931801e-05,  7.91129330e-05,  6.71427916e-04,
                    -5.70959923e-04, -2.58909882e-03, -9.33335202e-04,  1.93598260e-03,
                     7.60746961e-04,  4.15004609e-04,  4.36641271e-03,  2.21145698e-03,
                    -7.20910028e-03, -7.02173246e-03,  1.76766911e-03,  3.99628646e-04,
                    -2.14272655e-03,  1.16929931e-02,  1.65895581e-02, -7.02267493e-03,
                    -1.95269372e-02, -3.13295510e-03, -3.97275511e-03, -1.86108426e-02,
                     1.24296119e-02,  5.23266070e-02,  1.67401424e-02, -3.38624210e-02,
                    -1.10469411e-02, -9.24141343e-03, -8.52817347e-02, -4.53310789e-02,
                     1.66705930e-01,  2.07101144e-01, -8.93428693e-02,  6.99231544e-01,
                    -8.93428693e-02,  2.07101144e-01,  1.66705930e-01, -4.53310789e-02,
                    -8.52817347e-02, -9.24141343e-03, -1.10469411e-02, -3.38624210e-02,
                     1.67401424e-02,  5.23266070e-02,  1.24296119e-02, -1.86108426e-02,
                    -3.97275511e-03, -3.13295510e-03, -1.95269372e-02, -7.02267493e-03,
                     1.65895581e-02,  1.16929931e-02, -2.14272655e-03,  3.99628646e-04,
                     1.76766911e-03, -7.02173246e-03, -7.20910028e-03,  2.21145698e-03,
                     4.36641271e-03,  4.15004609e-04,  7.60746961e-04,  1.93598260e-03,
                    -9.33335202e-04, -2.58909882e-03, -5.70959923e-04,  6.71427916e-04,
                     7.91129330e-05,  4.91931801e-05,  6.02677592e-04]),         
       ),
       # FIR tipo 1 BE FIRLS
       ( test_module.fir_design_ls,
         70,
         [0, .2, .3, .5, .6, 1.0],
         [1, 1, 0, 0, 1, 1],
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'multiband' },
         np.array([  3.33382042e-04,  1.38105643e-04,  6.69749730e-05,  4.90105573e-04,
                    -3.93135051e-04, -2.04355186e-03, -8.48504401e-04,  1.56546889e-03,
                     7.23864165e-04,  4.50581191e-04,  3.91667785e-03,  2.11689430e-03,
                    -6.42636415e-03, -6.55154816e-03,  1.45437770e-03,  3.57365400e-04,
                    -2.00411785e-03,  1.10641449e-02,  1.60234658e-02, -6.42693297e-03,
                    -1.88215360e-02, -3.22363166e-03, -3.93995227e-03, -1.83095710e-02,
                     1.19657025e-02,  5.15946577e-02,  1.69381269e-02, -3.33195224e-02,
                    -1.09967245e-02, -9.13388300e-03, -8.50153953e-02, -4.56186888e-02,
                     1.66233940e-01,  2.07198881e-01, -8.90896872e-02,  6.99306823e-01,
                    -8.90896872e-02,  2.07198881e-01,  1.66233940e-01, -4.56186888e-02,
                    -8.50153953e-02, -9.13388300e-03, -1.09967245e-02, -3.33195224e-02,
                     1.69381269e-02,  5.15946577e-02,  1.19657025e-02, -1.83095710e-02,
                    -3.93995227e-03, -3.22363166e-03, -1.88215360e-02, -6.42693297e-03,
                     1.60234658e-02,  1.10641449e-02, -2.00411785e-03,  3.57365400e-04,
                     1.45437770e-03, -6.55154816e-03, -6.42636415e-03,  2.11689430e-03,
                     3.91667785e-03,  4.50581191e-04,  7.23864165e-04,  1.56546889e-03,
                    -8.48504401e-04, -2.04355186e-03, -3.93135051e-04,  4.90105573e-04,
                     6.69749730e-05,  1.38105643e-04,  3.33382042e-04]),         
       ),       
       # FIR tipo 3 derivador FIRPM
       ( test_module.fir_design_pm,
         70,
         [0, 0.9],
         [0, 1/np.pi] ,
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'd' },
         np.array([ -3.84615856e-06,  8.09130415e-06, -1.63356679e-05,  2.98237157e-05,
                    -5.07920312e-05,  8.21071421e-05, -1.27365025e-04,  1.90993138e-04,
                    -2.78351818e-04,  3.95833353e-04, -5.50958089e-04,  7.52468612e-04,
                    -1.01042244e-03,  1.33629138e-03, -1.74306781e-03,  2.24539981e-03,
                    -2.85975701e-03,  3.60468694e-03, -4.50114673e-03,  5.57304368e-03,
                    -6.84798663e-03,  8.35844876e-03, -1.01435098e-02,  1.22514963e-02,
                    -1.47440825e-02,  1.77027429e-02, -2.12392043e-02,  2.55130073e-02,
                    -3.07623614e-02,  3.73615367e-02, -4.59357232e-02,  5.76137357e-02,
                    -7.46597080e-02,  1.02367963e-01, -1.56642036e-01,  3.17046566e-01,
                     0.00000000e+00, -3.17046566e-01,  1.56642036e-01, -1.02367963e-01,
                     7.46597080e-02, -5.76137357e-02,  4.59357232e-02, -3.73615367e-02,
                     3.07623614e-02, -2.55130073e-02,  2.12392043e-02, -1.77027429e-02,
                     1.47440825e-02, -1.22514963e-02,  1.01435098e-02, -8.35844876e-03,
                     6.84798663e-03, -5.57304368e-03,  4.50114673e-03, -3.60468694e-03,
                     2.85975701e-03, -2.24539981e-03,  1.74306781e-03, -1.33629138e-03,
                     1.01042244e-03, -7.52468612e-04,  5.50958089e-04, -3.95833353e-04,
                     2.78351818e-04, -1.90993138e-04,  1.27365025e-04, -8.21071421e-05,
                     5.07920312e-05, -2.98237157e-05,  1.63356679e-05, -8.09130415e-06,
                     3.84615856e-06]),         
       ),
       # FIR tipo 3 derivador FIRLS
       ( test_module.fir_design_ls,
         70,
         [0, 0.9],
         [0, 1/np.pi] ,
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'd' },
         np.array([ -1.76027317e-06,  5.47464898e-06, -1.16425621e-05,  2.24752396e-05,
                    -3.97950820e-05,  6.63944944e-05, -1.05651972e-04,  1.61891342e-04,
                    -2.40325604e-04,  3.47287994e-04, -4.90245785e-04,  6.77977295e-04,
                    -9.20613589e-04,  1.22979141e-03, -1.61872359e-03,  2.10236744e-03,
                    -2.69756519e-03,  3.42330731e-03, -4.30105056e-03,  5.35524325e-03,
                    -6.61406163e-03,  8.11056129e-03, -9.88438245e-03,  1.19843902e-02,
                    -1.44727342e-02,  1.74312907e-02, -2.09720917e-02,  2.52548763e-02,
                    -3.05179217e-02,  3.71354364e-02, -4.57324104e-02,  5.74373260e-02,
                    -7.45138540e-02,  1.02255742e-01, -1.56565846e-01,  3.17008053e-01,
                     0.00000000e+00, -3.17008053e-01,  1.56565846e-01, -1.02255742e-01,
                     7.45138540e-02, -5.74373260e-02,  4.57324104e-02, -3.71354364e-02,
                     3.05179217e-02, -2.52548763e-02,  2.09720917e-02, -1.74312907e-02,
                     1.44727342e-02, -1.19843902e-02,  9.88438245e-03, -8.11056129e-03,
                     6.61406163e-03, -5.35524325e-03,  4.30105056e-03, -3.42330731e-03,
                     2.69756519e-03, -2.10236744e-03,  1.61872359e-03, -1.22979141e-03,
                     9.20613589e-04, -6.77977295e-04,  4.90245785e-04, -3.47287994e-04,
                     2.40325604e-04, -1.61891342e-04,  1.05651972e-04, -6.63944944e-05,
                     3.97950820e-05, -2.24752396e-05,  1.16425621e-05, -5.47464898e-06,
                     1.76027317e-06]),         
       ),       
       # FIR tipo 3 Hilbert FIRPM
       ( test_module.fir_design_pm,
         70,
         [0.05, 0.95],
         [1., 1.] ,
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'h' },
         np.array([  3.87456337e-16,  1.17533958e-03, -1.07441083e-16,  1.44354524e-03,
                     2.26628089e-16,  2.29271717e-03,  1.84231011e-16,  3.44140504e-03,
                     2.02831885e-16,  4.95299346e-03,  2.81517144e-16,  6.90112583e-03,
                     1.81243184e-16,  9.37973865e-03,  2.78796270e-16,  1.24980059e-02,
                     4.36566263e-16,  1.64039852e-02,  2.49270957e-16,  2.12998064e-02,
                     3.20849046e-16,  2.74881647e-02,  1.81865831e-16,  3.54480325e-02,
                     1.77209723e-16,  4.60053131e-02,  4.17230915e-16,  6.07215611e-02,
                     1.17334395e-16,  8.29625121e-02,  4.06375932e-16,  1.21514507e-01,
                    -5.13944729e-16,  2.08677424e-01,  3.78736171e-16,  6.35435950e-01,
                     0.00000000e+00, -6.35435950e-01, -3.78736171e-16, -2.08677424e-01,
                     5.13944729e-16, -1.21514507e-01, -4.06375932e-16, -8.29625121e-02,
                    -1.17334395e-16, -6.07215611e-02, -4.17230915e-16, -4.60053131e-02,
                    -1.77209723e-16, -3.54480325e-02, -1.81865831e-16, -2.74881647e-02,
                    -3.20849046e-16, -2.12998064e-02, -2.49270957e-16, -1.64039852e-02,
                    -4.36566263e-16, -1.24980059e-02, -2.78796270e-16, -9.37973865e-03,
                    -1.81243184e-16, -6.90112583e-03, -2.81517144e-16, -4.95299346e-03,
                    -2.02831885e-16, -3.44140504e-03, -1.84231011e-16, -2.29271717e-03,
                    -2.26628089e-16, -1.44354524e-03,  1.07441083e-16, -1.17533958e-03,
                    -3.87456337e-16]),         
       ),
       # FIR tipo 3 Hilbert FIRLS
       ( test_module.fir_design_ls,
         70,
         [0.05, 0.95],
         [1., 1.] ,
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'h' },
         np.array([ -7.07104564e-15,  8.27294033e-04, -3.74939926e-15,  1.23702988e-03,
                    -3.25451657e-15,  2.02063793e-03, -3.01636692e-15,  3.10895600e-03,
                    -3.12399514e-15,  4.56364015e-03, -3.35237992e-15,  6.46183581e-03,
                    -3.87437609e-15,  8.89624359e-03, -4.36550898e-15,  1.19815924e-02,
                    -4.97510224e-15,  1.58669070e-02, -5.73747964e-15,  2.07576329e-02,
                    -6.17124586e-15,  2.69565803e-02, -6.38360740e-15,  3.49437253e-02,
                    -6.56050964e-15,  4.55439944e-02, -6.27347362e-15,  6.03194193e-02,
                    -5.52744753e-15,  8.26337232e-02, -4.21835712e-15,  1.21271154e-01,
                    -3.49809584e-15,  2.08527791e-01, -1.33666332e-15,  6.35385567e-01,
                     0.00000000e+00, -6.35385567e-01,  1.33666332e-15, -2.08527791e-01,
                     3.49809584e-15, -1.21271154e-01,  4.21835712e-15, -8.26337232e-02,
                     5.52744753e-15, -6.03194193e-02,  6.27347362e-15, -4.55439944e-02,
                     6.56050964e-15, -3.49437253e-02,  6.38360740e-15, -2.69565803e-02,
                     6.17124586e-15, -2.07576329e-02,  5.73747964e-15, -1.58669070e-02,
                     4.97510224e-15, -1.19815924e-02,  4.36550898e-15, -8.89624359e-03,
                     3.87437609e-15, -6.46183581e-03,  3.35237992e-15, -4.56364015e-03,
                     3.12399514e-15, -3.10895600e-03,  3.01636692e-15, -2.02063793e-03,
                     3.25451657e-15, -1.23702988e-03,  3.74939926e-15, -8.27294033e-04,
                     7.07104564e-15]),         
       ),       
       # FIR tipo 4 Hilbert FIRPM
       ( test_module.fir_design_pm,
         71,
         [0.05, 1.],
         [1., 1.] ,
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'h' },
         np.array([ -8.14975386e-04, -5.06654763e-04, -6.61314337e-04, -8.44461309e-04,
                    -1.05908535e-03, -1.30882466e-03, -1.59795740e-03, -1.93047813e-03,
                    -2.31004728e-03, -2.74136080e-03, -3.23037275e-03, -3.78225870e-03,
                    -4.40149413e-03, -5.09819226e-03, -5.87745332e-03, -6.74977922e-03,
                    -7.72515099e-03, -8.81632127e-03, -1.00380608e-02, -1.14086630e-02,
                    -1.29513086e-02, -1.46950349e-02, -1.66772145e-02, -1.89475124e-02,
                    -2.15732109e-02, -2.46459875e-02, -2.82982736e-02, -3.27238758e-02,
                    -3.82214804e-02, -4.52740176e-02, -5.47175045e-02, -6.81325781e-02,
                    -8.89085702e-02, -1.25862077e-01, -2.11326724e-01, -6.36326109e-01,
                     6.36326109e-01,  2.11326724e-01,  1.25862077e-01,  8.89085702e-02,
                     6.81325781e-02,  5.47175045e-02,  4.52740176e-02,  3.82214804e-02,
                     3.27238758e-02,  2.82982736e-02,  2.46459875e-02,  2.15732109e-02,
                     1.89475124e-02,  1.66772145e-02,  1.46950349e-02,  1.29513086e-02,
                     1.14086630e-02,  1.00380608e-02,  8.81632127e-03,  7.72515099e-03,
                     6.74977922e-03,  5.87745332e-03,  5.09819226e-03,  4.40149413e-03,
                     3.78225870e-03,  3.23037275e-03,  2.74136080e-03,  2.31004728e-03,
                     1.93047813e-03,  1.59795740e-03,  1.30882466e-03,  1.05908535e-03,
                     8.44461309e-04,  6.61314337e-04,  5.06654763e-04,  8.14975386e-04]),         
       ),
       # FIR tipo 4 Hilbert FIRLS
       ( test_module.fir_design_ls,
         71,
         [0.05, 1.],
         [1., 1.] ,
         { 'grid_density': lgrid, 'fs':fs, 'filter_type': 'h' },
         np.array([ -4.88815650e-04, -4.32269141e-04, -5.56970712e-04, -7.21895081e-04,
                    -9.18900771e-04, -1.15394899e-03, -1.42652726e-03, -1.74452382e-03,
                    -2.10842745e-03, -2.52694482e-03, -3.00139968e-03, -3.54143799e-03,
                    -4.14949191e-03, -4.83655261e-03, -5.60671687e-03, -6.47306731e-03,
                    -7.44236247e-03, -8.53111244e-03, -9.75053974e-03, -1.11230367e-02,
                    -1.26676746e-02, -1.44174653e-02, -1.64060781e-02, -1.86869630e-02,
                    -2.13230235e-02, -2.44105983e-02, -2.80774405e-02, -3.25215780e-02,
                    -3.80374754e-02, -4.51120429e-02, -5.45770377e-02, -6.80167669e-02,
                    -8.88166021e-02, -1.25796680e-01, -2.11286494e-01, -6.36313232e-01,
                     6.36313232e-01,  2.11286494e-01,  1.25796680e-01,  8.88166021e-02,
                     6.80167669e-02,  5.45770377e-02,  4.51120429e-02,  3.80374754e-02,
                     3.25215780e-02,  2.80774405e-02,  2.44105983e-02,  2.13230235e-02,
                     1.86869630e-02,  1.64060781e-02,  1.44174653e-02,  1.26676746e-02,
                     1.11230367e-02,  9.75053974e-03,  8.53111244e-03,  7.44236247e-03,
                     6.47306731e-03,  5.60671687e-03,  4.83655261e-03,  4.14949191e-03,
                     3.54143799e-03,  3.00139968e-03,  2.52694482e-03,  2.10842745e-03,
                     1.74452382e-03,  1.42652726e-03,  1.15394899e-03,  9.18900771e-04,
                     7.21895081e-04,  5.56970712e-04,  4.32269141e-04,  4.88815650e-04]),         
       ),       
       
    ]
)
def test_parametrize_fir_design_pm_ls( func_ptr, order, band_edges, desired, par_dict, true_val ):

    
    b_coeffs = func_ptr(order, band_edges, desired, **par_dict)
    
    # firPM devuelve 3 resultados
    if len(b_coeffs) == 3:
        b_coeffs = b_coeffs[0]
    
    assert  np.allclose(true_val, b_coeffs)
        

@pytest.mark.parametrize(
    "func_ptr",
    [
       # FIR tipo 1 LP FIRPM
       test_module.fir_design_pm,
       test_module.fir_design_ls,
    ]
)
def test_parametrize_fir_design_invalid( func_ptr ):

    
    fs = 2.0
    ftype = 'multiband'
    lgrid = 16 # Aumentamos el orden para mejor selectividad

    band_edges = [0, 0.4, 0.6, 1.0]
    desired = [0, 0, 1, 1, 0, 0]

    order = '71 ' # Aumentamos el orden para mejor selectividad

    par_dict = { 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)
    
    order = -71  # Aumentamos el orden para mejor selectividad

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    order = 71  # Aumentamos el orden para mejor selectividad
   
    lgrid = '16 ' # Aumentamos el orden para mejor selectividad
    par_dict = { 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)
    
    lgrid = -16  # Aumentamos el orden para mejor selectividad
    par_dict = { 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    lgrid = 16  # Aumentamos el orden para mejor selectividad
    fs = '2'  # Aumentamos el orden para mejor selectividad
    par_dict = { 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    fs = -2  # Aumentamos el orden para mejor selectividad
    par_dict = { 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    fs = 2.  # Aumentamos el orden para mejor selectividad
    par_dict = { 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    ftype = '2'
    par_dict = { 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    ftype = 2
    par_dict = { 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    ftype = 'multiband'
    band_edges = 0

    par_dict = { 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)
        
    band_edges = [0]

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)
        
    band_edges = [0, 0.1, 0.3]

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    band_edges = [0, 0.4, 0.6, 1.0]

    desired = 0

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    desired = [0]

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)
        
    desired = [0, 1.]

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

        
    desired = [0, 0, 1, 1, 0, 0]
    par_dict = { 'weight' : 1. , 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    band_edges = [0, 0.4, 0.6, 1.0]

    par_dict = { 'weight' : [1., 1., 1.] , 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

    band_edges = [0, 0.4, 0.6, 1.0]

    ftype = 'h'
    
    par_dict = { 'weight' : [1., 1., 1.] , 'grid_density': lgrid, 'fs':fs, 'filter_type': ftype }

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(order, band_edges, desired, **par_dict)

