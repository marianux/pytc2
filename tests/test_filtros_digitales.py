#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
import numpy as np
from pytc2 import filtros_digitales as test_module

    
# if Ftype == 'd':

#     # fine tuning
    
#     # tipo 3
#     N = 16
#     Be = [0, 0.9]
#     D = [0, np.pi]
#     W = [1., 1.]
#     Be_sig = [0, 0, 0.9, 1.0]
#     D_sig = [1., 1.]
    
#     # tipo 4
#     # N = 15
#     # Be = [0, 1.]
#     # D = [0, np.pi]
#     # Be_sig = [0, 1.0]
#     # W = [1.]
#     # D_sig = [1.]

#     sig_ftype = 'differentiator'

# if Ftype == 'h':

#     W = [1., 1.]
#     D = [1., 1.]

#     Be_sig = [0, 0.1, 0.9, 1.0]
#     D_sig = [1., 1.]

#     # tipo 3
#     # N = 18
#     # Be = [0.1, 0.9]
#     # tipo 4
#     # N = 19
#     # Be = [0.1, 1.]

#     sig_ftype = 'hilbert'

@pytest.mark.parametrize(
    "true_parameters, order, band_edges, desired, par_dict",
    [
       # FIR tipo 1
       (  [ np.array([ 0.00698397, -0.00919773, -0.02593083,  0.00427416,  0.04494591,
                      -0.00732263, -0.0944305 ,  0.00803432,  0.31422038,  0.49157949,
                       0.31422038,  0.00803432, -0.0944305 , -0.00732263,  0.04494591,
                       0.00427416, -0.02593083, -0.00919773,  0.00698397]),
            0.025266406063163913,
            np.array([0.     , 0.1375 , 0.25625, 0.35625, 0.4    , 0.6    , 0.63125,
                      0.70625, 0.8    , 0.9    , 1.     ])
           ],
          18,                 # order
          [0, 0.4, 0.6, 1.0], # band_edges
          [1, 1, 0, 0],       # desired
          { 'weight': [1., 3.], 'filter_type': 'multiband' } # resto
       ),
       # FIR tipo 2
       (  [ np.array([ 3.21798503e-05, -2.19529039e-02, -1.57531136e-02,  3.11212630e-02,
                       3.21290665e-02, -6.01446218e-02, -7.62058008e-02,  1.52149872e-01,
                       4.41820490e-01,  4.41820490e-01,  1.52149872e-01, -7.62058008e-02,
                      -6.01446218e-02,  3.21290665e-02,  3.11212630e-02, -1.57531136e-02,
                      -2.19529039e-02,  3.21798503e-05]),
            0.033607138093049875,
            np.array([0.        , 0.12413793, 0.24827586, 0.35172414, 0.4       ,
                      0.6       , 0.63448276, 0.72413793, 0.82758621, 0.94482759])
           ],
          17,                 # order
          [0, 0.4, 0.6, 1.0], # band_edges
          [1, 1, 0, 0],       # desired
          { 'weight': [1., 3.], 'filter_type': 'multiband' } # resto
       )
    ]
)
def test_parametrize_fir_design_pm( true_parameters, order, band_edges, desired, par_dict ):
    
    hh_mio, Err, wext = test_module.fir_design_pm(order, band_edges, desired, **par_dict)
    
    assert  np.allclose(true_parameters[0], hh_mio)
    assert  true_parameters[1] == Err
    assert  np.allclose(true_parameters[2], wext)
        

