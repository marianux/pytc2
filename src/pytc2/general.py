#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:22:31 2023

Originally based on the work of Combination of 2011 Christopher Felton
Further modifications were added for didactic purposes
by Mariano Llamedo llamedom _at_ frba_utn_edu_ar

@author: marianux
"""

import sympy as sp


##########################################
#%% Variables para el análisis simbólico #
##########################################

# Laplace complex variable. s = σ + j.ω
s = sp.symbols('s', complex=True)
# Fourier real variable ω 
w = sp.symbols('w', complex=False)


def pp(z1, z2):
    '''
    Asocia en paralelo dos impedancias o en serie dos admitancias.

    Parameters
    ----------
    z1 : Symbolic
        Impedancia 1.
    z2 : Symbolic
        Impedancia 2.

    Returns
    -------
    zp : Symbolic
         Impedancia resultante.

    '''

    return(z1*z2/(z1+z2))
    




'''
Otras funciones

'''


def Chebyshev_polynomials(nn):
    
    Cn_pp = 1
    Cn_p = w
    
    if nn > 1:
        
        for ii in range(nn-1):
            
            Cn = 2 * w * Cn_p - Cn_pp
    
            Cn_pp = Cn_p
            Cn_p = Cn

    elif nn == 1:
        Cn = Cn_p
        
    else:
        Cn = 1
            
    return(sp.simplify(sp.expand(Cn)))

