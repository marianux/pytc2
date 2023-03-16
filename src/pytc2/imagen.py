#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:11:38 2023

@author: mariano
"""


import numpy as np

import sympy as sp



'''
    Bloque de funciones para parametros imagen
'''


def db2nepper(at_en_db):
    '''
    Convierte una matriz de parámetros scattering (S) simbólica 
    al modelo de parámetros transferencia de scattering (Ts).

    Parameters
    ----------
    Spar : Symbolic Matrix
        Matriz de parámetros S.

    Returns
    -------
    Ts : Symbolic Matrix
        Matriz de parámetros de transferencia scattering.

    '''
    
    return( at_en_db/(20*np.log10(np.exp(1))) )

def nepper2db(at_en_np):
    '''
    Convierte una matriz de parámetros scattering (S) simbólica 
    al modelo de parámetros transferencia de scattering (Ts).

    Parameters
    ----------
    Spar : Symbolic Matrix
        Matriz de parámetros S.

    Returns
    -------
    Ts : Symbolic Matrix
        Matriz de parámetros de transferencia scattering.

    '''
    
    return( at_en_np*(20*np.log10(np.exp(1))) )

    

def I2T(gamma, z01, z02 = None):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    if z02 is None:
        z02 = z01

    # if np.sqrt(z02/z01)
    
    TT = np.matrix([[np.cosh(gamma)*np.sqrt(z01/z02),
                     np.sinh(gamma)*np.sqrt(z01*z02)], 
                    [np.sinh(gamma)/np.sqrt(z01*z02),
                     np.cosh(gamma)*np.sqrt(z02/z01)]])
    
    return(TT)


def I2T_s(gamma, z01, z02 = None):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    if z02 is None:
        z02 = z01
    
    TT = sp.Matrix([[sp.cosh(gamma)*sp.sqrt(z01/z02),
                     sp.sinh(gamma)*sp.sqrt(z01*z02)], 
                    [sp.sinh(gamma)/sp.sqrt(z01*z02),
                     sp.cosh(gamma)*sp.sqrt(z02/z01)]])
    
    
    return(TT)
