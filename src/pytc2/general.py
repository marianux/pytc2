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
import numpy as np

from IPython.display import display, Math, Markdown

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

#########################
#%% Funciones generales #
#########################

def print_console_alert(strAux):
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
    
    strAux = '# ' + strAux + ' #\n'
    strAux1 =  '#' * (len(strAux)-1) + '\n' 
    
    print( '\n\n' + strAux1 + strAux + strAux1 )
    
def print_console_subtitle(strAux):
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
    
    strAux = strAux + '\n'
    strAux1 =  '-' * (len(strAux)-1) + '\n' 
    
    print( '\n\n' + strAux + strAux1 )
    
def print_subtitle(strAux):
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
    
    display(Markdown('#### ' + strAux))


def a_equal_b_latex_s( a, b):
    '''
    Convierte un símbolo en un string formateado para visualizarse en LaTex 

    Parameters
    ----------
    Spar : Symbolic Matrix
        Matriz de parámetros S.

    Returns
    -------
    Ts : Symbolic Matrix
        Matriz de parámetros de transferencia scattering.

    '''

    if isinstance(a, sp.Basic ):
        a_str = sp.latex(a)
    elif isinstance(a, str):
        a_str = a
    else:
        a_str = '??'
    
    
    return('$' + a_str + '=' + sp.latex(b) + '$')

def to_latex( unsimbolo ):
    '''
    Convierte un símbolo en un string formateado para visualizarse en LaTex 

    Parameters
    ----------
    Spar : Symbolic Matrix
        Matriz de parámetros S.

    Returns
    -------
    Ts : Symbolic Matrix
        Matriz de parámetros de transferencia scattering.

    '''
    
    return('$'+ sp.latex(unsimbolo) + '$')

def str_to_latex( unstr):
    '''
    Formatea un string para visualizarse en LaTex 

    Parameters
    ----------
    Spar : Symbolic Matrix
        Matriz de parámetros S.

    Returns
    -------
    Ts : Symbolic Matrix
        Matriz de parámetros de transferencia scattering.

    '''
    
    return('$'+ unstr + '$')



def print_latex(strAux):
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
    
    display(Math(strAux))

def Chebyshev_polynomials(nn):
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



'''
  ################################################
 ## Bloque de funciones para parametros imagen ##
################################################
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

    