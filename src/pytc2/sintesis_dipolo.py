#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:26:00 2023

@author: mariano
"""

import numpy as np

import sympy as sp

from .remociones import remover_polo_infinito, remover_valor_en_infinito, remover_polo_dc, remover_valor_en_dc, trim_func_s


##########################################
#%% Variables para el análisis simbólico #
##########################################

from .general import s



def cauer_RC( imm, remover_en_inf=True ):
    '''
    Description
    -----------
    Perform continued fraction expansion over imm following Cauer 2 synthesis method.

        imm = k0_0 / s + 1 / ( k0_1 + 1/ (k0_2 / s  + 1/ ... )) 

    Parameters
    ----------
    immittance : symbolic rational function
        La inmitancia a sintetizar.

    Returns
    -------
    A list k0 with the i-th k0_i resulted from continued fraction expansion.

    Ejemplo
    -------
    
    # Sea la siguiente función de excitación
    Imm = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    
    # Implementaremos Imm mediante Cauer 1 o remociones continuas en infinito
    imm_cauer_0, k0 = tc2.cauer_0(Imm)

    '''    
    
    ko = []

    if remover_en_inf:
        rem, koi = remover_polo_infinito(imm)
        bRemoverPolo = False

        if koi.is_zero:
            rem, koi = remover_valor_en_infinito(imm)
            bRemoverPolo = True
            
    else:

        rem, koi = remover_polo_dc(imm)
        bRemoverPolo = False

        if koi.is_zero:
            rem, koi = remover_valor_en_dc(imm)
            bRemoverPolo = True

    
        
    while not(rem.is_zero) and not(koi.is_zero):
        
        ko += [koi]
        rem = 1/rem

        if remover_en_inf:
            
            if bRemoverPolo:
                rem, koi = remover_polo_infinito(rem)
                bRemoverPolo = False
            else:
                rem, koi = remover_valor_en_infinito(rem)
                bRemoverPolo = True
        else:
            
            if bRemoverPolo:
                rem, koi = remover_polo_dc(rem)
                bRemoverPolo = False
            else:
                rem, koi = remover_valor_en_dc(rem)
                bRemoverPolo = True


    if koi.is_zero:
        # deshago para entender al resto de la misma 
        # naturaleza que el último elemento que retiró.
        rem = 1/rem
    else:
        ko += [koi]

    imm_as_cauer = koi
    
    for ii in np.flipud(np.arange(len(ko)-1)):

        imm_as_cauer = ko[ii] + 1/imm_as_cauer
        
    return(ko, imm_as_cauer, rem)

def cauer_LC( imm, remover_en_inf = True ):
    '''
    Description
    -----------
    Perform continued fraction expansion over imm following Cauer 1 synthesis method.

        imm = koo_0 * s + 1 / ( koo_1 * s + 1/ (koo_2 * s  + 1/ ... )) 

    Parameters
    ----------
    immittance : symbolic rational function
        La inmitancia a sintetizar.

    Returns
    -------
    A list koo with the i-th koo_i resulted from continued fraction expansion.

    Ejemplo
    -------
    
    # Sea la siguiente función de excitación
    Imm = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    
    # Implementaremos Imm mediante Cauer 1 o remociones continuas en infinito
    imm_cauer_oo, koo = tc2.cauer_oo(Imm)

    '''    
        
    rem = imm
    ko = []

    # a veces por problemas numéricos no hay cancelaciones de los términos 
    # de mayor o menor orden y quedan coeficientes muy bajos.
    rem = trim_func_s(sp.simplify(sp.expand(rem)))

    if remover_en_inf:
        rem, koi = remover_polo_infinito(rem)
    else:
        rem, koi = remover_polo_dc(rem)
        
    
    while not(rem.is_zero) and not(koi.is_zero):
        
        ko += [koi]
        rem = 1/rem

        # a veces por problemas numéricos no hay cancelaciones de los términos 
        # de mayor o menor orden y quedan coeficientes muy bajos.
        rem = trim_func_s(sp.simplify(sp.expand(rem)))

        if remover_en_inf:
            rem, koi = remover_polo_infinito(rem)
        else:
            rem, koi = remover_polo_dc(rem)

    if koi.is_zero:
        # deshago para entender al resto de la misma 
        # naturaleza que el último elemento que retiró.
        rem = 1/rem
    else:
        ko += [koi]

    imm_as_cauer = ko[-1] + rem

    for ii in np.flipud(np.arange(len(ko)-1)):
        
        imm_as_cauer = ko[ii] + 1/imm_as_cauer
        
    return(ko, imm_as_cauer, rem)



def foster_zRC2yRC( k0 = None, koo = None, ki_wi = None, kk = None, ZRC_foster = None ):
    '''
    Parameters
    ----------
    immittance : symbolic rational function
        La inmitancia a sintetizar.

    Returns
    -------
    Convierte una expansión disipativa de foster ZRC a YRC, cuando se 
    expande YRC/s para que quede de la forma ZRC.
        
        Imm = k0 / s + koo * s +  1 / ( k0_i / s + koo_i * s ) 


    imm_list = [ k0, koo, [k00, koo0], [k01, koo1], ..., [k0N, kooN]  ]
    
    Si algún elemento no está presente, su valor será de "None".

    Ejemplo
    -------
    
    # Sea la siguiente función de excitación
    Imm = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    
    # Implementaremos Imm mediante Foster
    k0, koo, ki = tc2.foster(Imm)


    '''    
    
    if koo is None:
    # koo tiene que ser None para ZRC ya que en inf habrá
    # o 0 o cte.
        
        if not(kk is None):
            koo = kk
            kk = None
            
        if not(k0 is None):
            kk = k0
            k0 = None
            
        if not(ki_wi is None):
            
            ki = ki_wi
            # ki = []
            # for this_ki_wi in ki_wi:
                
            #     ki += [[this_ki_wi[1], this_ki_wi[0]]]
            
            
        YRC_foster = sp.expand(ZRC_foster * s)

    return([k0, koo, ki, kk, YRC_foster])


def foster( imm ):
    '''
    Parameters
    ----------
    immittance : symbolic rational function
        La inmitancia a sintetizar.

    Returns
    -------
    Una lista imm_list con los elementos obtenidos de la siguientes expansión en 
    fracciones simples:
        
        Imm = k0 / s + koo * s +  1 / ( k0_i / s + koo_i * s ) 


    imm_list = [ k0, koo, [k00, koo0], [k01, koo1], ..., [k0N, kooN]  ]
    
    Si algún elemento no está presente, su valor será de "None".

    Ejemplo
    -------
    
    # Sea la siguiente función de excitación
    Imm = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    
    # Implementaremos Imm mediante Foster
    k0, koo, ki = tc2.foster(Imm)


    '''    

    num, den = imm.as_numer_denom()
    
    # grados de P y Q
    deg_P = sp.degree(num)
    deg_Q = sp.degree(den)
    
    imm_foster = sp.polys.partfrac.apart(imm)
    
    all_terms = imm_foster.as_ordered_terms()
    
    kk = None
    k0 = None
    koo = None
    ki = []
    ii = 0
    
    foster_form = sp.Rational(0)
    
    for this_term in all_terms:

        foster_form += this_term
        
        num, den = this_term.as_numer_denom()
        
        if sp.degree(num) == 1 and sp.degree(den) == 0:
        
            koo = num.as_poly(s).LC() / den
    
        elif sp.degree(den) == 1 and sp.degree(num) == 0:
            
            if den.as_poly(s).all_coeffs()[1] == 0:
                # red no disipativa
                k0 = num / den.as_poly(s).LC()
            else:
                # red disipativa - tanque RC-RL
                
                # kk_i, koo_i
                ki += [[(den / num).expand().as_poly(s).EC(), 
                        (den / num).expand().as_poly(s).LC() ]]
                ii += 1
    
        elif sp.degree(den) == 0 and sp.degree(num) == 0:
            # constante en redes disipativas
            kk = num / den.as_poly(s).LC()
    
        elif sp.degree(num) == 1 and sp.degree(den) == 2:
            # tanque
            tank_el = (den / num).expand().as_ordered_terms()
    
            koo_i = None
            k0_i = None
            
            for this_el in tank_el:
                
                num, den = this_el.as_numer_denom()
                
                if sp.degree(num) == 1 and sp.degree(den) == 0:
                
                    koo_i = num.as_poly(s).LC() / den
    
                elif sp.degree(den) == 1 and sp.degree(num) == 0:
                    
                    k0_i = num / den.as_poly(s).LC() 
                    
            
            ki += [[k0_i, koo_i]]
            ii += 1
            
        else:
            # error
            assert('Error al expandir en fracciones simples.')
    
    if ii == 0:
        ki = None

    return([k0, koo, ki, kk, foster_form])


