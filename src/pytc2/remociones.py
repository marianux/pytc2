#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:12:53 2023

@author: mariano
"""

import numpy as np

import sympy as sp

##########################################
#%% Variables para el análisis simbólico #
##########################################

from .general import s




def tanque_z( doska, omegasq ):
    '''
    Calcula los valores de L y C que componen un tanque resonante LC 
    (tanque Z), a partir del valor del residuo ($ k $) y la omega al cuadrado 
    ($ \omega^2 $) de la expresión de impedancia dada por:
        
        $$ Z_{LC} = \frac{2.k.s}{(s^2+\omega^2)} $$

    Parameters
    ----------
    doska : Symbolic
        Dos veces el residuo.
    omegasq : Symbolic
        Cuadrado de la omega a la que el tanque resuena.

    Returns
    -------
    L : Symbolic
        Valor de la admitancia
    C : Symbolic
        Valor de la capacidad

    '''
    
    return( [doska/omegasq, 1/doska] )

def tanque_y( doska, omegasq ):
    '''
    Calcula los valores de L y C que componen un tanque resonante LC 
    (tanque Z), a partir del valor del residuo ($ k $) y la omega al cuadrado 
    ($ \omega^2 $) de la expresión de impedancia dada por:
        
        $$ Y_{LC} = \frac{2.k.s}{(s^2+\omega^2)} $$

    Parameters
    ----------
    doska : Symbolic
        Dos veces el residuo.
    omegasq : Symbolic
        Cuadrado de la omega a la que el tanque resuena.

    Returns
    -------
    L : Symbolic
        Valor de la admitancia
    C : Symbolic
        Valor de la capacidad

    '''
    
    return( [1/doska, doska/omegasq] )



def trim_poly_s( this_poly, tol = 10**-6 ):
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

    all_terms = this_poly.as_poly(s).all_terms()
    
    poly_acc = 0
    
    for this_pow, this_coeff in all_terms:
    
        if np.abs(this_coeff) > tol:
            
            poly_acc = poly_acc + this_coeff * s**this_pow[0]


    return(poly_acc)

def trim_func_s( rat_func, tol = 10**-6 ):
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

    num, den = rat_func.as_numer_denom()
    
    num = trim_poly_s(num, tol)
    den = trim_poly_s(den, tol)
    
    return(num/den)

def modsq2mod_s( aa ):
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

    num, den = sp.fraction(aa)

    k = sp.poly(num,s).LC() / sp.poly(den,s).LC()
    
    # roots_num = num.as_poly(s).all_roots()
    
    # poly_acc = sp.Rational('1')
    
    roots_num = sp.roots(num)
    
    poly_acc = sp.Rational('1')

    for this_root in roots_num.keys():
        
        if sp.re(this_root) <= 0:
            
            # multiplicidad
            mult = roots_num[this_root]

            if mult > 1:
                poly_acc *= (s-this_root)**sp.Rational(mult/2)
            else:
                poly_acc *= (s-this_root)

    assert (len(num.as_poly(s).all_coeffs())-1)/2 == len(poly_acc.as_poly(s).all_coeffs())-1, 'Falló la factorización de modsq2mod_s. ¡Revisar!'

    num = sp.simplify(sp.expand(poly_acc))



    # poly_acc = 0
    
    # for each_term in num.as_poly(s).all_terms():
        
    #     poly_acc += np.abs(each_term[1]) * s**each_term[0][0]

    # num = poly_acc
    
    roots_num = sp.roots(den)
    
    poly_acc = sp.Rational('1')

    for this_root in roots_num.keys():
        
        if sp.re(this_root) <= 0:
            
            # multiplicidad
            mult = roots_num[this_root]

            if mult > 1:
                poly_acc *= (s-this_root)**sp.Rational(mult/2)
            else:
                poly_acc *= (s-this_root)
    
    den = sp.simplify(sp.expand(poly_acc))

    poly_acc = 0
    
    for each_term in den.as_poly(s).all_terms():
        
        poly_acc += np.abs(each_term[1]) * s**each_term[0][0]

    den = poly_acc

    return(sp.simplify(sp.expand(sp.sqrt(k) * num/den))) 


def modsq2mod( aa ):
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
    
    rr = np.roots(aa)
    bb = rr[np.real(rr) == 0]
    bb = bb[ :(bb.size//2)]
    bb = np.concatenate( [bb, rr[np.real(rr) < 0]])
    
    return np.flip(np.real(np.polynomial.polynomial.polyfromroots(bb)))


'''
    Bloque de funciones para la síntesis gráfica de imitancias
'''

def remover_polo_sigma( imm, sigma, isImpedance = True,  isRC = True,  sigma_zero = None ):
    '''
    Se removerá el residuo en sobre el eje $\sigma$ (sigma) de la impedancia (zz) 
    o admitancia (yy) de forma completa, o parcial en el caso que se especifique una 
    sigma_i.
    Como resultado de la remoción, quedará otra función racional definida
    como:
        
    $$ Z_{R}= Z - \frac{k_i}{s + \sigma_i} $$
    
    siendo 

    $$ k=\lim\limits _{s\to -\sigma_i} Z (s + \sigma_i) $$
    
    En cuanto se especifique sigma_i, la remoción parcial estará definida 
    como

    $$ Z_{R}\biggr\rfloor_{s=-\sigma_i}= 0 = Z - \frac{k_i}{s + \sigma_i}\biggr\rfloor_{s=-\sigma_i} $$
    
    siendo 
    
    $$ k = Z.(\frac{)s + \sigma_i)\biggr\rfloor_{s=-\sigma_i} $$
    

    Parameters
    ----------
    zz o yy: Symbolic
        Impedancia o admitancia que se utilizará para la remoción. Es una función racional 
        simbólica que tendrá un polo de orden 1 en \omega.
    omega_zero : Symbolic
        Frecuencia a la que la imitancia será cero luego de la remoción.

    Returns
    -------
    imit_r : Symbolic
        Imitancia luego de la remoción
    k : Symbolic
        Valor del residuo.
    '''

    if isImpedance:
        zz = imm
    else:
        yy = imm

    if sigma_zero is None:
        # remoción total
        
        if isImpedance:
            if isRC:
                kk = sp.limit(zz*(s + sigma), s, -sigma)
            else:
                # RL
                kk = sp.limit(zz*(s + sigma)/s, s, -sigma)
                
        else:
            if isRC:
                kk = sp.limit(yy*(s + sigma)/s, s, -sigma)
            else:
                kk = sp.limit(yy*(s + sigma), s, -sigma)
        
        if kk.is_negative:
            assert('Residuo negativo. Verificar Z/Y RC/RL')
        
    else:
        # remoción parcial
        if isImpedance:
            if isRC:
                kk = sp.simplify(sp.expand(zz*(s + sigma))).subs(s, -sigma_zero)
            else:
                kk = sp.simplify(sp.expand(zz*(s + sigma)/s)).subs(s, -sigma_zero)
            
        else:
            if isRC:
                kk = sp.simplify(sp.expand(yy*(s + sigma)/s)).subs(s, -sigma_zero)
            else:
                kk = sp.simplify(sp.expand(yy*(s + sigma))).subs(s, -sigma_zero)

        if kk.is_negative:
            assert('Residuo negativo. Verificar Z/Y RC/RL')
    
    # extraigo kk
    if isImpedance:
        if isRC:
            # Z_RC        
            R = kk/sigma
            CoL = 1/kk
            kk  = kk/(s+sigma)
        else:
            # Z_RL        
            R = kk
            CoL = kk/sigma
            kk  = kk*s/(s+sigma)
        
    else:

        if isRC:
            # Y_RC        
            CoL = kk/sigma
            R = 1/kk
            kk  = kk*s/(s+sigma)
        else:
            # Y_RL
            R = sigma/kk
            CoL = 1/kk
            kk  = kk/(s+sigma)
        

    if isImpedance:
        imit_r = sp.factor(sp.simplify(sp.expand(zz - kk)))
    
    else:
    
        imit_r = sp.factor(sp.simplify(sp.expand(yy - kk)))

    return( [imit_r, kk, R, CoL] )

def remover_polo_jw( imit, omega = None , isImpedance = True, omega_zero = None ):
    '''
    Se removerá el residuo en sobre el eje $j.\omega$ (omega) de la imitancia 
    $I$ (imit) de forma completa, o parcial en el caso que se especifique una 
    omega_zero.
    Como resultado de la remoción, quedará otra función racional definida
    como:
        
    $$ I_{R}=I-\frac{2.k.s}{s^{2}+\omega^{2}} $$
    
    siendo 

    $$ k=\lim\limits _{s^2\to-\omega^2}I\frac{2.k.s}{s^{2}+\omega^{2}} $$
    
    En cuanto se especifique omega_zero, la remoción parcial estará definida 
    como

    $$ I_{R}\biggr\rfloor_{s^{2}=-\omega_{z}^{2}}=0=I-\frac{2.k.s}{s^{2}+\omega^{2}}\biggr\rfloor_{s^{2}=-\omega_{z}^{2}} $$
    
    siendo 
    
    $$ 2.k^{'}=I.\frac{s^{2}+\omega^{2}}{s}\biggr\rfloor_{s^{2}=-\omega_z^{2}} $$
    

    Parameters
    ----------
    imit : Symbolic
        Imitancia que se utilizará para la remoción. Es una función racional 
        simbólica que tendrá un polo de orden 1 en \omega.
    omega_zero : Symbolic
        Frecuencia a la que la imitancia será cero luego de la remoción.

    Returns
    -------
    imit_r : Symbolic
        Imitancia luego de la remoción
    k_inf : Symbolic
        Valor del residuo en infinito
    '''

    if omega is None:
        # busco el primer polo finito en imit sobre el jw
        
        _, den = (imit).as_numer_denom()
        faux = sp.factor_list(den)
        
        if sp.degree(faux[1][0][0]) == 2:
            
            tt = faux[1][0][0].as_ordered_terms()
            
            # el último término sería omega**2. Cada factor sería
            # s**2 + omega**2
            omega = sp.sqrt(tt[-1])

    if omega_zero is None:
        # remoción total
        # kk = sp.limit(imit*(s**2+omega**2)/s, s**2, -omega**2)
        kk = sp.simplify(sp.expand(imit*(s**2+omega**2)/s)).subs(s**2, -(omega**2) )
        
    else:
        # remoción parcial
        kk = sp.simplify(sp.expand(imit*(s**2+omega**2)/s)).subs(s**2, -(omega_zero**2) )

    
    if isImpedance:
        # Z_LC
        L = kk/omega**2
        C = 1/kk
        
    else:
        # Y_LC
        C = kk/omega**2
        L = 1/kk

    kk = kk * s / (s**2+omega**2)
    
    # extraigo kk
    imit_r = sp.factor(sp.simplify(sp.expand(imit - kk)))

    return( [imit_r, kk, L, C] )

def remover_polo_dc( imit, omega_zero = None ):
    '''
    Se removerá el residuo en continua (s=0) de la imitancia ($I$) de forma 
    completa, o parcial en el caso que se especifique una omega_zero. 
    Como resultado de la remoción, quedará otra función racional definida
    como:
        
    $$ I_R = I - k_0/s  $$
    
    siendo 

    $$ k_0=\lim\limits _{s\to0}I.s $$
    
    En cuanto se especifique omega_zero, la remoción parcial estará definida 
    como

    $$ I_{R}\biggr\rfloor_{s^{2}=-\omega_z^{2}}=0=I-s.k_{0}^{'}\biggr\rfloor_{s^{2}=-\omega_z^{2}} $$
    
    siendo 
    
    $$ k_{0}^{'}=I.s\biggr\rfloor_{s^{2}=-\omega_z^{2}} $$
    

    Parameters
    ----------
    imit : Symbolic
        Imitancia que se utilizará para la remoción. Es una función racional 
        simbólica que tendrá un polo de orden 1 en 0, es decir la 
        diferencia de grados entre num y den será exactamente -1.
    omega_zero : Symbolic
        Frecuencia a la que la imitancia será cero luego de la remoción.

    Returns
    -------
    imit_r : Symbolic
        Imitancia luego de la remoción
    k_inf : Symbolic
        Valor del residuo en infinito
    '''

    if omega_zero is None:
        # remoción total
        k_cero = sp.limit(imit*s, s, 0)
        
    else:
        # remoción parcial
        k_cero = sp.simplify(sp.expand(imit*s)).subs(s**2, -(omega_zero**2) )

    k_cero = k_cero/s
    
    # extraigo C3
    imit_r = sp.factor(sp.simplify(sp.expand(imit - k_cero)))

    return( [imit_r, k_cero] )

def remover_polo_infinito( imit, omega_zero = None ):
    '''
    Se removerá el residuo en infinito de la imitancia ($I$) de forma 
    completa, o parcial en el caso que se especifique una omega_zero. 
    Como resultado de la remoción, quedará otra función racional definida
    como:
        
    $$ I_R = I - s.k_\infty  $$
    
    siendo 

    $$ k_{\infty}=\lim\limits _{s\to\infty}I.\nicefrac{1}{s} $$
    
    En cuanto se especifique omega_zero, la remoción parcial estará definida 
    como

    $$ I_{R}\biggr\rfloor_{s^{2}=-\omega_z^{2}}=0=I-s.k_{\infty}^{'}\biggr\rfloor_{s^{2}=-\omega_z^{2}} $$
    
    siendo 
    
    $$ k_{\infty}^{'}=I.\nicefrac{1}{s}\biggr\rfloor_{s^{2}=-\omega_z^{2}} $$
    

    Parameters
    ----------
    imit : Symbolic
        Imitancia que se utilizará para la remoción. Es una función racional 
        simbólica que tendrá un polo de orden 1 en infinito, es decir la 
        diferencia de grados entre num y den será exactamente 1.
    omega_zero : Symbolic
        Frecuencia a la que la imitancia será cero luego de la remoción.

    Returns
    -------
    imit_r : Symbolic
        Imitancia luego de la remoción
    k_inf : Symbolic
        Valor del residuo en infinito
    '''

    if omega_zero is None:
        # remoción total
        k_inf = sp.limit(imit/s, s, sp.oo)
        
    else:
        # remoción parcial
        k_inf = sp.simplify(sp.expand(imit/s)).subs(s**2, -(omega_zero**2) )

    k_inf = k_inf * s

    # extraigo C3
    imit_r = sp.factor(sp.simplify(sp.expand(imit - k_inf)))

    return( [imit_r, k_inf] )

def remover_valor( imit, sigma_zero):
    '''
    Se removerá un valor constante de la imitancia ($I$) de forma 
    que al removerlo, la imitancia luego de la remoción ($I_R$) tenga 
    un cero en sigma_zero. Es decir:

    $$ I_{R}\biggr\rfloor_{s = -\sigma_z} = 0 = (I - k_{\infty}^{'})\biggr\rfloor_{s = -\sigma_z} $$
    
    siendo 
    
    $$ k_{\infty}^{'}= I\biggr\rfloor_{s = -\sigma_z} $$

    Parameters
    ----------
    imit : Symbolic
        Imitancia que se utilizará para la remoción. Es una función racional 
        simbólica que tendrá un valor constante en infinito (mayor a su valor en s=0).
        
    omega_zero : Symbolic
        Frecuencia a la que la imitancia será cero luego de la remoción.

    Returns
    -------
    imit_r : Symbolic
        Imitancia luego de la remoción
    k_inf : Symbolic
        Valor del residuo en infinito
    '''

    # remoción parcial
    k_prima = sp.simplify(sp.expand(imit)).subs(s, -sigma_zero)
    
    # extraigo k_prima
    imit_r = sp.factor(sp.simplify(sp.expand(imit - k_prima)))

    return( [imit_r, k_prima] )

def remover_valor_en_infinito( imit ):
    '''
    Se removerá un valor constante en infinito de la imitancia ($I$) de forma 
    completa. 
    Como resultado de la remoción, quedará otra función racional definida
    como:
        
    $$ I_R = I - k_{\infty}  $$
    
    siendo 

    $$ k_{\infty}=\lim\limits _{s\to\infty}I $$

    Parameters
    ----------
    imit : Symbolic
        Imitancia que se utilizará para la remoción. Es una función racional 
        simbólica que tendrá un valor constante en infinito (mayor a su valor en s=0).

    Returns
    -------
    imit_r : Symbolic
        Imitancia luego de la remoción
    k_inf : Symbolic
        Valor del residuo en infinito
    '''

    # remoción total
    k_inf = sp.limit(imit, s, sp.oo)

    # extraigo k_inf
    imit_r = sp.factor(sp.simplify(sp.expand(imit - k_inf)))

    return( [imit_r, k_inf] )

def remover_valor_en_dc( imit ):
    '''
    Se removerá un valor constante en continua (s=0) de la imitancia ($I$) de forma 
    completa. 
    Como resultado de la remoción, quedará otra función racional definida
    como:
        
    $$ I_R = I - k_0  $$
    
    siendo 

    $$ k_0 = \lim\limits _{s \to 0}I $$
    
    Parameters
    ----------
    imit : Symbolic
        Imitancia que se utilizará para la remoción. Es una función racional 
        simbólica que tendrá un valor constante en infinito (mayor a su valor en s=0).
        
    Returns
    -------
    imit_r : Symbolic
        Imitancia luego de la remoción
    k_inf : Symbolic
        Valor del residuo en infinito
    '''

    # remoción total
    k0 = sp.limit(imit, s, 0)
        
    # extraigo k0
    imit_r = sp.factor(sp.simplify(sp.expand(imit - k0)))

    return( [imit_r, k0] )

