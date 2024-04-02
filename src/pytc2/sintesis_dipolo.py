#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:26:00 2023

@author: mariano
"""

import numpy as np

import sympy as sp

from .remociones import isFRP, remover_polo_infinito, remover_valor_en_infinito, remover_polo_dc, remover_valor_en_dc, trim_func_s


##########################################
#%% Variables para el análisis simbólico #
##########################################

from .general import s, expr_simb_expr, print_console_alert, print_latex
    

def cauer_RC( imm, remover_en_inf=True ):
    '''
    Realiza una expansión en fracciones continuas sobre una inmitancia (imm), 
    removiendo en DC o :math:`\\infty` dependiendo de (remover_en_inf). Este 
    procedimiento se conoce como métodos de Cauer I y II. En el ejemplo de 
    :math:`Z_{RC}` se remueve en DC y para el caso de :math:`Y_{RC}` 
    en :math:`\\infty`.

    .. math:: Z_{RC}(s)= \\frac{1}{s.C_1} + \\frac{1}{ \\frac{1}{R_1} + \\frac{1}{ \\frac{1}{s.C_2} + \\cdots } } = 
         R_1 + \\frac{1}{ s.C_1 + \\frac{1}{ R_2 + \\cdots } } 

    .. math:: Y_{RC}(s)= s.C_1 + \\frac{1}{ R_1 + \\frac{1}{ s.C_2 + \\cdots } } = 
         \\frac{1}{R_1} + \\frac{1}{ s.C_1 + \\frac{1}{ \\frac{1}{R_2} + \\cdots } } 

    Parameters
    ----------
    imm : symbolic rational function
        La inmitancia a expandir en fracciones continuas..
    remover_en_inf : boolean
        Determina en qué extremo se realiza la remoción.
        

    Returns
    -------
    A list k0 with the i-th k0_i resulted from continued fraction expansion.


    Raises
    ------
    ValueError
        Si y_exc y z_exc no son una instancia de sympy.Expr.


    See Also
    --------
    :func:`cauer_LC`
    :func:`dibujar_cauer_RC_RL`
    :func:`dibujar_cauer_LC`


    Example
    -------
    >>> import sympy as sp
    >>> from pytc2.sintesis_dipolo import cauer_RC
    >>> from pytc2.dibujar import dibujar_cauer_RC_RL
    >>> s = sp.symbols('s ', complex=True)
    >>> # Sea la siguiente función de excitación
    >>> ZRC = (s**2 + 4*s + 3)/(s**2 + 2*s)
    >>> # Implementaremos FF mediante Cauer 1 o remociones continuas en infinito
    >>> koo, ZRC_cauer_oo, rem = cauer_RC(ZRC, remover_en_inf=True)
    >>> # Tratamos a nuestra función inmitancia como una Z
    >>> dibujar_cauer_RC_RL(koo, z_exc = ZRC_cauer_oo)
    >>> # Tratamos a nuestra función inmitancia como una Y
    >>> dibujar_cauer_RC_RL(koo, y_exc = ZRC_cauer_oo)

    '''    

    if not isinstance(imm , sp.Expr):
        raise ValueError('Hay que definir imm como una expresión simbólica.')

    if not isinstance(remover_en_inf, bool):
        raise ValueError('remover_en_inf debe ser un booleano.')
    
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

    
    if isFRP(rem):
        
        bFRP = True
        
        while bFRP and not(rem.is_zero) and not(koi.is_zero):
            
            
            ko += [koi]
            rem = 1/rem
    
            if remover_en_inf:
                
                if bRemoverPolo:
                    rem_aux, koi = remover_polo_infinito(rem)
                    bRemoverPolo = False
                else:
                    rem_aux, koi = remover_valor_en_infinito(rem)
                    bRemoverPolo = True
            else:
                
                if bRemoverPolo:
                    rem_aux, koi = remover_polo_dc(rem)
                    bRemoverPolo = False
                else:
                    rem_aux, koi = remover_valor_en_dc(rem)
                    bRemoverPolo = True

            bFRP = isFRP(rem_aux)

            if bFRP:
                rem = rem_aux
   
            
        if koi.is_zero:
            # deshago para entender al resto de la misma 
            # naturaleza que el último elemento que retiró.
            rem = 1/rem
        else:
            if bFRP:
                ko += [koi]
    
        imm_as_cauer = ko[-1]
        
        for ii in np.flipud(np.arange(len(ko)-1)):
    
            imm_as_cauer = ko[ii] + 1/imm_as_cauer

    else:
        # no se pudo hacer ninguna remoción
        imm_as_cauer = imm
        ko = [s*0]
        rem = s*0
            
    if not (sp.simplify(sp.expand(imm_as_cauer - imm))).is_zero:
        # error
        print_console_alert('Fallo la expansión')
        print_latex(expr_simb_expr(imm, imm_as_cauer, ' \\neq '))
        RuntimeWarning('Fallo la expansión Cauer. Revisar!!')

    return(ko, imm_as_cauer, rem)

def cauer_LC( imm, remover_en_inf = True ):
    '''
    Dibuja una red escalera no disipativa, a partir de la expansión en fracciones 
    continuas (Método de Cauer). Se remueve en DC o :math:`\\infty` dependiendo 
    de *remover_en_inf*. En los siguientes ejemplos se expande tanto :math:`Z(s)`
    como :math:`Y(s)`, y se remueve a la izquierda en DC y a la derecha en :math:`\\infty`.
    La forma matemática será:

    .. math:: Z(s)= \\frac{1}{s.C_1} + \\frac{1}{ \\frac{1}{s.L_1} + \\frac{1}{ \\frac{1}{s.C_2} + \\cdots } } = 
             s.L_1 + \\frac{1}{ s.C_1 + \\frac{1}{ s.L_2 + \\cdots } } 

    .. math:: Y(s)= \\frac{1}{s.L_1} + \\frac{1}{ \\frac{1}{s.C_1} + \\frac{1}{ \\frac{1}{s.L_2} + \\cdots } } = 
             s.C_1 + \\frac{1}{ s.L_1 + \\frac{1}{ s.C_2 + \\cdots } }  


    Parameters
    ----------
    imm : symbolic rational function
        La inmitancia a expandir en fracciones continuas..
    remover_en_inf : boolean
        Determina en qué extremo se realiza la remoción.
        

    Returns
    -------
    ko : lista de expresiones simbólicas
        Conjunto de términos con los residuos de forma :math:`\\frac{k_0}{s}` y :math:`s.k_{\\infty}`
    imm_as_cauer : symbolic rational function
        La función inmitancia expandida en fracciones contínuas.
    rem : symbolic rational function
        0 en caso que la expansión sea exitosa, ó una función remanente que no 
        puede ser expresada en formato Cauer.


    Raises
    ------
    ValueError
        Si y_exc y z_exc no son una instancia de sympy.Expr.


    See Also
    --------
    :func:`cauer_LC`
    :func:`foster_zRC2yRC`
    :func:`dibujar_cauer_LC`

    
    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.sintesis_dipolo import cauer_LC
    >>> from pytc2.dibujar import dibujar_cauer_LC
    >>> s = sp.symbols('s ', complex=True)
    >>> # Sea la siguiente función de excitación
    >>> FF = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    >>> # Implementaremos FF mediante Cauer 1 o remociones continuas en infinito
    >>> koo, F_cauer_oo, rem = cauer_LC(FF, remover_en_inf=True)
    >>> # Tratamos a nuestra función inmitancia como una Z
    >>> dibujar_cauer_LC(koo, z_exc = F_cauer_oo)
    >>> # Tratamos a nuestra función inmitancia como una Y
    >>> dibujar_cauer_LC(koo, y_exc = F_cauer_oo)

    '''    

    if not isinstance(imm , sp.Expr):
        raise ValueError('Hay que definir imm como una expresión simbólica.')

    if not isinstance(remover_en_inf, bool):
        raise ValueError('remover_en_inf debe ser un booleano.')
    
        
    rem = imm
    ko = []

    # a veces por problemas numéricos no hay cancelaciones de los términos 
    # de mayor o menor orden y quedan coeficientes muy bajos.
    rem = trim_func_s(sp.simplify(sp.expand(rem)))

    if remover_en_inf:
        rem_aux, koi = remover_polo_infinito(rem)
    else:
        rem_aux, koi = remover_polo_dc(rem)

    bFRP = isFRP(rem_aux)

    if bFRP:
        
        rem = rem_aux
    
        while bFRP and not(rem.is_zero) and not(koi.is_zero):
            
            ko += [koi]
            rem = 1/rem
    
            # a veces por problemas numéricos no hay cancelaciones de los términos 
            # de mayor o menor orden y quedan coeficientes muy bajos.
            rem = trim_func_s(sp.simplify(sp.expand(rem)))
    
            if remover_en_inf:
                rem_aux, koi = remover_polo_infinito(rem)
            else:
                rem_aux, koi = remover_polo_dc(rem)
    
            bFRP = isFRP(rem_aux)

            if bFRP:
                rem = rem_aux

        if koi.is_zero:
            # deshago para entender al resto de la misma 
            # naturaleza que el último elemento que retiró.
            rem = 1/rem
        else:
            if bFRP:
                # si no salimos por rem NO FRP
               ko += [koi]
    
        imm_as_cauer = ko[-1] + rem
    
        for ii in np.flipud(np.arange(len(ko)-1)):
            
            imm_as_cauer = ko[ii] + 1/imm_as_cauer

    else:
        # no se pudo hacer ninguna remoción
        imm_as_cauer = imm
        ko = [s*0]
        rem = s*0
        
    if not (sp.simplify(sp.expand(imm_as_cauer - imm))).is_zero:
        # error
        print_console_alert('Fallo la expansión')
        print_latex(expr_simb_expr(imm, imm_as_cauer, ' \\neq '))
        
    return(ko, imm_as_cauer, rem)

# TODO: me gustaría documentar y probar mejor esta función
def foster_zRC2yRC( k0 = sp.Rational(0), koo = sp.Rational(0), ki_wi = sp.Rational(0), kk = sp.Rational(0), ZRC_foster = sp.Rational(0) ):
    '''
    Permite llegar a la forma foster de una inmitancia :math:`I(s)` (YRC - ZRL), 
    a partir de la propia función :func:`foster` de expansión en fracciones 
    simples, y una conversión término-a-término de cada residuo obtenido. 
    
    De esa manera se comienza con la expansión foster( I(s)/s ), para luego 
    realizarel siguiente mapeo de residuos:

    + :math:`k_\\infty = kk` 
    + :math:`k_k = k_0` 
    + :math:`k_i = ki_wi` 
    + :math:`I_F(s) = I(s)*s` 

    Parameters
    ----------
    k0:  simbólica, opcional
        Residuo de la función en DC o :math:`s \\to 0`. El valor predeterminado es 0.
    koo:  simbólica, opcional
        Residuo de la función en infinito o :math:`s \\to \\infty`. El valor predeterminado es 0.
    ki_wi:  simbólica, list o tuple opcional
        Residuo de la función en :math:`\\omega_i` o :math:`s^2 \\to -\\omega^2_i`. El valor predeterminado es 0.
    kk:  simbólica, opcional
        Residuo de la función en :math:`\\sigma_i` o :math:`\\omega \\to -\\omega_i`. El valor predeterminado es 0.
    ZRC_foster: simbólica
        Función inmitancia :math:`I(s)` a expresar como :math:`I_F(s)`

    Returns
    -------
    k0:  simbólica, opcional
        No está permitido para esta forma el residuo en 0.
    koo:  simbólica, opcional
        Residuo de la función en infinito o :math:`s \\to \\infty`.
    ki:  simbólica, list o tuple opcional
        Residuo de la función en :math:`\\omega_i` o :math:`s \\to -\\sigma_i`.
    kk:  simbólica, opcional
        Residuo de la función en :math:`\\sigma = 0`.
    YRC_foster: simbólica
        Función YRC expresada como :math:`I_F(s) = I(s)*s`    
    
    Raises
    ------
    ValueError
        Si cualquiera de los argumentos no son una instancia de sympy.Expr.


    See Also
    --------
    :func:`foster`
    :func:`foster_zRC2yRC`
    :func:`dibujar_foster_serie`

    
    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.sintesis_dipolo import foster, foster_zRC2yRC
    >>> from pytc2.dibujar import dibujar_foster_derivacion
    >>> s = sp.symbols('s ', complex=True)
    >>> # Sea la siguiente función de excitación
    >>> YRC = 2*(s**2 + 4*s + 3)/(s**2 + 8*s + 12)
    >>> k0, koo, ki_wi, kk, YRC_foster = foster(YRC/s)
    >>> k0, koo, ki_wi, kk, YRC_foster = foster_zRC2yRC(k0, koo, ki_wi, kk, YRC_foster)
    >>> dibujar_foster_derivacion(k0 = k0, koo = koo, ki = ki_wi, y_exc = YRC_foster)


    '''    
    
    if koo.is_zero:
    # koo tiene que ser 0 para ZRC ya que en inf habrá
    # o 0 o cte.
        
        if not(kk.is_zero):
            koo = kk
            kk = sp.Rational(0)
            
        if not(k0.is_zero):
            kk = k0
            k0 = sp.Rational(0)
            
        if (isinstance(ki_wi, sp.Expr) and not (ki_wi.is_zero)) or isinstance(ki_wi, list):
            
            ki = ki_wi
            # ki = []
            # for this_ki_wi in ki_wi:
                
            #     ki += [[this_ki_wi[1], this_ki_wi[0]]]
            
            
        YRC_foster = sp.expand(ZRC_foster * s)

    return([k0, koo, ki, kk, YRC_foster])

def foster( imm ):
    '''
    Expande una función inmitancia :math:`I(s)` en fracciones simples, de acuerdo al método 
    de Foster. La forma matemática es:

    .. math:: I(s)= \\frac{k_0}{s} + k_\\infty.s + \\sum_{i=1}^N\\frac{2.k_i.s}{s^2+\\omega_i^2}  

    Dependiendo la naturaleza de :math:`I(s)` como impedancia o admitancia, 
    resultará en los métodos de Foster serie, o paralelo. También existen 3 
    variantes 1) en caso que se trate de redes no disipativas (LC), y redes 
    disipativas compuestos solo por dos elementos circuitales: RC - RL. 2) Las 
    expresiones matemáticas para :math:`Z_{RC}` son las mismas que :math:`Y_{RL}`,
    mientras que 3) las de :math:`Z_{RL}` iguales a las de :math:`Y_{RC}`.


    Parameters
    ----------
    k0:  simbólica, opcional
        Residuo de la función en DC o :math:`s \\to 0`. El valor predeterminado es 0.
    koo:  simbólica, opcional
        Residuo de la función en infinito o :math:`s \\to \\infty`. El valor predeterminado es 0.
    ki:  simbólica, list o tuple opcional
        Residuo de la función en :math:`\\omega_i` o :math:`s^2 \\to -\\omega^2_i`. El valor predeterminado es 0.
    kk:  simbólica, opcional
        Residuo de la función en :math:`\\sigma_i` o :math:`\\omega \\to -\\omega_i`. El valor predeterminado es 0.
    

    Returns
    -------
    k0: simbólica, opcional
        El residuo en 0, expresado matemáticamente como :math:`\\frac{k_0}{s}`.
    koo: simbólica, opcional
        Residuo de la función en infinito o :math:`s \\to \\infty`, que se 
        corresponde al término :math:`k_\\infty*s`.
    ki: simbólica, list o tuple opcional
        Residuo de la función en :math:`s^2 \\to -\\omega_i^2`, matemáticamente
        :math:`\\frac{2.k_i.s}{s^2+\\omega_i^2}`.
    kk: simbólica, opcional
        Residuo de la función en :math:`\\sigma = 0`, para funciones disipativas.
    foster_form: simbólica
        Función YRC expresada como :math:`I_F(s) = I(s)*s`    


    Raises
    ------
    ValueError
        Si cualquiera de los argumentos no son una instancia de sympy.Expr.


    See Also
    --------
    :func:`foster`
    :func:`foster_zRC2yRC`
    :func:`dibujar_foster_paralelo`

    
    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.sintesis_dipolo import foster
    >>> from pytc2.dibujar import dibujar_foster_serie
    >>> s = sp.symbols('s ', complex=True)
    >>> # Sea la siguiente función de excitación
    >>> FF = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    >>> # Se expande FF a la Foster
    >>> k0, koo, ki_wi, _, FF_foster = foster(FF)
    >>> # Tratamos a nuestra función imitancia como una Z
    >>> dibujar_foster_serie(k0 = k0, koo = koo, ki = ki_wi, z_exc = FF)


    '''    

    if not isinstance(imm , sp.Expr):
        raise ValueError('Hay que definir imm como una expresión simbólica.')

    num, den = imm.as_numer_denom()
    
    # grados de P y Q
    # deg_P = sp.degree(num)
    # deg_Q = sp.degree(den)
    
    imm_foster = sp.polys.partfrac.apart(imm)
    
    all_terms = imm_foster.as_ordered_terms()
    
    kk = sp.Rational(0)
    k0 = sp.Rational(0)
    koo = sp.Rational(0)
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
    
            koo_i = sp.Rational(0)
            k0_i = sp.Rational(0)
            
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
        ki = sp.Rational(0)

    if not (sp.simplify(sp.expand(foster_form - imm))).is_zero:
        # error
        print_console_alert('Fallo la expansión')
        print_latex(expr_simb_expr(imm, foster_form, ' \\neq '))

    return([k0, koo, ki, kk, foster_form])


