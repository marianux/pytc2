#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:24:17 2023

@author: mariano
"""

import numpy as np

import sympy as sp

from IPython.display import display

from schemdraw import Drawing
from schemdraw.elements import  Resistor, ResistorIEC, Capacitor, Inductor, Line, Dot, Gap, Arrow


##########################################
#%% Variables para el análisis simbólico #
##########################################

from .general import s, to_latex, str_to_latex



########################################
#%% Funciones para dibujar cuadripolos #
########################################

def dibujar_Tee(ZZ, return_components=False):
    '''
    Dibuja una red Tee a partir de la matriz Z.

    Parameters
    ----------
    ZZ : sympy.Matrix
        Matriz de impedancia Z.
    return_components : bool, optional
        Indica si se deben devolver los componentes individuales de la red (Za, Zb, Zc). Por defecto es False.

    Returns
    -------
    list or None
        Si return_components es True, devuelve una lista con los componentes individuales de la red (Za, Zb, Zc). 
        Si return_components es False, no devuelve nada.

    Raises
    ------
    TypeError
        Si ZZ no es una instancia de sympy.Matrix.

    Examples
    --------
    >>> dibujar_Tee(sp.Matrix([[1, 2], [3, 4]]))
    # ... dibujo de la red Tee ...

    '''

    if not isinstance(ZZ, sp.Matrix):
        raise TypeError("ZZ debe ser una instancia de sympy.Matrix.")

    # Dibujo la red Tee
    d = Drawing(unit=4)

    d = dibujar_puerto_entrada(d, port_name='')

    Za = ZZ[0, 0] - ZZ[0, 1]
    Zb = ZZ[0, 1]
    Zc = ZZ[1, 1] - ZZ[0, 1]
    
    if isinstance(ZZ, sp.Basic):
        Za = sp.simplify(sp.expand(Za))
        Zb = sp.simplify(sp.expand(Zb))
        Zc = sp.simplify(sp.expand(Zc))

    if not Za.is_zero:
        d = dibujar_elemento_serie(d, ResistorIEC, Za)

    if not Zb.is_zero:
        d = dibujar_elemento_derivacion(d, ResistorIEC, Zb)

    if not Zc.is_zero:
        d = dibujar_elemento_serie(d, ResistorIEC, Zc)

    d = dibujar_puerto_salida(d, port_name='')

    display(d)

    if return_components:
        return [Za, Zb, Zc]
    
    
def dibujar_Pi(YY, return_components=False):
    '''
    Dibuja una red Pi a partir de la matriz Y.

    Breve descripción de la función

    Parameters
    ----------
    YY : Symbolic Matrix
        Matriz de admitancia Y.
    return_components : bool, optional
        Indica si se deben devolver los componentes individuales de la red (Ya, Yb, Yc). Por defecto es False.

    Returns
    -------
    None or list
        Si return_components es True, devuelve una lista con los componentes individuales de la red (Ya, Yb, Yc). 
        Si return_components es False, no devuelve nada.

    Raises
    ------
    None

    See Also
    --------
    :func:``

    Examples
    --------
    >>> dibujar_Pi(sp.Matrix([[2, -1], [-1, 3]]))
    [dibujo del red]
    >>> Ya, Yb, Yc = dibujar_Pi(sp.Matrix([[1, 0], [0, 1]]), return_components=True)
    [dibujo del red]

    '''

    # Comprobar el tipo de dato de YY
    if not isinstance(YY, sp.Matrix):
        raise TypeError("YY debe ser una matriz simbólica.")

    # Dibujo la red Pi
    d = Drawing(unit=4)

    d = dibujar_puerto_entrada(d, port_name='')

    Ya = YY[0, 0] + YY[0, 1]
    Yb = -YY[0, 1]
    Yc = YY[1, 1] + YY[0, 1]

    bSymbolic = isinstance(YY[0, 0], sp.Basic)

    if bSymbolic:
        Za = sp.simplify(sp.expand(1/Ya))
        Zb = sp.simplify(sp.expand(1/Yb))
        Zc = sp.simplify(sp.expand(1/Yc))
    else:
        Za = 1/Ya
        Zb = 1/Yb
        Zc = 1/Yc

    if (bSymbolic and (not Ya.is_zero) or (not bSymbolic) and Ya != 0):
        d = dibujar_elemento_derivacion(d, ResistorIEC, Za)

    if (bSymbolic and (not Yb.is_zero) or (not bSymbolic) and Yb != 0):
        d = dibujar_elemento_serie(d, ResistorIEC, Zb)

    if (bSymbolic and (not Yc.is_zero) or (not bSymbolic) and Yc != 0):
        d = dibujar_elemento_derivacion(d, ResistorIEC, Zc)

    d = dibujar_puerto_salida(d, port_name='')

    display(d)

    if return_components:
        return [Ya, Yb, Yc]
    
def dibujar_lattice(ZZ, return_components=False):
    '''
    Dibuja una red Lattice a partir de una matriz de parámetros Z.
    
    Parameters
    ----------
    ZZ : Matriz simbólica, opcional
        Parámetros Z de la red. Si no se proporciona, solo se genera el dibujo. El valor predeterminado es None.
    return_components : bool, opcional
        Indica si se deben devolver los componentes de la red Lattice simétrica (Za y Zb). El valor predeterminado es False.

    Returns
    -------
    list or None
        Si return_components es True, devuelve una lista con los componentes Za y Zb de la red Lattice simétrica.
        Si return_components es False, devuelve None.

    Raises
    ------
    ValueError
        Si ZZ no es una instancia de sympy.Matrix.
        Si ZZ no es de 2x2

    See Also
    --------
    :func:``

    Ejemplos
    --------
    >>> dibujar_lattice()
    >>> dibujar_lattice(sp.Matrix([[1, 0], [0, 1]]), return_components=True)
    '''

    if not isinstance(ZZ, (sp.Matrix, np.ndarray)):
        raise ValueError("ZZ debe ser una instancia de Symbolic Matrix")
    
    # Verificar que Spar tenga el formato correcto
    if ZZ.shape != (2, 2):
        raise ValueError("ZZ debe tener el formato [ [Z11, Z12], [Z21, Z22] ]")

    if ZZ is None:
        # Sin valores, solo el dibujo
        Za_lbl = 'Za'
        Zb_lbl = 'Zb'
        Za = 1
        Zb = 1
        bSymbolic = False
    else:
        # Calculo los valores de la matriz Z
        # z11 - z12
        Za = ZZ[0, 0] - ZZ[0, 1]
        Zb = ZZ[0, 0] + ZZ[0, 1]
        bSymbolic = isinstance(ZZ[0, 0], sp.Basic)
        
        if bSymbolic:
            Za = sp.simplify(Za)
            Zb = sp.simplify(Zb)
            Za_lbl = to_latex(Za)
            Zb_lbl = to_latex(Zb)
        else:
            Za_lbl = str_to_latex('{:3.3f}'.format(Za))
            Zb_lbl = str_to_latex('{:3.3f}'.format(Zb))

    # Dibujo la red Lattice
    with Drawing() as d:
        d.config(fontsize=16, unit=4)
        d = dibujar_puerto_entrada(d, port_name='')

        if (bSymbolic and (not Za.is_zero) or (not bSymbolic) and Za != 0):
            d += (Za_d := ResistorIEC().right().label(Za_lbl).dot().idot())
        else:
            d += (Za_d := Line().right().dot())

        d.push()
        d += Gap().down().label('')
        d += (line_down := Line(ls='dotted').left().dot().idot())
        cross_line_vec = line_down.end - Za_d.end
        d += Line(ls='dotted').endpoints(Za_d.end, Za_d.end + 0.25*cross_line_vec)
        d += Line(ls='dotted').endpoints(Za_d.end + 0.6*cross_line_vec, line_down.end)

        if (bSymbolic and (not Zb.is_zero) or (not bSymbolic) and Zb != 0):
            d += (Zb_d := ResistorIEC().label(Zb_lbl).endpoints(Za_d.start, line_down.start).dot())
        else:
            d += (Zb_d := Line().endpoints(Za_d.start, line_down.start).dot())

        d.pop()
        d = dibujar_puerto_salida(d, port_name='')

    if return_components:
        return [Za, Zb]
    
def dibujar_cauer_RC_RL(ki = None, y_exc = None, z_exc = None):
    '''
    Dibuja una red disipativa escalera (RC-RL) a partir de una expansión en 
    fracciones continuas (Método de Cauer). Dependiendo se especifique `z_exc`
    o `y_exc` y el tipo de residuos de `ki` se dibujará la red correspondiente.
    En caso que se trate de redes RC, la forma matemática será:

    .. math:: Z_{RC}(s)= \\frac{1}{s.C_1} + \\frac{1}{ \\frac{1}{R_1} + \\frac{1}{ \\frac{1}{s.C_2} + \\cdots } } = 
         R_1 + \\frac{1}{ s.C_1 + \\frac{1}{ R_2 + \\cdots } } 

    .. math:: Y_{RC}(s)= s.C_1 + \\frac{1}{ R_1 + \\frac{1}{ s.C_2 + \\cdots } } = 
         \\frac{1}{R_1} + \\frac{1}{ s.C_1 + \\frac{1}{ \\frac{1}{R_2} + \\cdots } } 

    Parameters
    ----------
    ki : lista con expresiones simbólicas
        Será una lista que contenga los residuos [k0, ki, koo ] como expresiones 
        simbólicas. Esta lista la provee la función :func:`cauer_RC`. El valor 
        predeterminado es None. Siendo:

        * k0  : Residuo de la función en DC o :math:`\\sigma \\to 0`.
        * koo : Residuo de la función en infinito o :math:`\\sigma \\to \\infty`.
        * ki  : Residuo de la función en :math:`\\sigma_i` o :math:`\\sigma \\to -\\sigma_i`

    Returns
    -------
    None

    Raises
    ------
    ValueError
        Si y_exc y z_exc no son una instancia de sympy.Expr.

    See Also
    --------
    :func:`cauer_RC`
    :func:`foster_zRC2yRC`
    :func:`dibujar_cauer_LC`

    Examples
    --------
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
    if not ( isinstance(y_exc , sp.Expr) and isinstance(z_exc , sp.Expr)):
        raise ValueError("'Hay que definir la función de excitación y_exc o z_exc como una expresión simbólica.'")

    if not(ki is None) or len(ki) > 0:
        # si hay algo para dibujar ...
        
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads

        d = dibujar_puerto_entrada(d,
                                       voltage_lbl = ('+', '$V$', '-'), 
                                       current_lbl = '$I$')

        if y_exc is None:
            
            bIsImpedance = True
            
            d, _ = dibujar_funcion_exc_abajo(d, 
                                                      'Z',  
                                                      z_exc, 
                                                      hacia_salida = True,
                                                      k_gap_width = 0.5)
        else:
            bIsImpedance = False
            
            d, _ = dibujar_funcion_exc_abajo(d, 
                                                      'Y',  
                                                      y_exc, 
                                                      hacia_salida = True,
                                                      k_gap_width = 0.5)
    
        if bIsImpedance:
            bSeries = True
        else:
            bSeries = False
        
        bComponenteDibujadoDerivacion = False

        for kii in ki:


            if bSeries:
                
                if sp.degree(kii*s) == 1:
                    d = dibujar_elemento_serie(d, Resistor, kii)
                elif sp.degree(kii*s) == 0:
                    d = dibujar_elemento_serie(d, Capacitor, 1/(s*kii))
                else:
                    d = dibujar_elemento_serie(d, Inductor, kii/s)
                    
                bComponenteDibujadoDerivacion = False

            else:

                if bComponenteDibujadoDerivacion:
                    
                    dibujar_espacio_derivacion(d)

                if sp.degree(kii*s) == 1:
                    d = dibujar_elemento_derivacion(d, Resistor, 1/kii)
                elif sp.degree(kii*s) == 2:
                    d = dibujar_elemento_derivacion(d, Capacitor, kii/s)
                else:
                    d = dibujar_elemento_derivacion(d, Inductor, 1/(s*kii))
                
                bComponenteDibujadoDerivacion = True

            bSeries = not bSeries

        if not bComponenteDibujadoDerivacion:
            
            d += Line().right().length(d.unit*.25)
            d += Line().down()
            d += Line().left().length(d.unit*.25)
        
        display(d)

    else:    
        
        print('Nada para dibujar')

def dibujar_cauer_LC(ki = None, y_exc = None, z_exc = None):
    '''
    Description
    -----------
    Dibuja una red disipativa escalera no disipativa a partir de una expansión en 
    fracciones continuas (Método de Cauer). Dependiendo se especifique `z_exc`
    o `y_exc` y el tipo de residuos de `ki` se dibujará la red correspondiente.
    La forma matemática será:

    .. math:: Z(s)= \\frac{1}{s.C_1} + \\frac{1}{ \\frac{1}{s.L_1} + \\frac{1}{ \\frac{1}{s.C_2} + \\cdots } } = 
             s.L_1 + \\frac{1}{ s.C_1 + \\frac{1}{ s.L_2 + \\cdots } } 

    .. math:: Y(s)= \\frac{1}{s.L_1} + \\frac{1}{ \\frac{1}{s.C_1} + \\frac{1}{ \\frac{1}{s.L_2} + \\cdots } } = 
             s.C_1 + \\frac{1}{ s.L_1 + \\frac{1}{ s.C_2 + \\cdots } }  

    Parameters
    ----------
    ki : lista con expresiones simbólicas
        Será una lista que contenga los residuos [k0, ki, koo ] como expresiones 
        simbólicas. Esta lista la provee la función :func:`cauer`. 
        El valor predeterminado es None. Siendo:

        * k0  : Residuo de la función en DC o :math:`s \\to 0`.
        * koo : Residuo de la función en infinito o :math:`s \\to \\infty`.
        * ki  : Residuo de la función en :math:`\\omega_i` o :math:`s^2 \\to -\\omega^2_i`
    

    Returns
    -------
    None

    Raises
    ------
    ValueError
        Si y_exc y z_exc no son una instancia de sympy.Expr.

    See Also
    --------
    :func:`cauer_RC`
    :func:`foster_zRC2yRC`
    :func:`dibujar_cauer_LC`

    
    Examples
    --------
    >>> s = sp.symbols('s ', complex=True)
    >>> # Sea la siguiente función de excitación
    >>> FF = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    >>> # Implementaremos FF mediante Cauer 1 o remociones continuas en infinito
    >>> koo, F_cauer_oo, rem = cauer_LC(FF, remover_en_inf=True)
    >>> print_latex(a_equal_b_latex_s(a_equal_b_latex_s('F(s)', FF)[1:-1], F_cauer_oo ))
    >>> # Tratamos a nuestra función inmitancia como una Z
    >>> dibujar_cauer_LC(koo, z_exc = F_cauer_oo)
    >>> # Tratamos a nuestra función inmitancia como una Y
    >>> dibujar_cauer_LC(koo, y_exc = F_cauer_oo)
    '''    
    if not ( isinstance(y_exc , sp.Expr) and isinstance(z_exc , sp.Expr)):
        raise ValueError("'Hay que definir la función de excitación y_exc o z_exc como una expresión simbólica.'")

    
    if y_exc is None and z_exc is None:

        assert('Hay que definir si se trata de una impedancia o admitancia')

    if not(ki is None) or len(ki) > 0:
        # si hay algo para dibujar ...
        
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads

        d = dibujar_puerto_entrada(d,
                                       voltage_lbl = ('+', '$V$', '-'), 
                                       current_lbl = '$I$')

        if y_exc is None:
            
            bIsImpedance = True
            
            d, _ = dibujar_funcion_exc_abajo(d, 
                                                      'Z',  
                                                      z_exc, 
                                                      hacia_salida = True,
                                                      k_gap_width = 0.5)
        else:
            bIsImpedance = False
            
            d, _ = dibujar_funcion_exc_abajo(d, 
                                                      'Y',  
                                                      y_exc, 
                                                      hacia_salida = True,
                                                      k_gap_width = 0.5)
    
        if bIsImpedance:
            bSeries = True
        else:
            bSeries = False

        # 1/s me da orden 1, atenti.
        if sp.degree(ki[0]*s) == 2 :
            bCauer1 = True
        else:
            bCauer1 = False
        
        
        bComponenteDibujadoDerivacion = False

        for kii in ki:


            if bSeries:
                
                if bCauer1:
                    d = dibujar_elemento_serie(d, Inductor, kii/s)
                else:
                    d = dibujar_elemento_serie(d, Capacitor, 1/(s*kii))
                    
                bComponenteDibujadoDerivacion = False

            else:

                if bComponenteDibujadoDerivacion:
                    
                    dibujar_espacio_derivacion(d)

                if bCauer1:
                    d = dibujar_elemento_derivacion(d, Capacitor, kii/s)
                else:
                    d = dibujar_elemento_derivacion(d, Inductor, 1/(s*kii))
                
                bComponenteDibujadoDerivacion = True

            bSeries = not bSeries

        if not bComponenteDibujadoDerivacion:
            
            d += Line().right().length(d.unit*.25)
            d += Line().down()
            d += Line().left().length(d.unit*.25)
        
        display(d)

    else:    
        
        print('Nada para dibujar')

def dibujar_foster_derivacion(k0 = None, koo = None, ki = None, kk = None, y_exc = None):
    '''
    Description
    -----------
    Draws a parallel non-disipative admitance following Foster synthesis method.

        Y = k0 / s + koo * s +  1 / ( k0_i / s + koo_i * s ) 
    
    Parameters
    ----------
    k0 : symbolic positive real number. The residue value at DC or s->0.
        
    koo : symbolic positive real number. The residue value at inf or s->oo.
        
    ki : symbolic positive real array of numbers. A list of residue pairs at 
         each i-th finite pole or s**2->-(w_i**2). The first element of the pair
         is the k0_i value (capacitor), while the other is the koo_i (inductor)
         value.

    Returns
    -------
    The drawing object.
    
    Ejemplo
    -------

    # Sea la siguiente función de excitación
    Imm = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    
    # Implementaremos Imm mediante Foster
    k0, koo, ki = tc2.foster(Imm)
    
    # Tratamos a nuestra función imitancia como una Z
    tc2.dibujar_foster_derivacion(k0, koo, ki, y_exc = Imm)

    '''    

    if not(k0 is None and koo is None and ki is None and kk is None):
        
        if kk is None:
            bDisipativo = False
        else:
            bDisipativo = True
        
        # si hay algo para dibujar ...
        
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads

        bComponenteDibujado = False

        d = dibujar_puerto_entrada(d,
                                       voltage_lbl = ('+', '$V$', '-'), 
                                       current_lbl = '$I$')

        if not(y_exc is None):
            d, _ = dibujar_funcion_exc_abajo(d, 
                                                      'Y',  
                                                      y_exc, 
                                                      hacia_salida = True,
                                                      k_gap_width = 0.5)


        if not(kk is None):
            
            d = dibujar_elemento_derivacion(d, Resistor, 1/kk)

            bComponenteDibujado = True

        if not(k0 is None):

            if bComponenteDibujado:
                
                dibujar_espacio_derivacion(d)

            d = dibujar_elemento_derivacion(d, Inductor, 1/k0)
            
            bComponenteDibujado = True
            
            
        if not(koo is None):
        
            if bComponenteDibujado:
                
                dibujar_espacio_derivacion(d)
                    
            d = dibujar_elemento_derivacion(d, Capacitor, koo)

            bComponenteDibujado = True
            
        if not(ki is None):

            for un_tanque in ki:

                if bComponenteDibujado:
                    
                    dibujar_espacio_derivacion(d)
                
                if bDisipativo:
                    
                    if k0 is None:
                        d = dibujar_tanque_RC_derivacion(d, capacitor_lbl = 1/un_tanque[0], sym_R_label = un_tanque[1] )
                        bComponenteDibujado = True
                    else:
                        d = dibujar_tanque_RL_derivacion(d, sym_ind_label = un_tanque[1], sym_R_label = un_tanque[0] )
                        bComponenteDibujado = True
                        
                else:    
                
                    d = dibujar_tanque_derivacion(d, inductor_lbl = un_tanque[1], capacitor_lbl = 1/un_tanque[0])
                    bComponenteDibujado = True

        
        display(d)

    else:    
        
        print('Nada para dibujar')

def dibujar_foster_serie(k0 = None, koo = None, ki = None, kk = None, z_exc = None):
    '''
    Description
    -----------
    Draws a series non-disipative impedance following Foster synthesis method.

        Z = k0 / s + koo * s +  1 / ( k0_i / s + koo_i * s ) 
    
    Parameters
    ----------
    k0 : symbolic positive real number. The residue value at DC or s->0.
        
    koo : symbolic positive real number. The residue value at inf or s->oo.
        
    ki : symbolic positive real array of numbers. A list of residue pairs at 
         each i-th finite pole or s**2->-(w_i**2). The first element of the pair
         is the k0_i value (inductor), while the other is the koo_i (capacitor)
         value.

    Returns
    -------
    The drawing object.
    
    Ejemplo
    -------

    # Sea la siguiente función de excitación
    Imm = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    
    # Implementaremos Imm mediante Foster
    k0, koo, ki = tc2.foster(Imm)
    
    # Tratamos a nuestra función imitancia como una Z
    tc2.dibujar_foster_serie(k0, koo, ki, z_exc = Imm)

    '''    

    if not(k0 is None and koo is None and ki is None and kk is None):
        
        
        if kk is None:
            bDisipativo = False
        else:
            bDisipativo = True
        
        # si hay algo para dibujar ...
        
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads

        d = dibujar_puerto_entrada(d,
                                       voltage_lbl = ('+', '$V$', '-'), 
                                       current_lbl = '$I$')

        if not(z_exc is None):
            d, z5_lbl = dibujar_funcion_exc_abajo(d, 
                                                      'Z',  
                                                      z_exc, 
                                                      hacia_salida = True,
                                                      k_gap_width = 0.5)

        if not(kk is None):
            
            d = dibujar_elemento_serie(d, Resistor, kk)
            
        if not(k0 is None):
        
            d = dibujar_elemento_serie(d, Capacitor, 1/k0)
            
        if not(koo is None):
        
            d = dibujar_elemento_serie(d, Inductor, koo)
            
        if not(ki is None):

            for un_tanque in ki:
                
                if bDisipativo:
                    
                    if k0 is None:
                        d = dibujar_tanque_RL_serie(d, sym_ind_label = 1/un_tanque[0], sym_R_label = 1/un_tanque[1] )
                    else:
                        d = dibujar_tanque_RC_serie(d, sym_R_label = 1/un_tanque[0], capacitor_lbl = un_tanque[1] )
                        
                else:    
                    d = dibujar_tanque_serie(d, sym_ind_label = 1/un_tanque[0], sym_cap_label = un_tanque[1] )

                dibujar_espacio_derivacion(d)


        d += Line().right().length(d.unit*.25)
        d += Line().down()
        d += Line().left().length(d.unit*.25)
        
        display(d)

    else:    
        
        print('Nada para dibujar')

##################################################
#%% Funciones para dibujar redes de forma bonita #
##################################################

def dibujar_puerto_entrada(d, port_name = None, voltage_lbl = None, current_lbl = None):
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
    
    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    d += Dot(open=True)
    
    if voltage_lbl is None:
        d += Gap().down().label( '' )
    else:
        d += Gap().down().label( voltage_lbl, fontsize=16)
    
    d.push()

    if not(port_name is None):
        d += Gap().left().label( '' ).length(d.unit*.35)
        d += Gap().up().label( port_name, fontsize=22)
        d.pop()
        
    d += Dot(open=True)
    d += Line().right().length(d.unit*.5)
    d += Gap().up().label( '' )
    d.push()
    
    if current_lbl is None:
        d += Line().left().length(d.unit*.5)
    else:
        d += Line().left().length(d.unit*.25)
        d += Arrow(reverse=True).left().label( current_lbl, fontsize=16).length(d.unit*.25)
    
    d.pop()

    return(d)

def dibujar_puerto_salida(d, port_name = None, voltage_lbl = None, current_lbl = None):
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
    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    if current_lbl is None:
        d += Line().right().length(d.unit*.5)
    else:
        d += Line().right().length(d.unit*.25)
        d += Arrow(reverse=True).right().label( current_lbl, fontsize=16).length(d.unit*.25)
    
    d += Dot(open=True)
    
    d.push()

    if voltage_lbl is None:
        d += Gap().down().label( '' )
    else:
        d += Gap().down().label( voltage_lbl, fontsize=16)


    if not(port_name is None):
        d.push()
        d += Gap().right().label( '' ).length(d.unit*.35)
        d += Gap().up().label( port_name, fontsize=22)
        d.pop()

    d += Dot(open=True)
    d += Line().left().length(d.unit*.5)

    d.pop()

    return(d)

def dibujar_espaciador( d ):
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

    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads

    d += Line().right().length(d.unit*.5)

    d.push()

    d += Gap().down().label( '' )

    d += Line().left().length(d.unit*.5)

    d.pop()

    return(d)

def dibujar_funcion_exc_abajo(d, func_label, sym_func, k_gap_width=0.5, hacia_salida  = False, hacia_entrada  = False ):
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

    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads

    half_width = d.unit*k_gap_width/2
    
    d += Line().right().length(half_width)
    d.push()
    d += Gap().down().label('')
    d.push()
    
    if isinstance(sym_func, sp.Basic ):
        sym_func = '$ ' + func_label + ' = ' + sp.latex(sym_func) + ' $'
    elif isinstance(sym_func, np.number):
        sym_func =  '$ ' + func_label + ' = ' + '{:3.3f}'.format(sym_func) + ' $'
    elif isinstance(sym_func, str):
        sym_func = '$ ' + func_label + ' = ' +  sym_func + ' $'
    else:
        sym_func = '$ ' + func_label + ' = ?? $'
    
    lbl = d.add(Gap().down().label( sym_func, fontsize=22 ).length(0.5*half_width))
    d += Gap().down().label('').length(0.5*half_width)
    d.pop()
    d.push()
    d += Line().up().at( (d.here.x, d.here.y - .2 * half_width) ).length(half_width).linewidth(1)
    
    if( hacia_salida ):
        d.push()
        d += Arrow().right().length(.5*half_width).linewidth(1)
        d.pop()
        
    if( hacia_entrada ):
        d += Arrow().left().length(.5*half_width).linewidth(1)
        
    d.pop()
    d.push()
    d += Line().left().length(half_width)
    d.pop()
    d += Line().right().length(half_width)
    d.pop()
    d += Line().right().length(half_width)

    return([d, lbl])

def dibujar_funcion_exc_arriba(d, func_label, sym_func, k_gap_width=0.5, hacia_salida = False, hacia_entrada = False ):
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
    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads

    half_width = d.unit*k_gap_width/2
    
    d += Line().right().length(half_width)
    d.push()
    
    if isinstance(sym_func, sp.Basic ):
        sym_func = '$ ' + func_label + ' = ' + sp.latex(sym_func) + ' $'
    elif isinstance(sym_func, np.number):
        sym_func =  '$ ' + func_label + ' = ' + '{:3.3f}'.format(sym_func) + ' $'
    elif isinstance(sym_func, str):
        sym_func = '$ ' + func_label + ' = ' +  sym_func + ' $'
    else:
        sym_func = '$ ' + func_label + ' = ?? $'

    
    lbl = d.add(Gap().up().label( sym_func, fontsize=22 ).length(3* half_width))
    d.pop()
    d.push()
    d += Line().down().at( (d.here.x, d.here.y + .2 * half_width) ).length(half_width).linewidth(1)
    
    if( hacia_salida ):
        d.push()
        d += Arrow().right().length(.5*half_width).linewidth(1)
        d.pop()
        
    if( hacia_entrada ):
        d += Arrow().left().length(.5*half_width).linewidth(1)
        
    d.pop()
    d.push()
    d += Gap().down().label('')
    d.push()
    d += Line().left().length(half_width)
    d.pop()
    d += Line().right().length(half_width)
    d.pop()
    d += Line().right().length(half_width)



    return([d, lbl])

def dibujar_elemento_serie(d, elemento, sym_label=''):
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

    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    if isinstance(sym_label, sp.Basic ):
        sym_label = to_latex(sym_label)
    elif isinstance(sym_label, np.number):
        sym_label = str_to_latex('{:3.3f}'.format(sym_label))
    elif isinstance(sym_label, str):
        sym_label = str_to_latex(sym_label)
    else:
        sym_label = '$ ?? $'

    
    d += elemento().right().label(sym_label, fontsize=16)
    d.push()
    d += Gap().down().label( '' )
    d += Line().left()
    d.pop()

    return(d)

def dibujar_espacio_derivacion(d):
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
    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads

    d += Line().right().length(d.unit*.5)
    d.push()
    d += Gap().down().label( '' )
    d += Line().left().length(d.unit*.5)
    d.pop()

    return(d)

        
def dibujar_cierre(d):
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
    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads

    d += Line().right().length(d.unit*.5)
    d.push()
    d += Line().down()
    d += Line().left().length(d.unit*.5)
    d.pop()

    return(d)


def dibujar_elemento_derivacion(d, elemento, sym_label=''):
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
    
    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    if isinstance(sym_label, sp.Basic ):
        sym_label = to_latex(sym_label)
    elif isinstance(sym_label, np.number):
        sym_label = str_to_latex('{:3.3f}'.format(sym_label))
    elif isinstance(sym_label, str):
        sym_label = str_to_latex(sym_label)
    else:
        sym_label = '$ ?? $'
    
    d += Dot()
    d.push()
    d += elemento().down().label(sym_label, fontsize=16)
    d += Dot()
    d.pop()

    return(d)

def dibujar_tanque_RC_serie(d, sym_R_label='', capacitor_lbl=''):
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
    
    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    if isinstance(sym_R_label, sp.Basic ):
        sym_R_label = to_latex(sym_R_label)
    else:
        sym_R_label = str_to_latex(sym_R_label)
    
    if isinstance(capacitor_lbl, sp.Basic ):
        capacitor_lbl = to_latex(capacitor_lbl)
    else:
        capacitor_lbl = str_to_latex(capacitor_lbl)
    
    d.push()
    d += Dot()
    d += Capacitor().right().label(capacitor_lbl, fontsize=16)
    d.pop()
    d += Line().up().length(d.unit*.5)
    d += Resistor().right().label(sym_R_label, fontsize=16)
    d += Line().down().length(d.unit*.5)
    d += Dot()
    d.push()
    d += Gap().down().label( '' )
    d += Line().left()
    d.pop()

    return(d)

def dibujar_tanque_RC_derivacion(d, sym_R_label='', capacitor_lbl=''):
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

    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    if isinstance(sym_R_label, sp.Basic ):
        sym_R_label = to_latex(sym_R_label)
    else:
        sym_R_label = str_to_latex(sym_R_label)
    
    if isinstance(capacitor_lbl, sp.Basic ):
        capacitor_lbl = to_latex(capacitor_lbl)
    else:
        capacitor_lbl = str_to_latex(capacitor_lbl)
    
    d.push()
    d += Dot()
    d += Capacitor().down().label(capacitor_lbl, fontsize=16).length(d.unit*.5)
    d += Resistor().down().label(sym_R_label, fontsize=16).length(d.unit*.5)
    d += Dot()
    d.pop()

    return(d)

def dibujar_tanque_RL_serie(d, sym_R_label='', sym_ind_label=''):
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

    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    if isinstance(sym_R_label, sp.Basic ):
        sym_R_label = to_latex(sym_R_label)
    else:
        sym_R_label = str_to_latex(sym_R_label)
    
    if isinstance(sym_ind_label, sp.Basic ):
        sym_ind_label = to_latex(sym_ind_label)
    else:
        sym_ind_label = str_to_latex(sym_ind_label)
    
    d.push()
    d += Dot()
    d += Inductor().right().label(sym_ind_label, fontsize=16)
    d.pop()
    d += Line().up().length(d.unit*.5)
    d += Resistor().right().label(sym_R_label, fontsize=16)
    d += Line().down().length(d.unit*.5)
    d += Dot()
    d.push()
    d += Gap().down().label( '' )
    d += Line().left()
    d.pop()

    return(d)

def dibujar_tanque_RL_derivacion(d, sym_R_label='', sym_ind_label=''):
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

    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    if isinstance(sym_R_label, sp.Basic ):
        sym_R_label = to_latex(sym_R_label)
    else:
        sym_R_label = str_to_latex(sym_R_label)
    
    if isinstance(sym_ind_label, sp.Basic ):
        sym_ind_label = to_latex(sym_ind_label)
    else:
        sym_ind_label = str_to_latex(sym_ind_label)
    
    d.push()
    d += Dot()
    d += Inductor().down().label(sym_ind_label, fontsize=16).length(d.unit*.5)
    d += Resistor().down().label(sym_R_label, fontsize=16).length(d.unit*.5)
    d += Dot()
    d.pop()

    return(d)

def dibujar_tanque_serie(d, sym_ind_label='', sym_cap_label=''):
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

    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    if isinstance(sym_cap_label, sp.Basic ):
        sym_cap_label = to_latex(sym_cap_label)
    else:
        sym_cap_label = str_to_latex(sym_cap_label)
    
    if isinstance(sym_ind_label, sp.Basic ):
        sym_ind_label = to_latex(sym_ind_label)
    else:
        sym_ind_label = str_to_latex(sym_ind_label)
    
    d.push()
    d += Dot()
    d += Inductor().right().label(sym_ind_label, fontsize=16)
    d.pop()
    d += Line().up().length(d.unit*.5)
    d += Capacitor().right().label(sym_cap_label, fontsize=16)
    d += Line().down().length(d.unit*.5)
    d += Dot()
    d.push()
    d += Gap().down().label( '' )
    d += Line().left()
    d.pop()

    return(d)

def dibujar_tanque_derivacion(d, inductor_lbl='', capacitor_lbl=''):
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

    if not isinstance(d, Drawing):
        d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    if isinstance(inductor_lbl, sp.Basic ):
        inductor_lbl = to_latex(inductor_lbl)
    else:
        inductor_lbl = str_to_latex(inductor_lbl)
    
    if isinstance(capacitor_lbl, sp.Basic ):
        capacitor_lbl = to_latex(capacitor_lbl)
    else:
        capacitor_lbl = str_to_latex(capacitor_lbl)
    
    d.push()
    d += Dot()
    d += Capacitor().down().label(capacitor_lbl, fontsize=16).length(d.unit*.5)
    d += Inductor().down().label(inductor_lbl, fontsize=16).length(d.unit*.5)
    d += Dot()
    d.pop()

    return(d)

