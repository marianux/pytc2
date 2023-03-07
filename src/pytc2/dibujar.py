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

# Laplace complex variable. s = σ + j.ω
s = sp.symbols('s', complex=True)
# Fourier real variable ω 
w = sp.symbols('w', complex=False)



#########################
#%% Funciones generales #
#########################



def to_latex( unsimbolo ):
    '''
    Convierte un símbolo en un string formateado para visualizarse en LaTex 
    '''
    
    return('$'+ sp.latex(unsimbolo) + '$')

def str_to_latex( unstr):
    '''
    Formatea un string para visualizarse en LaTex 
    '''
    
    return('$'+ unstr + '$')


########################################
#%% Funciones para dibujar cuadripolos #
########################################

def dibujar_Tee(ZZ):
    '''
    Dibuja una red Tee a partir de la matriz Z.

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
    
    # Dibujo la red Tee
    
    d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    d = dibujar_puerto_entrada(d,
                                   port_name = 'In' )
    
    Za = ZZ[0,0] - ZZ[0,1] 
    Zb = ZZ[0,1] 
    Zc = ZZ[1,1] - ZZ[0,1] 
    
    d = dibujar_elemento_serie(d, ResistorIEC, Za )
    d = dibujar_elemento_derivacion(d, ResistorIEC, Zb )
    d = dibujar_elemento_serie(d, ResistorIEC, Zc )
    
    d = dibujar_puerto_salida(d, 
                                  port_name = 'Out')

    display(d)        
    
    return([Za,Zb,Zc])


def dibujar_Pi(YY):
    '''
    Dibuja una red Pi a partir de la matriz Y.

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
    
    # Dibujo la red Tee
    
    d = Drawing(unit=4)  # unit=2 makes elements have shorter than normal leads
    
    d = dibujar_puerto_entrada(d,
                                   port_name = 'In')
    
    Ya = YY[0,0] + YY[0,1]
    Yb = -YY[0,1]
    Yc = YY[1,1] + YY[0,1]
    
    if isinstance(YY[0,0], sp.Symbol):
        
        Za = sp.simplify(sp.expand(1/Ya))
        Zb = sp.simplify(sp.expand(1/Yb))
        Zc = sp.simplify(sp.expand(1/Yc))
        
    else:

        Za = 1/(YY[0,0] + YY[0,1])
        Zb = 1/(-YY[0,1])
        Zc = 1/(YY[1,1] + YY[0,1])
    
    d = dibujar_elemento_derivacion(d, ResistorIEC, Za )
    d = dibujar_elemento_serie(d, ResistorIEC, Zb )
    d = dibujar_elemento_derivacion(d, ResistorIEC, Zc )
    
    d = dibujar_puerto_salida(d, 
                                  port_name = 'Out')
    
    display(d)        
    
    return([Ya, Yb, Yc])



def dibujar_cauer_RC_RL(ki = None, y_exc = None, z_exc = None):
    '''
    Description
    -----------
    Draws a parallel non-disipative admitance following Foster synthesis method.

        YorZ = ki / s +  1 / ( ki_i / s + koo_i * s ) 
    
    Parameters
    ----------
    ki : symbolic positive real number. The residue value at DC or s->0.
        
    koo : symbolic positive real number. The residue value at inf or s->oo.
        
    ki : symbolic positive real array of numbers. A list of residue pairs at 
         each i-th finite pole or s**2->-(w_i**2). The first element of the pair
         is the ki_i value (capacitor), while the other is the koo_i (inductor)
         value.

    Returns
    -------
    The drawing object.
    
    Ejemplo
    -------

    # Sea la siguiente función de excitación
    Imm = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    
    # Implementaremos Imm mediante Foster
    ki, koo, ki = tc2.foster(Imm)
    
    # Tratamos a nuestra función imitancia como una Z
    tc2.dibujar_foster_derivacion(ki, koo, ki, y_exc = Imm)

    '''    
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
    Draws a parallel non-disipative admitance following Foster synthesis method.

        YorZ = ki / s +  1 / ( ki_i / s + koo_i * s ) 
    
    Parameters
    ----------
    ki : symbolic positive real number. The residue value at DC or s->0.
        
    koo : symbolic positive real number. The residue value at inf or s->oo.
        
    ki : symbolic positive real array of numbers. A list of residue pairs at 
         each i-th finite pole or s**2->-(w_i**2). The first element of the pair
         is the ki_i value (capacitor), while the other is the koo_i (inductor)
         value.

    Returns
    -------
    The drawing object.
    
    Ejemplo
    -------

    # Sea la siguiente función de excitación
    Imm = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
    
    # Implementaremos Imm mediante Foster
    ki, koo, ki = tc2.foster(Imm)
    
    # Tratamos a nuestra función imitancia como una Z
    tc2.dibujar_foster_derivacion(ki, koo, ki, y_exc = Imm)

    '''    
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



def dibujar_foster_derivacion(k0 = None, koo = None, ki = None, y_exc = None):
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

    if not(k0 is None and koo is None and ki is None):
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

        if not(k0 is None):
        
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
                
                d = dibujar_tanque_derivacion(d, inductor_lbl = un_tanque[1], capacitor_lbl = 1/un_tanque[0])

                bComponenteDibujado = True

        
        display(d)

    else:    
        
        print('Nada para dibujar')


def dibujar_foster_serie(k0 = None, koo = None, ki = None, z_exc = None):
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

    if not(k0 is None and koo is None and ki is None):
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

        if not(k0 is None):
        
            d = dibujar_elemento_serie(d, Capacitor, 1/k0)
            
        if not(koo is None):
        
            d = dibujar_elemento_serie(d, Inductor, koo)
            
        if not(ki is None):

            for un_tanque in ki:
                
                d = dibujar_tanque_serie(d, sym_ind_label = 1/un_tanque[0], sym_cap_label = un_tanque[1] )

                dibujar_espacio_derivacion(d)


        d += Line().right().length(d.unit*.25)
        d += Line().down()
        d += Line().left().length(d.unit*.25)
        
        display(d)
        
        return(d)

    else:    
        
        print('Nada para dibujar')





##################################################
#%% Funciones para dibujar redes de forma bonita #
##################################################


def dibujar_puerto_entrada(d, port_name = None, voltage_lbl = None, current_lbl = None):
    
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

    d += Line().right().length(d.unit*.5)

    d.push()

    d += Gap().down().label( '' )

    d += Line().left().length(d.unit*.5)

    d.pop()

    return(d)


def dibujar_funcion_exc_abajo(d, func_label, sym_func, k_gap_width=0.5, hacia_salida  = False, hacia_entrada  = False ):
# ok

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
    
    if isinstance(sym_label, sp.Number ):
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

    d += Line().right().length(d.unit*.25)
    d.push()
    d += Gap().down().label( '' )
    d += Line().left().length(d.unit*.25)
    d.pop()

    return(d)

def dibujar_elemento_derivacion(d, elemento, sym_label=''):
    
    if isinstance(sym_label, sp.Number ):
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
    
    if isinstance(sym_R_label, sp.Number ):
        sym_R_label = to_latex(sym_R_label)
    else:
        sym_R_label = str_to_latex(sym_R_label)
    
    if isinstance(capacitor_lbl, sp.Number ):
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
    
    if isinstance(sym_R_label, sp.Number ):
        sym_R_label = to_latex(sym_R_label)
    else:
        sym_R_label = str_to_latex(sym_R_label)
    
    if isinstance(capacitor_lbl, sp.Number ):
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
    
    if isinstance(sym_R_label, sp.Number ):
        sym_R_label = to_latex(sym_R_label)
    else:
        sym_R_label = str_to_latex(sym_R_label)
    
    if isinstance(sym_ind_label, sp.Number ):
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
    
    if isinstance(sym_R_label, sp.Number ):
        sym_R_label = to_latex(sym_R_label)
    else:
        sym_R_label = str_to_latex(sym_R_label)
    
    if isinstance(sym_ind_label, sp.Number ):
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
    
    if isinstance(sym_cap_label, sp.Number ):
        sym_cap_label = to_latex(sym_cap_label)
    else:
        sym_cap_label = str_to_latex(sym_cap_label)
    
    if isinstance(sym_ind_label, sp.Number ):
        sym_ind_label = to_latex(sym_ind_label)
    else:
        sym_ind_label = str_to_latex(sym_ind_label)
    
    d.push()
    d += Dot()
    d += Inductor().right().label(sym_ind_label, fontsize=16)
    d.pop()
    d += Line().up().length(d.unit*.5)
    d += Resistor().right().label(sym_cap_label, fontsize=16)
    d += Line().down().length(d.unit*.5)
    d += Dot()
    d.push()
    d += Gap().down().label( '' )
    d += Line().left()
    d.pop()

    return(d)

# def dibujar_tanque_RL_derivacion(d, sym_R_label='', inductor_lbl=''):
    
#     if isinstance(sym_R_label, sp.Number ):
#         sym_R_label = to_latex(sym_R_label)
#     else:
#         sym_R_label = str_to_latex(sym_R_label)
    
#     if isinstance(inductor_lbl, sp.Number ):
#         inductor_lbl = to_latex(inductor_lbl)
#     else:
#         inductor_lbl = str_to_latex(inductor_lbl)
    
#     d.push()
#     d += Dot()
#     d += Inductor().down().label(inductor_lbl, fontsize=16).length(d.unit*.5)
#     d += Resistor().down().label(sym_R_label, fontsize=16).length(d.unit*.5)
#     d += Dot()
#     d.pop()

#     return(d)

# def dibujar_tanque_serie(d, inductor_lbl='', capacitor_lbl=''):
    
#     if isinstance(inductor_lbl, sp.Number ):
#         inductor_lbl = to_latex(inductor_lbl)
#     else:
#         inductor_lbl = str_to_latex(inductor_lbl)
    
#     if isinstance(capacitor_lbl, sp.Number ):
#         capacitor_lbl = to_latex(capacitor_lbl)
#     else:
#         capacitor_lbl = str_to_latex(capacitor_lbl)
    
#     d.push()
#     d += Dot()
#     d += Capacitor().right().label(capacitor_lbl, fontsize=16)
#     d.pop()
#     d += Line().up().length(d.unit*.5)
#     d += Inductor().right().label(inductor_lbl, fontsize=16)
#     d += Line().down().length(d.unit*.5)
#     d += Dot()
#     d.push()
#     d += Gap().down().label( '' )
#     d += Line().left()
#     d.pop()

#     return(d)

def dibujar_tanque_derivacion(d, inductor_lbl='', capacitor_lbl=''):
    
    if isinstance(inductor_lbl, sp.Number ):
        inductor_lbl = to_latex(inductor_lbl)
    else:
        inductor_lbl = str_to_latex(inductor_lbl)
    
    if isinstance(capacitor_lbl, sp.Number ):
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
