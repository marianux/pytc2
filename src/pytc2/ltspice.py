#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:22:31 2023

@author: mariano
"""

import numpy as np

import sympy as sp


from os import path
import shutil

import time

from numbers import Real

############################################
#%% Variables para la interfaz con LTspice #
############################################

ltux = 16 
"""
Unidades X para dibujar en la hoja de LTspice
"""
ltuy = 16
"""
Unidades Y para dibujar en la hoja de LTspice
"""

filename_eq_base = 'ltspice_equalizador_base'
"""
Archivo marco contenedor de las redes sintetizadas como ecualizadores/filtros
"""

cap_num = 1
"""
cuenta de capacitores
"""
res_num = 1
"""
cuenta de resistores
"""
ind_num = 1
"""
cuenta de inductores
"""
node_num = 1
"""
cuenta de nodos
"""


cur_x = 0
"""
cursor X para la localización de componentes
"""
cur_y = 0
"""
cursor Y para la localización de componentes
"""

lt_wire_length = 4 # ltux/ltuy unidades normalizadas
"""
tamaño estandard del cable
"""

#####
# Palabras clave del LTspice para disponer los componentes en el
# esquemático.

# 

res_der_str = [ 'SYMBOL res {:d} {:d} R0\n', # posición absoluta X-Y en el esquemático
                'WINDOW 0 48 43 Left 2\n', # posiciones relativas de etiquetas
                'WINDOW 3 47 68 Left 2\n', # posiciones relativas de etiquetas
                'SYMATTR InstName {:s}\n', # etiqueta que tendrá
                'SYMATTR Value {:3.5f}\n' # valor que tendrá
               ]
"""
resistor en derivacion
"""


ind_der_str = [ 'SYMBOL ind {:d} {:d} R0\n', # posición absoluta X-Y en el esquemático
                'WINDOW 0 47 34 Left 2\n', # posiciones relativas de etiquetas
                'WINDOW 3 43 65 Left 2\n', # posiciones relativas de etiquetas
                'SYMATTR InstName {:s}\n', # etiqueta que tendrá
                'SYMATTR Value {:3.5f}\n' # valor que tendrá
               ]
"""
inductor en derivacion
"""

cap_der_str = [ 'SYMBOL cap {:d} {:d} R0\n', # posición absoluta X-Y en el esquemático
                'WINDOW 0 48 18 Left 2\n', # posiciones relativas de etiquetas
                'WINDOW 3 45 49 Left 2\n', # posiciones relativas de etiquetas
                'SYMATTR InstName {:s}\n', # etiqueta que tendrá
                'SYMATTR Value {:3.5f}\n' # valor que tendrá
               ]
"""
capacitor en derivacion
"""

res_ser_str = [ 'SYMBOL res {:d} {:d} R90\n', # posición absoluta X-Y en el esquemático
                'WINDOW 0 -7 86 VBottom 2\n', # posiciones relativas de etiquetas
                'WINDOW 3 -36 24 VTop 2\n', # posiciones relativas de etiquetas
                'SYMATTR InstName {:s}\n', # etiqueta que tendrá
                'SYMATTR Value {:3.5f}\n' # valor que tendrá
               ]
"""
resistor en serie
"""

ind_ser_str = [ 'SYMBOL ind {:d} {:d} R270\n', # posición absoluta X-Y en el esquemático
                'WINDOW 0 40 19 VTop 2\n', # posiciones relativas de etiquetas
                'WINDOW 3 67 100 VBottom 2\n', # posiciones relativas de etiquetas
                'SYMATTR InstName {:s}\n', # etiqueta que tendrá
                'SYMATTR Value {:3.5f}\n' # valor que tendrá
               ]
"""
inductor en serie
"""

cap_ser_str = [ 'SYMBOL cap {:d} {:d} R90\n', # posición absoluta X-Y en el esquemático
                'WINDOW 0 -8 55 VBottom 2\n', # posiciones relativas de etiquetas
                'WINDOW 3 -37 0 VTop 2\n', # posiciones relativas de etiquetas
                'SYMATTR InstName {:s}\n', # etiqueta que tendrá
                'SYMATTR Value {:3.5f}\n' # valor que tendrá
               ]  
"""
capacitor en serie
"""


#############################################
#%% Funciones para dibujar redes en LTspice #
#############################################

def ltsp_nuevo_circuito(circ_name=None):
    '''
    Se genera un circuito nuevo en LTspice de nombre *circ_name*.

    Parameters
    ----------
    circ_name : string
        Nombre del circuito.


    Returns
    -------
    circ_hdl : archivo de texto
        Handle al archivo de texto


    Raises
    ------
    TypeError
        Si ZZ no es una instancia de sympy.Matrix.


    See Also
    --------
    :func:`ltsp_capa_derivacion`
    :func:`ltsp_ind_serie`


    Examples
    --------
    >>> from pytc2.ltspice import ltsp_nuevo_circuito, ltsp_etiquetar_nodo, ltsp_ind_serie, ltsp_capa_derivacion, ltsp_etiquetar_nodo
    >>> circ_hdl = ltsp_nuevo_circuito('prueba1')
    >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vi')
    >>> ltsp_ind_serie(circ_hdl, 1.0) 
    >>> ltsp_capa_derivacion(circ_hdl, 2.0) 
    >>> ltsp_ind_serie(circ_hdl, 1.0) 
    >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vo')
    >>> R01 = 1.0
    >>> R02 = 1.0
    >>> circ_hdl.writelines('TEXT -48 304 Left 2 !.param RG={:3.3f} RL={:3.3f}'.format(R01, R02))
    >>> circ_hdl.close()
    [ Buscar el archivo "ltsp_prueba.asc" en LTspice ]


    '''
    
    if not isinstance(circ_name, (str, type(None)) ):
        raise ValueError('El nombre del circuito debe ser un string u omitirse.')

    global cap_num, res_num, ind_num, cur_x, cur_y

    if circ_name is None:

        timestr = time.strftime("%Y%m%d-%H%M%S")
        circ_name = 'NN-' + timestr

    circ_hdl = None
    
    src_fname = path.join('.', filename_eq_base + '.asc' )
    
    if path.isfile(src_fname):
        
        dst_fname = 'pyltspice_{:s}.asc'.format(circ_name)
        
        shutil.copy(src_fname, dst_fname)
    
        # configuración de los gráficos standard S11 / S21
        src_fname = path.join('.', filename_eq_base + '.plt' )
        
        dst_fname = 'pyltspice_{:s}.plt'.format(circ_name)
        
        shutil.copy(src_fname, dst_fname)
        
        circ_hdl = open(dst_fname, 'a')
        
        cap_num = 1
        res_num = 1
        ind_num = 1

        cur_x = 0
        cur_y = 0
        
    else:
        print("No se encontró el archivo {}".format(src_fname))
    
    return(circ_hdl)

def ltsp_capa_derivacion(circ_hdl, cap_value, cap_label=None):
    '''
    Incorpora un capacitor en derivación a un circuito en LTspice.

    Parameters
    ----------
    circ_hdl : archivo de texto LTspice
        Handle al archivo LTspice.
    cap_value : float o numéro simbólico
        Valor del capacitor.
    cap_label : string o None
        Etiqueta para identificar al capacitor en el circuito.


    Returns
    -------
    None


    Raises
    ------
    ValueError
        Si cap_value no es numérico o el valor no es positivo.


    See Also
    --------
    :func:`ltsp_capa_derivacion`
    :func:`ltsp_ind_serie`


    Examples
    --------
    >>> from pytc2.ltspice import ltsp_nuevo_circuito, ltsp_etiquetar_nodo, ltsp_ind_serie, ltsp_capa_derivacion, ltsp_etiquetar_nodo
    >>> circ_hdl = ltsp_nuevo_circuito('prueba1')
    >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vi')
    >>> ltsp_ind_serie(circ_hdl, 1.0) 
    >>> ltsp_capa_derivacion(circ_hdl, 2.0) 
    >>> ltsp_ind_serie(circ_hdl, 1.0) 
    >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vo')
    >>> R01 = 1.0
    >>> R02 = 1.0
    >>> circ_hdl.writelines('TEXT -48 304 Left 2 !.param RG={:3.3f} RL={:3.3f}'.format(R01, R02))
    >>> circ_hdl.close()
    [ Buscar el archivo "ltsp_prueba.asc" en LTspice ]

    '''
    
    if not ( isinstance(cap_value, (Real, sp.Number) ) and cap_value > 0 ):
        raise ValueError('Se espera un valor numérico positivo para el capacitor.')
    
    global cap_der_str, cap_num
    
    if cap_label is None:
        
        cap_label = 'C{:d}'.format(cap_num)
        cap_num += 1

    this_cap_str = cap_der_str.copy()
    
    element_xy = [cur_x - ltux, cur_y + lt_wire_length*ltuy]
    
    this_cap_str[0] = this_cap_str[0].format(element_xy[0], element_xy[1])
    this_cap_str[3] = this_cap_str[3].format(cap_label)
    this_cap_str[4] = this_cap_str[4].format(cap_value)

    # conectamos el elemento en derivación con el cursor actual.
    wire_str = 'WIRE {:d} {:d} {:d} {:d}\n'.format(cur_x, cur_y, element_xy[0] + ltux, element_xy[1] )
    # y el otro extremo a referencia GND
    gnd_str = 'FLAG {:d} {:d} 0\n'.format(element_xy[0] + ltux, element_xy[1] + 4*ltuy )

    circ_hdl.writelines(wire_str)
    circ_hdl.writelines(this_cap_str)
    circ_hdl.writelines(gnd_str)
    
    return()



def ltsp_ind_serie(circ_hdl, ind_value, ind_label=None):
    '''
    Incorpora un inductor en serie a un circuito en LTspice.

    Parameters
    ----------
    circ_hdl : archivo de texto LTspice
        Handle al archivo LTspice.
    ind_value : float o numéro simbólico
        Valor del inductor.
    ind_label : string o None
        Etiqueta para identificar al inductor en el circuito.


    Returns
    -------
    None


    Raises
    ------
    ValueError
        Si cap_value no es numérico o el valor no es positivo.


    See Also
    --------
    :func:`ltsp_capa_derivacion`
    :func:`ltsp_ind_serie`


    Examples
    --------
    >>> from pytc2.ltspice import ltsp_nuevo_circuito, ltsp_etiquetar_nodo, ltsp_ind_serie, ltsp_capa_derivacion, ltsp_etiquetar_nodo
    >>> circ_hdl = ltsp_nuevo_circuito('prueba1')
    >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vi')
    >>> ltsp_ind_serie(circ_hdl, 1.0) 
    >>> ltsp_capa_derivacion(circ_hdl, 2.0) 
    >>> ltsp_ind_serie(circ_hdl, 1.0) 
    >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vo')
    >>> R01 = 1.0
    >>> R02 = 1.0
    >>> circ_hdl.writelines('TEXT -48 304 Left 2 !.param RG={:3.3f} RL={:3.3f}'.format(R01, R02))
    >>> circ_hdl.close()
    [ Buscar el archivo "ltsp_prueba.asc" en LTspice ]

    '''

    if not ( isinstance(ind_value, (Real, sp.Number) ) and ind_value > 0 ):
        raise ValueError('Se espera un valor numérico positivo para el inductor.')
    
    global ind_ser_str, cap_num, cur_x, cur_y
    
    if ind_label is None:
        
        ind_label = 'C{:d}'.format(cap_num)
        cap_num += 1


    this_ind_str = ind_ser_str.copy()
    
    element_xy = [cur_x + lt_wire_length*ltux, cur_y + ltuy]
    
    this_ind_str[0] = this_ind_str[0].format(element_xy[0], element_xy[1])
    this_ind_str[3] = this_ind_str[3].format(ind_label)
    this_ind_str[4] = this_ind_str[4].format(ind_value)

    # conectamos el elemento en serie con el cursor actual, y el otro extremo
    # al siguiente elemento.
    
    next_x = element_xy[0] + 6*ltux + lt_wire_length*ltux
    next_y = element_xy[1] - ltuy
    
    wire_str = ['WIRE {:d} {:d} {:d} {:d}\n'.format(cur_x, cur_y, element_xy[0] + ltux, element_xy[1] - ltuy), 
                'WIRE {:d} {:d} {:d} {:d}\n'.format(element_xy[0] + 6*ltux, element_xy[1] - ltuy, next_x, next_y) ]

    # actualizamos cursor.    
    cur_x = next_x
    cur_y = next_y

    circ_hdl.writelines(wire_str)
    circ_hdl.writelines(this_ind_str)
    
    return()


def ltsp_etiquetar_nodo(circ_hdl, node_label=None):
    '''
    Asigna una etiqueta a un nodo de un circuito en LTspice.

    Parameters
    ----------
    circ_hdl : archivo de texto LTspice
        Handle al archivo LTspice.
    node_label : string o None
        Etiqueta para identificar al nodo en el circuito.


    Returns
    -------
    None


    Raises
    ------
    ValueError
        Si cap_value no es numérico o el valor no es positivo.


    See Also
    --------
    :func:`ltsp_capa_derivacion`
    :func:`ltsp_ind_serie`


    Examples
    --------
    >>> from pytc2.ltspice import ltsp_nuevo_circuito, ltsp_etiquetar_nodo, ltsp_ind_serie, ltsp_capa_derivacion, ltsp_etiquetar_nodo
    >>> circ_hdl = ltsp_nuevo_circuito('prueba1')
    >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vi')
    >>> ltsp_ind_serie(circ_hdl, 1.0) 
    >>> ltsp_capa_derivacion(circ_hdl, 2.0) 
    >>> ltsp_ind_serie(circ_hdl, 1.0) 
    >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vo')
    >>> R01 = 1.0
    >>> R02 = 1.0
    >>> circ_hdl.writelines('TEXT -48 304 Left 2 !.param RG={:3.3f} RL={:3.3f}'.format(R01, R02))
    >>> circ_hdl.close()
    [ Buscar el archivo "ltsp_prueba.asc" en LTspice ]

    '''
    
    global cap_der_str, node_num, cur_x, cur_y
    
    if node_label is None:
        
        node_label = 'v{:d}'.format(node_num)
        node_num += 1

    flag_str = ['FLAG {:d} {:d} {:s}\n'.format(cur_x, cur_y, node_label) ]

    circ_hdl.writelines(flag_str)
    
    return()


