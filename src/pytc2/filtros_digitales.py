 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mariano
"""

import numpy as np
from numbers import Integral, Real

import matplotlib.pyplot as plt

from scipy.signal import find_peaks

import warnings


#%%
   ##########################################
  ## Variables para el análisis simbólico #
 ##########################################
#%%

# from .general import s, small_val
"""
Variable compleja de Laplace s = σ + j.ω
En caso de necesitar usarla, importar el símbolo fr_desiredde este módulo.
"""

#%%
  ##############################################
 ## Variables para el funcionamiento general ##
##############################################
#%%

# phase_change_thr = 3/5*np.pi
"""
Representa la máxima variación en una función de fase de muestra a muestra.
Es útil para detectar cuando una función atraviesa un cero y se produce un
salto de :math:`\\pi` radianes. Por cuestiones de muestreo y de variaciones 
muy rápidas de fase, a veces conviene que sea un poco menor a :math:`\\pi`.
"""


#%%
  #########################
 ## Funciones generales ##
#########################
#%%

def fir_design_ls(order, band_edges, desired, weight = None, grid_density = 16, 
                  fs = 2.0, filter_type = 'multiband'):
    """
    Algoritmo de Parks-McClellan para el diseño de filtros FIR de fase lineal
    utilizando un criterio minimax. El algoritmo está basado en RERMEZ_FIR de 
    :ref:`Tapio Saramaki y Lars Whannamar <DSPMatlab20>` y el detallado trabajo
    en el material suplementario de :ref:`Thomas Holton <holton21>`. La imple_
    mentación del algoritmo ha sido ampliamente modificada con fines didácticos
    respecto a la version original de Saramaki y Parks McClellan.
    
    Parameters
    -----------
    order : TransferFunction
        Orden del filtro a diseñar. El tamaño del filtro será de *orden+1*.
    band_edges : array_like
        Los límites de cada banda indicada en la plantilla de diseño del filtro.
        Habrá dos valores, principio y fin, por cada banda definida en *fr_desiredired*.
        Ej: [0., 0.3, 0.7, 1.] Para un pasabajos con corte en 0.3
    fr_desiredired : array_like
        El valor numérico fr_desiredado por cada banda. Ej: [1.0, 0.] para un pasabajos.
    weight : array_like
        Un valor postivo que pesará cada banda al momento de calcular el error.
    grid_density : int, numeric
        Un entero que indicará por cuanto interpolar la respuesta del filtro al
        calcular el error del filtro. El valor de interpolación se calcula 
        *aproximadamente* por grid_density*orden/2. Por defecto se usa 16.
    fs : float, numeric
        Frecuencia de muestreo a la que se implementará el filtro digital. Por
        defecto se usa 2.0, es decir se normaliza a la f. de Nyquist.
    filter_type : string, 
        Un string que identifica el filtro que se diseñará. Se admiten tres 
        posibilidafr_desired: 'multiband' o 'm'. Filtros FIR tipo 1 o 2 de propósitos 
        generales. 'differentiator' o 'd', se utilizará para diseñar filtro FIR 
        derivadores de tipo 3 o 4 dependiendo el orden. Finalmente, 'hilbert' o
        'h' para implementar filtros FIR que permiten calcular la parte 
        imaginaria de una señal analítica. Es decir tener una transferencia 
        aproximadamente constante y una rotación constante de pi/2 para todas 
        las frecuencias.
    max_iter : int, numeric
        Cantidad máxima de iteraciones del algoritmo de Remez para hallar las 
        frecuencias extremas.
    debug : boolean
        Un valor booleano para activar la depuración de la propia función.
    
    order - filter order
    band_edges  - specifies the upper and lower band_edgess of the bands under consideration.
            The program, however, uses band efr_desired in terms of fractions of pi rad.
    	    band_edges = band_edges/pi;
    fr_desiredired -    specifies the fr_desiredired values at the band_edgess of each band.

    Returns
    --------
    h_coeffs : array_like
        Los coeficientes de la respuesta al impulso del filtro FIR diseñado.
    err : float, numeric
        Error máximo obtenido de la iteración del algoritmo Remez.
    w_extremas : array_like
        Las frecuencias extremas obtenidas de la iteración del algoritmo Remez.

    Raises
    ------
    ValueError
        Si no se cumple con el formato y valores indicados en la documentación.

    See Also
    -----------
    :func:``
    :func:``

    Examples
    --------
    >>> 
    >>> 
    >>> 
    >>> 
    >>> 

    Notes:
    -------
    .. _pm73:
        
    J. H. McClellan, T. W. Parks, and L. R. Rabiner, "A computer program for fr_desiredigning optimum FIR linear phase digital filters," IEEE Transactions on Audio and Electroacoustics, vol. AU-21, no. 6, pp. 506 - 526, December 1973.
    .. _DSPMatlab20:

    L. Wanhammar, T. Saramäki. Digital Filters Using MATLAB. Springer 2020.
    M. Ahsan and T. Saramäki, "A MATLAB based optimum multiband FIR filters fr_desiredign program following the original idea of the Remez multiple exchange algorithm," in Proc. 2011 IEEE International Symposium on Circuits and Systems, Rio de Janeiro, Brazil, May 15-17, 2011, pp. 137-140. 
    .. _holton21:

    T. Holton, Digital Signal Processing: Principles and Applications. Cambridge University Press, 2021. 
    	
    """

    if not (isinstance(order, (Integral, Real)) and order > 0 ):
        raise ValueError("El argumento 'order' debe ser un número positivo.")
                  
    if not (isinstance(grid_density, (Integral, Real)) and grid_density > 0 ):
        raise ValueError("El argumento 'grid_density' debe ser un número positivo.")

    if not (isinstance(fs, Real) and fs > 0 ):
        raise ValueError("El argumento 'fs' debe ser un número positivo.")

    valid_filters = ['multiband', 'lowpass', 'highpass', 'bandpass', 
                     'bandstop', 'notch',
                     'h', 'd', 'm', 'lp', 'hp', 'lp', 'bp',
                     'differentiator', 'hilbert']

    if not isinstance(filter_type, str):
        raise ValueError("El argumento 'filter_type' debe ser un string de %s" % (valid_filters))

	#==========================================================================
	#  Find out jtype that was used in the PM code.
	#  This not necessary but simplifies the undertanding of this code snippet.
	#==========================================================================
    if filter_type.lower().startswith('d'):
        jtype = 2  # Differentiator
    elif 'hilbert' == filter_type.lower() or 'h' == filter_type.lower():
        jtype = 3  # Hilbert transformer
    else: 
        jtype = 1  # Multiband filter

	#==========================================================================
	# Determine the filter cases and cant_bases, the number of basis functions to be 
	# used in the Remez algorithm 
	# In the below, filtercase=1,2,3,4 is used for making it easier to 
	# understand this code snippet.   
	#==========================================================================
    # Determine the filter cases and cant_bases
    if jtype == 1:
        if order % 2 == 0:
            filtercase = 1  # Even order and even symmetry multiband filter
        else:
            filtercase = 2  # Odd order and even symmetry multiband filter 
    else:
        if order % 2 == 0:
            # Even order and odd symmetry -> a Hilbert transforer or a 
            # differentiator (jtype indicates)
            filtercase = 3  
        else:
            # Odd order and odd symmetry -> a Hilbert transforer or a 
            # differentiator (jtype indicates)
            filtercase = 4  

    if filter_type not in valid_filters:
        raise ValueError('El tipo de filtro debe ser uno de %s, no %s'
                         % (valid_filters, filter_type))

    if not isinstance(band_edges, (list, np.ndarray)):
        raise ValueError("El argumento 'band_edges' debe ser una lista o un array de numpy.")

    if not isinstance(desired, (list, np.ndarray)):
        raise ValueError("El argumento 'fr_desiredired' debe ser una lista o un array de numpy.")
    
    if not isinstance(weight, (type(None), list, np.ndarray)):
        raise ValueError("El argumento 'weight' debe ser una lista o un array de numpy.")

    # Chequear si la plantilla de requerimientos del filtro está bien armada.
    ndesired = len(desired)
    nedges = len(band_edges)
    nbands = nedges // 2

    if isinstance(weight ,type(None)):
        weight = np.ones(nbands)
        
    if isinstance(weight, list):
        weight = np.array(weight)

    nweights = len(weight)

    if ndesired != nedges:
        raise ValueError(f"Debe haber tantos elementos en 'fr_desired' {ndesired} como en 'band_edges' {nedges}")

    if jtype == 1:
        # multibanda 
        if nweights != nbands:
            raise ValueError(f"Debe haber tantos elementos en 'weight' {nweights} como cantidad de bandas {nbands}")

    if jtype == 2 or jtype == 3:
        # derivador y hilbert
        if nbands != 1:
            raise ValueError(f"Debe haber en una sola banda definida para FIR tipo {filter_type}, hay {nbands} bandas")

    # normalizar respecto a Nyquist
    band_edges = np.array(band_edges) / (fs/2)

    desired = np.array(desired)

            
	# cant_bases - number of basis functions 
    cant_coeffs = order + 1
    
    if filtercase == 1 or filtercase == 3:
        M = (cant_coeffs-1) // 2
        # cantidad de frecuencias extremas.
        cant_bases = M + 1
        
    if filtercase == 2 or filtercase == 4:
        M = cant_coeffs // 2
        # cantidad de frecuencias extremas.
        cant_bases = M 

    
    # propuesta original de Tapio
    # cant_bases = (np.fix((order + 1)/ 2)).astype(int)
    # if filtercase == 1:
    #     cant_bases += 1
    
    
	#=========================================================================
	# DETERMINE fr_grid, fr_desired, and fr_weight 
	#========================================================================
	# Compared with the PM code, there are the following key differences:
	# (1) The upper band_edges for each band under consideration is automatically 
	#     included in fr_grid. This somehow increases the accuracy. 
	# (2) Since the frequency range is now from 0 to 1, freq_resolution has been increased
	#     by a factor of 2.
	# (3) The change of fr_desired and fr_weight depending on the filter type is peformed 
	#     before using the (modified) Remez algorithm.
	# (4) The removal of problematic angular frequencies at 0 and pi is 
	#     performed simultaneously for all filter types. Now the remomal is
	#     is performed while generating fr_grid.
	#=========================================================================
    
    # Determine fr_grid, fr_desired, and fr_weight
    freq_resolution = 1.0 / (grid_density * cant_bases)
    # full resolution (fr) fr_grid, desired and wieight arrays
    fr_grid = []
    fr_desired = []
    fr_weight = []
    # indexes of the band-edges corresponding to the fr freq. fr_grid array
    band_edges_idx = []

    for ll in range(nbands):
        number_fr_grid = int(np.ceil((band_edges[2 * ll + 1] - band_edges[2 * ll]) / freq_resolution))
        fr_grid_more = np.linspace(band_edges[2 * ll], band_edges[2 * ll + 1], number_fr_grid + 1)
        
        # Adjust fr_grid for harmful frequencies at omega = 0 
        if ll == 0 and (filtercase == 3 or filtercase == 4) and fr_grid_more[0] < freq_resolution:
            fr_grid_more = fr_grid_more[1:]
            number_fr_grid -= 1

        # Adjust fr_grid for harmful frequencies at omega = 1
        if ll == nbands - 1 and (filtercase == 2 or filtercase == 3) and fr_grid_more[-1] > 1 - freq_resolution:
            fr_grid_more = fr_grid_more[:-1]
            number_fr_grid -= 1

        #
        band_edges_idx.extend([len(fr_grid)])
        fr_grid.extend(fr_grid_more)
        band_edges_idx.extend([len(fr_grid)-1])

        if jtype == 2:
            # differentiator
            
            des_more = desired[2*ll+1] * fr_grid_more * np.pi
            
            if np.abs(desired[2*ll]) < 1.0e-3:
                wt_more = weight[ll] * np.ones(number_fr_grid + 1)
            else:
                wt_more = weight[ll] / (fr_grid_more * np.pi)
        else:
            # others

            wt_more = weight[ll] * np.ones(number_fr_grid + 1)
            if desired[2 * ll + 1] != desired[2 * ll]:
                des_more = np.linspace(desired[2 * ll], desired[2 * ll + 1], number_fr_grid + 1)
            else:
                des_more = desired[2 * ll] * np.ones(number_fr_grid + 1)

        fr_desired.extend(des_more)
        fr_weight.extend(wt_more)

    fr_grid = np.array(fr_grid)
    fr_desired = np.array(fr_desired)
    fr_weight = np.array(fr_weight)
    band_edges_idx = np.array(band_edges_idx)

	#==========================================================================
	# Modify fr_desired and fr_weight depending on the filter case
	#========================================================================== 
    # Este es un elegante truco para hacer una sola función de optimización
    # de Remez para todos los tipos de FIRs. 
    # Ver :ref:`Thomas Holton supplimentary material <holton21>`.
    # 
    
    if filtercase == 2:
        fr_desired /= np.cos(np.pi * fr_grid / 2)
        fr_weight *= np.cos(np.pi * fr_grid / 2)
    if filtercase == 4:
        fr_desired /= np.sin(np.pi * fr_grid / 2)
        fr_weight *= np.sin(np.pi * fr_grid / 2)
    if filtercase == 3:
        fr_desired /= np.sin(np.pi * fr_grid)
        fr_weight *= np.sin(np.pi * fr_grid)

    #==========================================================================
    # Resolvemos el sistema mediante LS  (CC'*WW*)a = CC'*WW*D
    # ver I. Selesnick 713 Lecture Notes: "LINEAR-PHASE FIR FILTER DESIGN BY 
    # LEAST SQUARES"
	#==========================================================================
    
    # cantidad de puntos donde se calcula la diferencia entre
    # la respuesta deseada y la obtenida
    R = len(fr_grid)
    
    # Construir la matriz de diseño A
    CC = np.zeros((R, cant_bases))
    
    for i,f in enumerate(fr_grid):
        CC[i, :] = np.cos( np.pi * f * np.arange(cant_bases) )
    
    WW = np.diag(np.sqrt(fr_weight))
    
    # Resolver el sistema de ecuaciones para los coeficientes únicos
    a_coeffs = np.linalg.lstsq( np.matmul(np.matmul(CC.transpose(), WW),CC) , np.matmul(np.matmul(CC.transpose(), WW), fr_desired), rcond=None)[0]


    #======================================================
    # Construir el filtro a partir de los coeficientes "a"
	#======================================================
    
    cant_acoeffs = len(a_coeffs)

    # convertir los coeficientes según el tipo de FIR
    if filtercase == 1:
        
        a_coeffs [1:] = a_coeffs[1:]/2
        h_coeffs = np.concatenate((a_coeffs[::-1], a_coeffs[1:]))
    
    if filtercase == 2:
        
        last_coeff = cant_acoeffs
        cant_hcoeff = 2*cant_acoeffs
        h_coeffs = np.zeros(cant_hcoeff)
        h_coeffs[cant_hcoeff-1] = a_coeffs[last_coeff-1]/4
        h_coeffs[last_coeff] = a_coeffs[0] /2 + a_coeffs[1]/4
        h_coeffs[last_coeff+1:cant_hcoeff-1]= (a_coeffs[1:last_coeff-1] + a_coeffs[2:last_coeff])/4
            
        h_coeffs[:last_coeff] = h_coeffs[last_coeff:][::-1]

        
    if filtercase == 3:
        
        cant_hcoeff = 2*cant_acoeffs+1
        h_coeffs = np.zeros(cant_hcoeff)
        last_coeff = cant_acoeffs # punto de simetría, demora del filtro


        h_coeffs[0:2] = a_coeffs[last_coeff-2:][::-1]/4
        h_coeffs[2:last_coeff-1] = ((a_coeffs[1:last_coeff-2] - a_coeffs[3:last_coeff])/4)[::-1]
        h_coeffs[last_coeff-1] = a_coeffs[0]/2 - a_coeffs[2]/4
        
        h_coeffs[last_coeff+1:] = (-1.)*h_coeffs[:last_coeff][::-1]

    if filtercase == 4:
        
        last_coeff = cant_acoeffs
        cant_hcoeff = 2*cant_acoeffs
        h_coeffs = np.zeros(2*cant_acoeffs)
        h_coeffs[cant_hcoeff-1] = a_coeffs[last_coeff-1]/4
        h_coeffs[last_coeff] = a_coeffs[0]/2 - a_coeffs[1]/4
        h_coeffs[last_coeff+1:cant_hcoeff-1]= (a_coeffs[1:last_coeff-1] - a_coeffs[2:last_coeff])/4
            
        h_coeffs[:last_coeff] = -1. * h_coeffs[last_coeff:][::-1]
    
    return h_coeffs

def fir_design_pm(order, band_edges, desired, weight = None, grid_density = 16, 
                  fs = 2.0, filter_type = 'multiband', max_iter = 25, debug=False):
    """
    Algoritmo de Parks-McClellan para el diseño de filtros FIR de fase lineal
    utilizando un criterio minimax. El algoritmo está basado en RERMEZ_FIR de 
    :ref:`Tapio Saramaki y Lars Whannamar <DSPMatlab20>` y el detallado trabajo
    en el material suplementario de :ref:`Thomas Holton <holton21>`. La imple_
    mentación del algoritmo ha sido ampliamente modificada con fines didácticos
    respecto a la version original de Saramaki y Parks McClellan.
    
    Parameters
    -----------
    order : TransferFunction
        Orden del filtro a diseñar. El tamaño del filtro será de *orden+1*.
    band_edges : array_like
        Los límites de cada banda indicada en la plantilla de diseño del filtro.
        Habrá dos valores, principio y fin, por cada banda definida en *fr_desiredired*.
        Ej: [0., 0.3, 0.7, 1.] Para un pasabajos con corte en 0.3
    fr_desiredired : array_like
        El valor numérico fr_desiredado por cada banda. Ej: [1.0, 0.] para un pasabajos.
    weight : array_like
        Un valor postivo que pesará cada banda al momento de calcular el error.
    grid_density : int, numeric
        Un entero que indicará por cuanto interpolar la respuesta del filtro al
        calcular el error del filtro. El valor de interpolación se calcula 
        *aproximadamente* por grid_density*orden/2. Por defecto se usa 16.
    fs : float, numeric
        Frecuencia de muestreo a la que se implementará el filtro digital. Por
        defecto se usa 2.0, es decir se normaliza a la f. de Nyquist.
    filter_type : string, 
        Un string que identifica el filtro que se diseñará. Se admiten tres 
        posibilidafr_desired: 'multiband' o 'm'. Filtros FIR tipo 1 o 2 de propósitos 
        generales. 'differentiator' o 'd', se utilizará para diseñar filtro FIR 
        derivadores de tipo 3 o 4 dependiendo el orden. Finalmente, 'hilbert' o
        'h' para implementar filtros FIR que permiten calcular la parte 
        imaginaria de una señal analítica. Es decir tener una transferencia 
        aproximadamente constante y una rotación constante de pi/2 para todas 
        las frecuencias.
    max_iter : int, numeric
        Cantidad máxima de iteraciones del algoritmo de Remez para hallar las 
        frecuencias extremas.
    debug : boolean
        Un valor booleano para activar la depuración de la propia función.
    
    order - filter order
    band_edges  - specifies the upper and lower band_edgess of the bands under consideration.
            The program, however, uses band efr_desired in terms of fractions of pi rad.
    	    band_edges = band_edges/pi;
    fr_desiredired -    specifies the fr_desiredired values at the band_edgess of each band.

    Returns
    --------
    h_coeffs : array_like
        Los coeficientes de la respuesta al impulso del filtro FIR diseñado.
    err : float, numeric
        Error máximo obtenido de la iteración del algoritmo Remez.
    w_extremas : array_like
        Las frecuencias extremas obtenidas de la iteración del algoritmo Remez.

    Raises
    ------
    ValueError
        Si no se cumple con el formato y valores indicados en la documentación.

    See Also
    -----------
    :func:``
    :func:``

    Examples
    --------
    >>> 
    >>> 
    >>> 
    >>> 
    >>> 

    Notes:
    -------
    .. _pm73:
        
    J. H. McClellan, T. W. Parks, and L. R. Rabiner, "A computer program for fr_desiredigning optimum FIR linear phase digital filters," IEEE Transactions on Audio and Electroacoustics, vol. AU-21, no. 6, pp. 506 - 526, December 1973.
    .. _DSPMatlab20:

    L. Wanhammar, T. Saramäki. Digital Filters Using MATLAB. Springer 2020.
    M. Ahsan and T. Saramäki, "A MATLAB based optimum multiband FIR filters fr_desiredign program following the original idea of the Remez multiple exchange algorithm," in Proc. 2011 IEEE International Symposium on Circuits and Systems, Rio de Janeiro, Brazil, May 15-17, 2011, pp. 137-140. 
    .. _holton21:

    T. Holton, Digital Signal Processing: Principles and Applications. Cambridge University Press, 2021. 
    	
    """

    if not (isinstance(order, (Integral, Real)) and order > 0 ):
        raise ValueError("El argumento 'order' debe ser un número positivo.")
                  
    if not (isinstance(grid_density, (Integral, Real)) and grid_density > 0 ):
        raise ValueError("El argumento 'grid_density' debe ser un número positivo.")
                  
    if not (isinstance(max_iter, (Integral, Real)) and max_iter > 0 ):
        raise ValueError("El argumento 'max_iter' debe ser un número positivo.")

    if not isinstance(debug, bool):
        raise ValueError('displaystr debe ser un booleano')

    if not (isinstance(fs, Real) and fs > 0 ):
        raise ValueError("El argumento 'fs' debe ser un número positivo.")

    valid_filters = ['multiband', 'lowpass', 'highpass', 'bandpass', 
                     'bandstop', 'notch',
                     'h', 'd', 'm', 'lp', 'hp', 'lp', 'bp',
                     'differentiator', 'hilbert']

    if not isinstance(filter_type, str):
        raise ValueError("El argumento 'filter_type' debe ser un string de %s" % (valid_filters))

	#==========================================================================
	#  Find out jtype that was used in the PM code.
	#  This not necessary but simplifies the undertanding of this code snippet.
	#==========================================================================
    if filter_type.lower().startswith('d'):
        jtype = 2  # Differentiator
    elif 'hilbert' == filter_type.lower() or 'h' == filter_type.lower():
        jtype = 3  # Hilbert transformer
    else: 
        jtype = 1  # Multiband filter

	#==========================================================================
	# Determine the filter cases and cant_bases, the number of basis functions to be 
	# used in the Remez algorithm 
	# In the below, filtercase=1,2,3,4 is used for making it easier to 
	# understand this code snippet.   
	#==========================================================================
    # Determine the filter cases and cant_bases
    if jtype == 1:
        if order % 2 == 0:
            filtercase = 1  # Even order and even symmetry multiband filter
        else:
            filtercase = 2  # Odd order and even symmetry multiband filter 
    else:
        if order % 2 == 0:
            # Even order and odd symmetry -> a Hilbert transforer or a 
            # differentiator (jtype indicates)
            filtercase = 3  
        else:
            # Odd order and odd symmetry -> a Hilbert transforer or a 
            # differentiator (jtype indicates)
            filtercase = 4  

    if filter_type not in valid_filters:
        raise ValueError('El tipo de filtro debe ser uno de %s, no %s'
                         % (valid_filters, filter_type))

    if not isinstance(band_edges, (list, np.ndarray)):
        raise ValueError("El argumento 'band_edges' debe ser una lista o un array de numpy.")

    if not isinstance(desired, (list, np.ndarray)):
        raise ValueError("El argumento 'fr_desiredired' debe ser una lista o un array de numpy.")
    
    if not isinstance(weight, (type(None), list, np.ndarray)):
        raise ValueError("El argumento 'weight' debe ser una lista o un array de numpy.")

    # Chequear si la plantilla de requerimientos del filtro está bien armada.
    ndesired = len(desired)
    nedges = len(band_edges)
    nbands = nedges // 2

    if isinstance(weight ,type(None)):
        weight = np.ones(nbands)
        
    if isinstance(weight, list):
        weight = np.array(weight)

    nweights = len(weight)

    if ndesired != nedges:
        raise ValueError(f"Debe haber tantos elementos en 'fr_desired' {ndesired} como en 'band_edges' {nedges}")

    if jtype == 1:
        # multibanda 
        if nweights != nbands:
            raise ValueError(f"Debe haber tantos elementos en 'weight' {nweights} como cantidad de bandas {nbands}")

    if jtype == 2 or jtype == 3:
        # derivador y hilbert
        if nbands != 1:
            raise ValueError(f"Debe haber en una sola banda definida para FIR tipo {filter_type}, hay {nbands} bandas")

    # normalizar respecto a Nyquist
    band_edges = np.array(band_edges) / (fs/2)

            
	# cant_bases - number of basis functions 
    cant_coeffs = order + 1
    
    if filtercase == 1 or filtercase == 3:
        M = (cant_coeffs-1) // 2
        # cantidad de frecuencias extremas.
        cant_bases = M + 1
        
    if filtercase == 2 or filtercase == 4:
        M = cant_coeffs // 2
        # cantidad de frecuencias extremas.
        cant_bases = M 

    
    # propuesta original de Tapio
    # cant_bases = (np.fix((order + 1)/ 2)).astype(int)
    # if filtercase == 1:
    #     cant_bases += 1
    
    
	#=========================================================================
	# DETERMINE fr_grid, fr_desired, and fr_weight 
	#========================================================================
	# Compared with the PM code, there are the following key differences:
	# (1) The upper band_edges for each band under consideration is automatically 
	#     included in fr_grid. This somehow increases the accuracy. 
	# (2) Since the frequency range is now from 0 to 1, freq_resolution has been increased
	#     by a factor of 2.
	# (3) The change of fr_desired and fr_weight depending on the filter type is peformed 
	#     before using the (modified) Remez algorithm.
	# (4) The removal of problematic angular frequencies at 0 and pi is 
	#     performed simultaneously for all filter types. Now the remomal is
	#     is performed while generating fr_grid.
	#=========================================================================
    
    # Determine fr_grid, fr_desired, and fr_weight
    freq_resolution = 1.0 / (grid_density * cant_bases)
    # full resolution (fr) fr_grid, desired and wieight arrays
    fr_grid = []
    fr_desired = []
    fr_weight = []
    # indexes of the band-edges corresponding to the fr freq. fr_grid array
    band_edges_idx = []

    for ll in range(nbands):
        number_fr_grid = int(np.ceil((band_edges[2 * ll + 1] - band_edges[2 * ll]) / freq_resolution))
        fr_grid_more = np.linspace(band_edges[2 * ll], band_edges[2 * ll + 1], number_fr_grid + 1)
        
        # Adjust fr_grid for harmful frequencies at omega = 0 
        if ll == 0 and (filtercase == 3 or filtercase == 4) and fr_grid_more[0] < freq_resolution:
            fr_grid_more = fr_grid_more[1:]
            number_fr_grid -= 1

        # Adjust fr_grid for harmful frequencies at omega = 1
        if ll == nbands - 1 and (filtercase == 2 or filtercase == 3) and fr_grid_more[-1] > 1 - freq_resolution:
            fr_grid_more = fr_grid_more[:-1]
            number_fr_grid -= 1

        #
        band_edges_idx.extend([len(fr_grid)])
        fr_grid.extend(fr_grid_more)
        band_edges_idx.extend([len(fr_grid)-1])

        if jtype == 2:
            # differentiator
            
            des_more = desired[2*ll+1] * fr_grid_more * np.pi
            
            if np.abs(desired[2*ll]) < 1.0e-3:
                wt_more = weight[ll] * np.ones(number_fr_grid + 1)
            else:
                wt_more = weight[ll] / (fr_grid_more * np.pi)
        else:
            # others

            wt_more = weight[ll] * np.ones(number_fr_grid + 1)
            if desired[2 * ll + 1] != desired[2 * ll]:
                des_more = np.linspace(desired[2 * ll], desired[2 * ll + 1], number_fr_grid + 1)
            else:
                des_more = desired[2 * ll] * np.ones(number_fr_grid + 1)

        fr_desired.extend(des_more)
        fr_weight.extend(wt_more)

    fr_grid = np.array(fr_grid)
    fr_desired = np.array(fr_desired)
    fr_weight = np.array(fr_weight)
    band_edges_idx = np.array(band_edges_idx)

	#==========================================================================
	# Modify fr_desired and fr_weight depending on the filter case
	#========================================================================== 
    # Este es un elegante truco para hacer una sola función de optimización
    # de Remez para todos los tipos de FIRs. 
    # Ver :ref:`Thomas Holton supplimentary material <holton21>`.
    # 
    if filtercase == 2:
        fr_desired /= np.cos(np.pi * fr_grid / 2)
        fr_weight *= np.cos(np.pi * fr_grid / 2)
    if filtercase == 4:
        fr_desired /= np.sin(np.pi * fr_grid / 2)
        fr_weight *= np.sin(np.pi * fr_grid / 2)
    if filtercase == 3:
        fr_desired /= np.sin(np.pi * fr_grid)
        fr_weight *= np.sin(np.pi * fr_grid)

    #==========================================================================
	# CALL THE REMEZ ALGORITHM 
	#==========================================================================

    a_coeffs, err, w_extremas = _remez_exchange_algorithm(cant_bases, fr_grid, fr_desired, fr_weight, band_edges_idx, max_iter = max_iter, debug=debug)
    
    # por ahora si aparecen más picos los limito a la hora de construir el filtro
    # a_coeffs = a_coeffs[:cant_bases]
    
    cant_acoeffs = len(a_coeffs)

    # convertir los coeficientes según el tipo de FIR
    if filtercase == 1:
        
        a_coeffs [1:] = a_coeffs[1:]/2
        h_coeffs = np.concatenate((a_coeffs[::-1], a_coeffs[1:]))
    
    if filtercase == 2:
        
        last_coeff = cant_acoeffs
        cant_hcoeff = 2*cant_acoeffs
        h_coeffs = np.zeros(cant_hcoeff)
        h_coeffs[cant_hcoeff-1] = a_coeffs[last_coeff-1]/4
        h_coeffs[last_coeff] = a_coeffs[0] /2 + a_coeffs[1]/4
        h_coeffs[last_coeff+1:cant_hcoeff-1]= (a_coeffs[1:last_coeff-1] + a_coeffs[2:last_coeff])/4
            
        h_coeffs[:last_coeff] = h_coeffs[last_coeff:][::-1]

        
    if filtercase == 3:
        
        cant_hcoeff = 2*cant_acoeffs+1
        h_coeffs = np.zeros(cant_hcoeff)
        last_coeff = cant_acoeffs # punto de simetría, demora del filtro


        h_coeffs[0:2] = a_coeffs[last_coeff-2:][::-1]/4
        h_coeffs[2:last_coeff-1] = ((a_coeffs[1:last_coeff-2] - a_coeffs[3:last_coeff])/4)[::-1]
        h_coeffs[last_coeff-1] = a_coeffs[0]/2 - a_coeffs[2]/4
        
        h_coeffs[last_coeff+1:] = (-1.)*h_coeffs[:last_coeff][::-1]

    if filtercase == 4:
        
        last_coeff = cant_acoeffs
        cant_hcoeff = 2*cant_acoeffs
        h_coeffs = np.zeros(2*cant_acoeffs)
        h_coeffs[cant_hcoeff-1] = a_coeffs[last_coeff-1]/4
        h_coeffs[last_coeff] = a_coeffs[0]/2 - a_coeffs[1]/4
        h_coeffs[last_coeff+1:cant_hcoeff-1]= (a_coeffs[1:last_coeff-1] - a_coeffs[2:last_coeff])/4
            
        h_coeffs[:last_coeff] = -1. * h_coeffs[last_coeff:][::-1]

    err = np.abs(err)
    
    return h_coeffs, err, w_extremas

   ########################
  ## Funciones internas #
 ########################
#%%


# Función para filtrar los extremos consecutivos de mismo signo y mantener el de mayor módulo absoluto
def _filter_extremes(Ew, peaks):
    filtered_peaks = []
    current_sign = np.sign(Ew[peaks[0]])
    max_peak = peaks[0]
    
    for peak in peaks[1:]:
        peak_sign = np.sign(Ew[peak])
        
        # Si el signo del siguiente extremo es el mismo, conservamos el de mayor módulo absoluto
        if peak_sign == current_sign:
            if np.abs(Ew[peak]) > np.abs(Ew[max_peak]):
                max_peak = peak  # Actualizamos el pico con el mayor valor absoluto
        else:
            filtered_peaks.append(max_peak)  # Guardamos el pico de mayor valor absoluto del grupo
            max_peak = peak  # Empezamos a comparar en el nuevo grupo
            current_sign = peak_sign
    
    # Agregar el último extremo
    filtered_peaks.append(max_peak)
    
    return np.array(filtered_peaks)


def _remez_exchange_algorithm(cant_bases, fr_grid, fr_desired, fr_weight, band_edges_idx, max_iter = 250, error_tol = 10e-4, debug = False):
	# 	Function REMEZ_EX_MLLS implements the Remez exchange algorithm for the weigthed 
	#	Chebyshev approximation of a continous function with a sum of cosines.
	# Inputs
	#     cant_bases - number of basis functions 
	#     fr_grid - frequency fr_grid between 0 and 1
	#     fr_desired - fr_desiredired function on frequency fr_grid
	#     fr_weight - weight function on frequency fr_grid
	# Outputs
	#     h - coefficients of the filtercase = 1 filter
	#     dev - the resulting value of the weighted error function
	#     w_extremas - indices of extremal frequencies
    
    # Initializations
    nfr_grid = len(fr_grid)
    # l_ove = np.arange(nfr_grid)

    # Definir frecuencias extremas iniciales
    omega_scale = (nfr_grid - 1) / cant_bases
    jj = np.arange(cant_bases)
    omega_ext_iniciales_idx = np.concatenate((np.fix(omega_scale * jj), [nfr_grid-1])).astype(int)

    
    # aseguro que siempre haya una omega extrema en los band_edgess.
    aux_idx = np.array([np.argmin(np.abs(fr_grid[omega_ext_iniciales_idx] - fr_grid[ii])) for ii in band_edges_idx])
    omega_ext_iniciales_idx[aux_idx] = band_edges_idx

    cant_edges = len(band_edges_idx) 

    ## Debug

    fs = 2.0
    fft_sz = 512
    half_fft_sz = fft_sz//2
    frecuencias = np.arange(start=0, stop=fs, step=fs/fft_sz )

    if debug:
        ## Debug
        plt.figure(1)
        plt.clf()
        plt.figure(2)
        plt.clf()
        plt.figure(3)
        plt.clf()
        D_ext = np.interp(frecuencias[:half_fft_sz], fr_grid, fr_desired)
        plt.plot(frecuencias[:half_fft_sz], D_ext, label='D($\Omega$)')
        ## Debug
    
    niter = 1

    omega_ext_idx = omega_ext_iniciales_idx
    omega_ext_prev_idx = np.zeros_like(omega_ext_idx)

    cant_extremos_esperados = cant_bases+1
    cant_extremos = cant_extremos_esperados

    prev_error_target = np.finfo(np.float64).max
    
    # Remez loop
    while niter < max_iter:

        # Construir el sistema de ecuaciones a partir de la matriz de diseño A.
        A = np.zeros((cant_extremos, cant_extremos))
        for ii, omega_idx in enumerate(omega_ext_idx):
            A[ii,:] = np.hstack((np.cos( np.pi * fr_grid[omega_idx] * np.arange(cant_extremos-1)), (-1)**ii/fr_weight[omega_idx]))

        # Resolver el sistema de ecuaciones para los coeficientes únicos
        xx = np.linalg.solve(A, fr_desired[omega_ext_idx])
        
        # los primeros resultados están realacionados a los coeficientes del filtro
        a_coeffs_half = xx[:-1]
        # el último es el error cometido en la aproximación
        this_error_target = np.abs(xx[-1])

        # Construimos la respuesta interpolada en "fr_grid" para refinar las 
        # frecuencias extremas
        Aw_fr_grid = np.zeros(nfr_grid)
        for ii in range(cant_extremos-1):
            Aw_fr_grid  += a_coeffs_half[ii] * np.cos( ii * np.pi * fr_grid )

        # Calculamos la secuencia de error pesado: nos van a interesar los 
        # signos en las omega extremas para filtrar aquellas omega que NO 
        # alternan.
        Ew = fr_weight*(fr_desired - Aw_fr_grid)
        # también el módulo para verificar que ninguno esté por encima del 
        # error cometido "this_error_target"
        Ew_abs = np.abs(Ew)
        
        # procedemos a filtrar las omega extremas.
        peaks_pos , _ = find_peaks(Ew, height= 0.0)
        peaks_neg , _ = find_peaks(-Ew, height= 0.0)
        peaks = np.sort(np.concatenate((peaks_pos,peaks_neg)))
        
        # Aplicar el filtro a los picos encontrados
        peaks = _filter_extremes(Ew, peaks)

        omega_ext_idx = np.unique(np.concatenate((band_edges_idx, peaks)))

        omega_ext_idx = _filter_extremes(Ew, omega_ext_idx)
        
        cant_extremos = len(omega_ext_idx)

        # probamos si converge exitosamente
        if np.std(Ew_abs[omega_ext_idx] - this_error_target) < np.max(Ew_abs[omega_ext_idx]) * error_tol:
            
            print("Convergencia exitosa!")
            break

        # Problemas en la convergencia: sin cambios en el error ni las frecuencias extremas 
        elif this_error_target  == prev_error_target and np.array_equal(omega_ext_idx, omega_ext_prev_idx):
            warnings.warn("Problemas de convergencia: El error no disminuyó y ni cambiaron las frecuencias extremas.", UserWarning)
            break
        
        # Problemas en la convergencia: más extremos de los esperados
        elif cant_extremos > cant_extremos_esperados:
            # warnings.warn(f"Encontramos más extremos {cant_extremos}, de los que se esperaban {cant_extremos_esperados}. Extrarriple?", UserWarning)
            
            cant_extra = cant_extremos - cant_extremos_esperados

            if cant_extra % 2 == 1:
                # impar
                if Ew_abs[omega_ext_idx[0]] > Ew_abs[omega_ext_idx[-1]]:
                    #descarto el último
                    omega_ext_idx = omega_ext_idx[:-1]
                else:
                    #descarto el primero
                    omega_ext_idx = omega_ext_idx[1:]

            cant_extremos = len(omega_ext_idx)
            cant_extra = cant_extremos - cant_extremos_esperados

            while cant_extra > 0:

                Ew_abs_comp = np.hstack((Ew_abs[omega_ext_idx[:-2]],Ew_abs[omega_ext_idx[1:]]))
                
                if np.max( (Ew_abs[omega_ext_idx[0]], Ew_abs[omega_ext_idx[-1]])) <= np.min(Ew_abs_comp):
                    # descarto los extremos
                    omega_ext_idx = omega_ext_idx[1:-1]
                else:
                    # descarto el mínimo y su adyacente para no romper la alternancia de Remez.
                    min_idx = np.argmin(Ew_abs_comp)

                    omega_ext_idx = np.concatenate( (omega_ext_idx[:min_idx], omega_ext_idx[min_idx+2:] ) )
                          
                cant_extremos = len(omega_ext_idx)
                cant_extra = cant_extremos - cant_extremos_esperados
            

        if debug:
            ## Debug
            # Graficar la respuesta en frecuencia
            plt.figure(1)
            # plt.clf()
            # plt.plot(frecuencias[:half_fft_sz], Aw_ext, label=f'Aw_ext {niter}')
            # plt.plot(fr_grid[omega_ext_idx], Aw, 'ob')
            # plt.plot(frecuencias[:half_fft_sz], W_err_orig, label=f'orig {niter}')
        
            # plt.plot(fr_grid, Ew, label=f'$E_{niter}$')
            plt.plot(fr_grid, Ew)
            plt.plot(fr_grid[omega_ext_prev_idx], Ew[omega_ext_prev_idx], 'or')
            # plt.plot(frecuencias[:half_fft_sz], w_err_ext, label=f'Ew_ext {niter}')
            plt.plot(fr_grid[omega_ext_idx], Ew[omega_ext_idx], 'xb')
            plt.plot([ 0, 1], [0, 0], '-k', lw=0.8)
            plt.plot([ 0, 1], [this_error_target, this_error_target], ':k', lw=0.8, label=f'{cant_extremos} $\delta_{niter}=$ {this_error_target:3.3f}')
            plt.plot([ 0, 1], [-this_error_target, -this_error_target], ':k', lw=0.8)
        
            plt.title("Error pesado: $E(\Omega) = W(\Omega) \cdot [D(\Omega) - H_R(\Omega)]$")
            plt.xlabel("Frecuencia Normalizada")
            plt.ylabel("Magnitud")
            plt.legend()
        
            a_coeffs_half = xx[:-1]
            a_coeffs_half[1:] = a_coeffs_half[1:]/2
            h_coeffs = np.concatenate((a_coeffs_half[::-1], a_coeffs_half[1:]))
        
            H = np.fft.fft(h_coeffs, fft_sz)
        
            plt.figure(2)
            plt.plot(frecuencias[:half_fft_sz], 20*np.log10(np.abs(H[:half_fft_sz])), label=f'Iter: {niter}')
    
            plt.title("Respuesta en frecuencia de módulo: $ \\left|H(\Omega)\\right| $")
            plt.xlabel("Frecuencia Normalizada")
            plt.ylabel("$\\left|H(\Omega)\\right|_{{dB}}$")
            plt.legend()
        
            plt.figure(3)
            Aw_ext = np.interp(frecuencias[:half_fft_sz], fr_grid, Aw_fr_grid)
            plt.plot(frecuencias[:half_fft_sz], Aw_ext, label=f'$H_{{R{niter}}}$')
            plt.legend()
            plt.show()
            pass
    
            ## Debug

        # continuamos buscando la convergencia
        omega_ext_prev_idx = omega_ext_idx
        prev_error_target = this_error_target
        niter += 1


    # coeficientes del filtro        
    a_coeffs_half = xx[:-1]
    aux_val = a_coeffs_half.copy()
    aux_val [1:] = aux_val[1:]/2

    h_coeffs = np.concatenate((aux_val [::-1], aux_val [1:]))

    ## Debug
    if debug:
        # Graficar la respuesta en frecuencia
        plt.figure(1)
        # plt.clf()
        # plt.plot(frecuencias[:half_fft_sz], Aw_ext, label=f'Aw_ext {niter}')
        # plt.plot(fr_grid[omega_ext_idx], Aw, 'ob')
        # plt.plot(frecuencias[:half_fft_sz], W_err_orig, label=f'orig {niter}')
    
        # plt.plot(fr_grid, Ew, label=f'$E_{niter}$')
        plt.plot(fr_grid, Ew)
        plt.plot(fr_grid[omega_ext_prev_idx], Ew[omega_ext_prev_idx], 'or')
        # plt.plot(frecuencias[:half_fft_sz], w_err_ext, label=f'Ew_ext {niter}')
        plt.plot(fr_grid[omega_ext_idx], Ew[omega_ext_idx], 'xb')
        plt.plot([ 0, 1], [0, 0], '-k', lw=0.8)
        plt.plot([ 0, 1], [this_error_target, this_error_target], ':k', lw=0.8, label=f'{cant_extremos} $\delta_{niter}=$ {this_error_target:3.3f}')
        plt.plot([ 0, 1], [-this_error_target, -this_error_target], ':k', lw=0.8)
    
        plt.title("Error pesado: $E(\Omega) = W(\Omega) \cdot [D(\Omega) - H_R(\Omega)]$")
        plt.xlabel("Frecuencia Normalizada")
        plt.ylabel("Magnitud")
        plt.legend()
    
        H = np.fft.fft(h_coeffs, fft_sz)
    
        plt.figure(2)
        plt.plot(frecuencias[:half_fft_sz], 20*np.log10(np.abs(H[:half_fft_sz])), label=f'Iter: {niter}')

        plt.title("Respuesta en frecuencia de módulo: $ \\left|H(\Omega)\\right| $")
        plt.xlabel("Frecuencia Normalizada")
        plt.ylabel("$\\left|H(\Omega)\\right|_{{dB}}$")
        plt.legend()
    
        plt.figure(3)
        Aw_ext = np.interp(frecuencias[:half_fft_sz], fr_grid, Aw_fr_grid)
        plt.plot(frecuencias[:half_fft_sz], Aw_ext, label=f'$H_{{R{niter}}}$')
        plt.legend()
        plt.show()
        pass
        ## Debug

    return a_coeffs_half, this_error_target, fr_grid[omega_ext_idx]

