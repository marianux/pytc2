 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mariano
"""

import numpy as np
from numbers import Integral, Real

import matplotlib.pyplot as plt

from scipy.signal import find_peaks, welch, csd, freqz

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
    
    order = int(order)
              
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
    total_freq_bins = grid_density * cant_bases
    freq_resolution = 1.0 / total_freq_bins
    # full resolution (fr) fr_grid, desired and wieight arrays
    fr_grid = []
    fr_desired = []
    fr_weight = []
    # indexes of the band-edges corresponding to the fr freq. fr_grid array
    band_edges_idx = []

    min_number_fr_grid = int(np.min( (cant_bases * .1, 20) ))

    for ll in range(nbands):
        
        number_fr_grid = np.max( (min_number_fr_grid, int(np.ceil((band_edges[2 * ll + 1] - band_edges[2 * ll]) / freq_resolution))) )
        
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
        
    order = int(order)
                  
    if not (isinstance(grid_density, (Integral, Real)) and grid_density > 0 ):
        raise ValueError("El argumento 'grid_density' debe ser un número positivo.")
                  
    if not (isinstance(max_iter, (Integral, Real)) and max_iter > 0 ):
        raise ValueError("El argumento 'max_iter' debe ser un número positivo.")

    if not isinstance(debug, bool):
        raise ValueError('debug debe ser un booleano')

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
    total_freq_bins = grid_density * cant_bases
    freq_resolution = 1.0 / total_freq_bins
    # full resolution (fr) fr_grid, desired and wieight arrays
    fr_grid = []
    fr_desired = []
    fr_weight = []
    # indexes of the band-edges corresponding to the fr freq. fr_grid array
    band_edges_idx = []

    min_number_fr_grid = int(np.min( (cant_bases * .1, 20) ))

    for ll in range(nbands):
        
        number_fr_grid = np.max( (min_number_fr_grid, int(np.ceil((band_edges[2 * ll + 1] - band_edges[2 * ll]) / freq_resolution))) )
        
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

class DC_removal_recursive_filter:
    
    def __init__(self, samp_avg = 16, upsample = 1):
        """
        Introducción:
        -------------
        
        Implementación recursiva de un filtro FIR de fase lineal y 
        transferencia:
        
        Tₘₐ(z) = (1 - z^{-D*U})/(D*(1 - z^{-U}))
        
        
        Parámetros:
        -----------
            
        samp_avg : int
            La cantidad de muestras a promediar (D), que para el filtro
            peine, un D más grande redunda en un una mayor selectividad 
            de los "notches" o ceros de transmisión del filtro. 
            (Default = 16 muestras)
        
        upsample : int
            Factor de sobremuestreo (entero U) de la transferencia. Default=1

        Referencias: 
            
        [1] Lyons, R., "Linear-phase DC Removal Filter", dsprelated.com, 
        March 2008, https://www.dsprelated.com/showarticle/58.php
        
        """

        if not (isinstance(samp_avg, (Integral, Real)) and samp_avg > 1 ):
            raise ValueError("El argumento 'samp_avg' debe ser un número positivo.")

        if not (isinstance(upsample, (type(None), Integral, Real)) ):
            raise ValueError("El argumento 'upsample' debe ser un número positivo.")
        
        
        self.samp_avg = samp_avg
        self.upsample = upsample
        self.effective_D = samp_avg * upsample
        
        # se asumen condiciones iniciales nulas
        self.yy_ci = np.zeros(upsample)
        self.xx_ci = np.zeros(self.effective_D)
        
        self.kk_offset = 0

    def reset(self):
        """
        Reseteo del filtro y sus variables internas. Es necesario resetear el 
        filtro cuando se procesa una señal diferente o un bloque disjunto.
        
        """
        # se asumen condiciones iniciales nulas
        self.yy_ci = np.zeros(self.upsample)
        self.xx_ci = np.zeros(self.effective_D)

    def set_initial_conditions(self, xx_ci, yy_ci):
        """
        Configuración de los parámetros iniciales del filtro recursivo.
        
        
        Parámetros:
        -----------
        
        xx_ci  : array_like
            La matriz de datos de entrada para las (samp_avg * upsample)-ésimas
            muestras anteriores a la primer muestra a procesar.
        
        yy_ci  : array_like
            La matriz de datos de entrada para las (upsample)-ésimas muestras 
            anteriores a la primer muestra a procesar.
            
        
        """

        if not isinstance(xx_ci, np.ndarray ):
            raise ValueError("El argumento 'xx_ci' debe ser una un array de numpy.")
        
        if not isinstance(yy_ci, np.ndarray ):
            raise ValueError("El argumento 'yy_ci' debe ser una un array de numpy.")

        # se asumen condiciones iniciales nulas
        self.yy_ci = yy_ci
        self.xx_ci = xx_ci

    def process(self, xx):
        """
        Función que ejecuta la recursión y calcula la salida a partir de la 
        entrada "xx".
        
        Parámetros:
        -----------
        
        xx  : array_like
            La matriz de datos de entrada. Será un vector o array de N muestras.
        
        Retorna:
        --------
        
        yy  : array_like
            La salida del filtro para cada muestra de xx.
       
        """
        
        NN = xx.shape[0]
    
        
        # if xx_ci is None:
        #     xx_ci = self.xx_ci

        # if yy_ci is None:
        #     yy_ci = self.yy_ci

    
        # resultaron ser importante las condiciones iniciales
        yy = np.zeros_like(xx)
        # yy = np.ones_like(xx) * xx[0] * self.effective_D
    
        # para todos los bloques restantes salvo el primero
           
        min_ups = np.min((self.upsample, NN))
        
        for kk in range( min_ups ):
    
            # Calcula la salida según la ecuación recursiva
            yy[kk] = xx[kk] \
                      - self.xx_ci[kk] \
                      + self.yy_ci[kk]

        min_efd = np.min((self.effective_D, NN))
        
        for kk in range(self.upsample, min_efd ):

            # Calcula la salida según la ecuación recursiva
            yy[kk] = xx[kk] \
                      - self.xx_ci[kk] \
                      + yy[(kk - self.upsample)]
    
        #
        # kk += 1
        
        # for kk in range(NN):
        for kk in range(self.effective_D, NN):
    
            # Calcula la salida según la ecuación recursiva
            yy[kk] = xx[kk]  \
                      - xx[kk - self.effective_D] \
                      + yy[kk - self.upsample]
        
            # if(kk == 204):
           #     print(f' yy[{kk}] : {yy[204]}  = xx[{kk}]:{xx[kk]} - xx[{kk - self.effective_D} ]:{xx[kk - self.effective_D]} + yy[{kk - self.upsample} ]:{yy[kk - self.upsample]}')
                
        
        # calculo las condiciones iniciales del siguiente bloque
        
        if self.effective_D < NN:
            xx_ci = xx[(NN - self.effective_D):]
        else:
            xx_ci = np.concatenate((self.xx_ci[:self.effective_D-min_efd], xx[(NN - min_efd):])) 

        if self.upsample < NN:
            yy_ci = yy[(NN - self.upsample):]
        else:
            yy_ci = np.concatenate((self.yy_ci[:self.upsample-min_ups], yy[(NN - min_ups):]))
        

        self.xx_ci = xx_ci.copy()
        self.yy_ci = yy_ci.copy()
    
        # escalo y devuelvo
        return( yy.copy()/self.samp_avg )

    def impulse_response(self, length=None):
        """
        Calcula la respuesta al impulso empírica del filtro
        
        Parámetros:
        -----------
        
        length : int 
            Longitud de la respuesta a calcular, si no se define se 
            calcula automáticamente en función de los parámetros. 
            Default = None
        
        Retorna:
        --------
            
        y      : ndarray 
            Respuesta al impulso
            
        """

        if self.upsample is None:

            raise ValueError("El factor de sobremuestreo 'upsample' debe ser definido." )
            
        if not isinstance(length, (type(None), Integral)):
            raise ValueError("La longitud, puede omitirse o debe ser un número entero" )

        if length is None:

            length = self.upsample * self.samp_avg
            
        impulse = np.zeros(length)
        impulse[0] = 1
        
        self.reset()
        
        return self.process(impulse)
    
    def frequency_response(self, n_freq=None, bTeorica = False):
        """
        Calcula la respuesta en frecuencia teórica/empírica excitando el 
        filtro con ruido blanco
        
        Parámetros:
        -----------
        
        n_freq   : int, None 
            Número de puntos en frecuencia a evaluar.
        
        bTeorica : Bool
            Calcular la respuesta teórica o empírica. 
        
        Retorna:
        --------
            
        w         : ndarray
            Frecuencias normalizadas (0 a π)
        
        frec_resp : ndarray
            Respuesta en frecuencia compleja
            
        """
        
        if not isinstance(n_freq, (type(None), Integral)):
            raise ValueError("La longitud del espectro, puede omitirse o debe ser un número entero" )
        
        if n_freq is None:

            n_freq = 2048
        
        if bTeorica:
        
            # Forma directa (no recursiva) para verificación
            b = np.zeros(self.effective_D + 1)
            b[0] = 1
            b[self.effective_D] = -1
            b = b / self.samp_avg
            
            a = np.zeros(self.upsample + 1)
            a[0] = 1
            a[self.upsample] = -1
            
            w, frec_resp = freqz(b, a, worN=n_freq)
            
        else:
                
            # Generar señal de prueba (ruido blanco)
            # np.random.seed(42)  # Para reproducibilidad
            x = np.random.normal(1, 1, n_freq)
            
            # print(f'Freq Response')
            
            # Procesar señal a través del filtro
            y = self.process(x)
            
            # Eliminar transitorios iniciales
            discard = min(self.samp_avg * self.upsample, n_freq // 10)
            x = x[discard:]
            y = y[discard:]
            
            welch_avg = 10
            # Calcular densidad espectral cruzada y autoespectro
            w, Pxy = csd(x, y, nperseg=n_freq//welch_avg, nfft=n_freq)
            w, Pxx = welch(x, nperseg=n_freq//welch_avg, nfft=n_freq)
            
            # Respuesta en frecuencia estimada
            frec_resp = Pxy / Pxx
        
        return w, frec_resp

class DC_PWL_removal_recursive_filter:
    
    def __init__(self, fs = 2, fpwl = None, samp_avg = 16, cant_ma = 2, upsample = None, batch = None ):
        """
        Introducción:
        -----------
            
        Implementación recursiva de un filtro FIR de fase lineal, para remover:
        
        * DC y baja frecuencia. 
        * Interferencias de frecuencia de línea (power-line) y sus armónicos.
            
        La transferencia del filtro es:
            
        T(z) = z^{-(D-1)*U} - Tₘₐ^{cant_ma}(z)

        siendo
        
        Tₘₐ(z) = (1 - z^{-D*U})/(D*(1 - z^{-U}))
        
        Esta transferencia posee las ventajas de un FIR de fase lineal, y la 
        posibilidad de implementarse de forma recursiva. Debido a la gran 
        cantidad de 0's en su respuesta al impulso, su implementación puede 
        realizarse con muy poco costo computacional, a expensas de su demora.
        Este filtro es en consecuencia, una referencia obligada para el 
        preprocesamiento de señales del ámbito de las bioseñales, como el 
        electrocardiograma (ECG).
        
        Parámetros:
        -----------
            
        fs       : float
            la frecuencia de muestreo. Default = 2 (norm. a Nyquist)
        
        fpwl     : float, None
            la frecuencia de línea de alimentación (típ. 50/60 Hz).
            En caso que no se especifique, el algoritmo intentará la
            detección automática (Default = None).
        
        samp_avg : int
            La cantidad de muestras a promediar (D), que para el filtro
            peine, un D más grande redunda en un una mayor selectividad 
            de los "notches" o ceros de transmisión del filtro. 
            (Default = 16 muestras)
                   
        cant_ma  : int
            La cantiad de promdiadores en cascada que se implementarán. 
            Al igual que D, a mayor cantidad de etapas, más selectividad.
            Debe ser un número par, ya que cada promediador recursivo 
            tiene demora no entera. Default = 2 secciones
        
        upsample : int
            Factor de sobremuestreo (entero U) de la transferencia. Debe 
            ajustarse de forma tal que los notches coincidan con "fpwl".
            Debe verificarse que se cumpla que fpwl = fs / U.  
            En caso que no se especifica, se intenta ajustar de acuerdo 
            a fpwl y a fs. Default=None
        
        batch    : int
            Para registros muy largos, puede especificarse un número de 
            muestras para el procesamiento por bloques. La clase ya prevé la 
            posibilidad de calcular las condiciones iniciales de los bloques 
            adyacentes. Default: None (toda la señal)
        
        Referencias: 
        -----------
            
        [1] Lyons, R., "60-Hz Noise and Baseline Drift Reduction in ECG Signal Processing", 
            dsprelated.com, January 23, 2021, https://www.dsprelated.com/showarticle/1383.php
        
        """

        if not (isinstance(fs, (Integral, Real)) ):
            raise ValueError("La frecuencia de muestreo 'fs' debe ser un número positivo" )
        else:
            if fs > 0: 
                self.fs = fs
            else:
                raise ValueError("La frecuencia de muestreo 'fs' debe ser un número positivo" )
        
        if not (isinstance(fpwl, (type(None), Integral, Real)) ):
            raise ValueError(f"La frecuencia de línea 'fpwl' (typ. 50/60 Hz) debe ser un número positivo menor a fs/2 ({fs/2})" )

        if not (isinstance(samp_avg, (Integral, Real)) and samp_avg > 1 ):
            raise ValueError("El argumento 'samp_avg' debe ser un número positivo.")

        if not (isinstance(upsample, (type(None), Integral, Real)) ):
            raise ValueError("El argumento 'upsample' debe ser un número positivo.")
        
        if not (isinstance(cant_ma, (Integral, Real)) and cant_ma > 1 ):
            raise ValueError("El argumento 'cant_ma' debe ser un número positivo y par.")
        
        if not (isinstance(batch, (type(None), Integral, Real)) ):
            raise ValueError("El argumento 'batch' debe ser un número positivo.")

        
        self.samp_avg = samp_avg
        # forzar un par
        self.cant_ma = (cant_ma//2) * 2
        self.batch = batch

        self.upsample = upsample

        self.fpwl = fpwl
        
        self.t_ma = None

    def reset(self):
        """
        Reseteo del filtro y sus variables internas. Es necesario resetear el 
        filtro cuando se procesa una señal diferente o un bloque disjunto.

        """
            
        self.t_ma = None
        
    def process(self, xx):
        """
        Función que ejecuta la recursión y calcula la salida a partir de la 
        entrada "xx".
        
        
        Parámetros:
        -----------
        
        xx    : array_like
            La matriz de datos de entrada. Puede ser un conjunto de 
            señales, por lo general será una matriz de NxM, siendo N la 
            cantidad de muestras y M la cantidad de señales.
        
        Retorna:
        --------
        
        frec: int
            Frecuencia de línea estimada 50 ó 60 Hz.
       
        """
            
        
        if not isinstance(xx, np.ndarray ):
            raise ValueError("El argumento 'xx_ci' debe ser una un array de numpy.")
        
        # Convertir a array numpy por si acaso es una lista
        xx = np.asarray(xx)
        bFlatreshape = False
        
        # Manejar arrays de 1 dimensión (vectores)
        if xx.ndim == 1:
            # Asumimos que es un vector de N muestras de 1 señal (N, 1)
            xx = xx.reshape(-1, 1)
            bFlatreshape = True
        
        # Manejar arrays de 2 dimensiones
        elif xx.ndim == 2:
            # Si tiene forma (M, N) transponer
            if xx.shape[0] < xx.shape[1]:  # vector fila (1, N) - 1 señal de N muestras
                xx = xx.transpose()
        
        else:
            raise ValueError("El array 'xx' debe ser de 1 o 2 dimensiones.")
        
        NN, cant_sigs = xx.shape
        
        yy = np.zeros_like(xx)


        if self.batch is None:
            
            self.batch = NN
        else:
            
            self.batch = np.clip(self.batch, a_min= 0, a_max=NN )


        if self.fpwl is None:
            
            if self.fs == 2:

                if self.upsample is None:
                    self.fpwl = 1/10
                else:
                    # se setea upsample, adecuo fpwl
                    self.fpwl = self.fs / self.upsample
                    
                print(f"Se asume fpwl = {self.fpwl:3.3f} Hz")
                
            else:
                    
                self.fpwl = self.__detectar_interferencia_linea(xx[:,0], self.fs)
                print(f"Se detectó en fpwl = {self.fpwl} Hz")
            
        else:

            if self.fpwl > self.fs/2:
                raise ValueError(f"La frecuencia de línea 'fpwl' (typ. 50/60 Hz) debe ser un número positivo menor a fs/2 ({self.fs/2})" )

            elif not np.any(self.fpwl == np.array((50, 60))):
                print(f"Se configuró fpwl = {self.fpwl:3.3f} Hz")

        if self.upsample is None:
            
            self.upsample = self.fs // self.fpwl

            if self.upsample * self.fpwl != self.fs:
                print(f"La frecuencia de línea ({self.fpwl} Hz) no es un múltiplo entero de fs ({self.fs} Hz)" )
                print(f"El filtro estará eliminando múltiplos enteros de ({self.fs/self.upsample:.2f} Hz)" )
                warnings.warn(f"Considere remuestrear la señal a {self.fpwl * self.upsample:.2f} Hz, o cualquier múltiplo entero superior a {self.fpwl} Hz para usar este filtro", UserWarning)

        else:

            if self.upsample  < 2:
                
                raise ValueError("El argumento 'upsample' debe ser un número mayor a 1.")

        self.samp_avgdelay = int( self.cant_ma/2*(self.samp_avg - 1) * self.upsample )
        self.demora =        int((self.samp_avg-1)/2*self.cant_ma * self.upsample)

        if self.t_ma is None:
    
            self.t_ma = [[DC_removal_recursive_filter(samp_avg=self.samp_avg, upsample=self.upsample) 
              for _ in range(self.cant_ma)] 
              for _ in range(cant_sigs)]


        # print('Enter')
        
        for ss in range(0, cant_sigs):
        
            # se procesa cada bloque por separado y se concatena la salida
            for jj in range(0, NN, self.batch):
    
                yy_aux = self.t_ma[ss][0].process(xx[jj:jj+self.batch, ss])
        
                yy[jj:jj+self.batch, ss] = yy_aux
        
        
            # print(f' ma: 0 yy[204] : {yy[204]}')
        
            # cascadeamos MA_stages-1 más
            for ii in range(1, self.cant_ma):
        
                # se procesa cada bloque por separado y se concatena la salida
                for jj in range(0, NN, self.batch):
        
                    yy_aux = self.t_ma[ss][ii].process(yy[jj:jj+self.batch, ss] )
            
                    yy[jj:jj+self.batch, ss] = yy_aux
    
                # print(f' ma: {ii} yy[204] : {yy[204]}')

        # Aplicar retardo (D-1)*U
        delayed_x = np.roll(xx, self.samp_avgdelay, axis = 0)
        delayed_x[:self.samp_avgdelay, :] = 0
        
        yy =  delayed_x - yy
        
        if bFlatreshape:
            yy = yy.flatten()
            
        return yy

    def __detectar_interferencia_linea(xx, fs, ancho_banda=4.0, prominencia_min=3.0):
        """
        Detecta interferencia de frecuencia de línea (50/60 Hz) en una señal.
        
        Parámetros:
        -----------
        
        xx : array_like
            Señal de entrada.
        
        fs : float
            Frecuencia de muestreo (en Hz).
        
        ancho_banda : float, opcional
            Ancho de banda alrededor de 50/60 Hz para buscar interferencia (por defecto 4 Hz).
        
        prominencia_min : float, opcional
            Prominencia mínima para considerar un pico significativo (por defecto 3 dB).
        
        Retorna:
        --------
        frec: int
            Frecuencia de línea estimada 50 ó 60 Hz.
        """
        
        # Calculamos el espectro de potencia usando Welch
        freqs, Pxx = welch(xx, fs = fs, nperseg=min(4096, len(xx)), scaling='spectrum')
        Pxx_db = 10 * np.log10(Pxx)  # Convertir a dB
        
        # Buscamos picos en las bandas de 50 Hz y 60 Hz
        # resultados = {'frecuencia_interferencia': None, 'amplitud_pico': None, 'SNR_estimado': None}
        
        snr = np.zeros(2)
        fpwlines = np.array([50, 60])
        
        for ii, f_linea in enumerate(fpwlines):
            # Definimos la banda de búsqueda alrededor de la frecuencia de línea
            banda = (freqs >= f_linea - ancho_banda) & (freqs <= f_linea + ancho_banda)
            
            if not any(banda):
                continue
                
            # Encontramos el pico más prominente en la banda
            peaks, properties = find_peaks(Pxx_db[banda], prominence=prominencia_min)
            
            if len(peaks) > 0:
                idx_max = peaks[np.argmax(properties['prominences'])]
                idx_global = np.where(banda)[0][idx_max]
                
                # Estimamos el SNR (diferencia entre el pico y el percentil 50 del espectro)
                nivel_ruido = np.percentile(Pxx_db[banda], 50)
                snr[ii] = Pxx_db[idx_global] - nivel_ruido
       
        idx_max = np.argmax(snr)
        
        return fpwlines[idx_max]
    
    def impulse_response(self, length=None):
        """
        Calcula la respuesta al impulso empírica del filtro
        
        Parámetros:
        -----------
        
        length : int 
            Longitud de la respuesta a calcular, si no se define se 
            calcula automáticamente en función de los parámetros. 
            Default = None
        
        Retorna:
        --------
            
        y      : ndarray 
            Respuesta al impulso
            
        """

        if self.fpwl is None:
            
            raise ValueError("La frecuencia de línea 'fpwl' (typ. 50/60 Hz) debe ser definida." )

        if self.upsample is None:

            raise ValueError("El factor de sobremuestreo 'upsample' debe ser definido." )
            
        if not isinstance(length, (type(None), Integral)):
            raise ValueError("La longitud, puede omitirse o debe ser un número entero" )

        if length is None:

            length = self.upsample * self.samp_avg
            
            
        if self.batch is None:
            
            self.batch = length

        impulse = np.zeros(length)
        impulse[0] = 1
        
        self.reset()
        
        # print(f'Impulse R')
        
        return self.process(impulse)
    
    def frequency_response(self, n_freq=None, bTeorica = False):
        """
        Calcula la respuesta en frecuencia teórica/empírica excitando el 
        filtro con ruido blanco
        
        Parámetros:
        -----------
        
        n_freq   : int, None 
            Número de puntos en frecuencia a evaluar.
        
        bTeorica : Bool
            Calcular la respuesta teórica o empírica. 
        
        Retorna:
        --------
            
        w         : ndarray
            Frecuencias normalizadas (0 a π)
        
        frec_resp : ndarray
            Respuesta en frecuencia compleja
            
        """
        
        
        if not isinstance(n_freq, (type(None), Integral)):
            raise ValueError("La longitud del espectro, puede omitirse o debe ser un número entero" )
        
        if n_freq is None:

            n_freq = 2048
        
        if bTeorica:
            
            w, h_ma_sq = self.t_ma[0][0].frequency_response(n_freq, bTeorica = True)
            h_delay = np.exp(-1j * w * self.samp_avgdelay)
            
            frec_resp = h_delay - h_ma_sq**self.cant_ma
            
        else:
            
            # # respuesta empírica
            # if self.fpwl is None:
                
            #     raise ValueError("La frecuencia de línea 'fpwl' (typ. 50/60 Hz) debe ser definida." )
    
            # if self.upsample is None:
    
            #     raise ValueError("El factor de sobremuestreo 'upsample' debe ser definido." )
            
            # Eliminar transitorios iniciales
            discard = min(self.samp_avg * self.upsample, n_freq // 10)
                
            # Generar señal de prueba (ruido blanco)
            # np.random.seed(42)  # Para reproducibilidad
            x = np.random.normal(1, 1, n_freq+discard)
            
            # print(f'Freq Response')
            
            # Procesar señal a través del filtro
            y = self.process(x)
            
            # Eliminar transitorios iniciales
            x = x[discard:]
            y = y[discard:]
            
            welch_avg = 10
            # Calcular densidad espectral cruzada y autoespectro
            w, Pxy = csd(x, y, nperseg=n_freq//welch_avg, nfft=(n_freq*2)-1, fs = self.fs)
            w, Pxx = welch(x, nperseg=n_freq//welch_avg, nfft=(n_freq*2)-1, fs = self.fs)
            
            # Respuesta en frecuencia estimada
            frec_resp = Pxy / Pxx
            
        
        return w, frec_resp


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

    if niter == max_iter:
        print(f"No convergió! probar aumentar max_iter {max_iter}")

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

