 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:15:44 2023

@author: mariano
"""

import numpy as np
from numbers import Integral, Real

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import rgb2hex
from collections import defaultdict
from scipy.signal import tf2zpk, TransferFunction, zpk2tf
import sympy as sp

from IPython.display import display, Math

from fractions import Fraction

#%%
   ##########################################
  ## Variables para el análisis simbólico #
 ##########################################
#%%

from .general import s

def tfcascade(tfa, tfb):
    """
    Realiza la cascada de dos funciones de transferencia.

    Esta función toma dos funciones de transferencia, 'tfa' y 'tfb', y calcula la función de transferencia resultante
    de su cascada. La cascada de dos funciones de transferencia es el producto de ambas en el dominio de Laplace.

    Parameters
    -----------
    tfa : TransferFunction
        Primera función de transferencia.
    tfb : TransferFunction
        Segunda función de transferencia.

    Returns
    --------
    TransferFunction
        Función de transferencia resultante de la cascada.

    Raises
    ------
    ValueError
        Si 'tfa' o 'tfb' no son instancias de TransferFunction.

    See Also
    -----------
    :func:`tfadd`
    :func:`pretty_print_lti`


    Examples
    --------
    >>> from scipy.signal import TransferFunction
    >>> tfa = TransferFunction([1, 2], [1, 3, 2])
    >>> tfb = TransferFunction([1], [1, 4])
    >>> tfc = tfcascade(tfa, tfb)
    >>> print(tfc)
    TransferFunction([1, 2], [1, 7, 14, 8])


    Notes
    -----
    Esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.
        
    """
    if not isinstance(tfa, TransferFunction) or not isinstance(tfb, TransferFunction):
        raise ValueError("Los argumentos deben ser instancias de la clase TransferFunction.")

    # Calcula el numerador y el denominador de la función de transferencia resultante
    numerador_resultante = np.polymul(tfa.num, tfb.num)
    denominador_resultante = np.polymul(tfa.den, tfb.den)
    
    # Crea la función de transferencia resultante
    tfc = TransferFunction(numerador_resultante, denominador_resultante)

    return tfc

def tfadd(tfa, tfb):
    """
    Suma dos funciones de transferencia.

    Esta función toma los coeficientes de dos funciones de transferencia, 'tfa' y 'tfb', y devuelve los coeficientes
    de la función de transferencia resultante que es la suma de las dos funciones de transferencia dadas. La función
    resultante se devuelve como un objeto de tipo TransferFunction.

    Parameters
    -----------
    tfa : TransferFunction
        Coeficientes de la primera función de transferencia.
    tfb : TransferFunction
        Coeficientes de la segunda función de transferencia.

    Returns
    --------
    TransferFunction
        Función de transferencia resultante.

    Raises
    ------
    ValueError
        Si 'tfa' o 'tfb' no son instancias de TransferFunction.

    See Also
    -----------
    :func:`tfcascade`
    :func:`pretty_print_lti`

    Examples
    --------
    >>> from scipy.signal import TransferFunction
    >>> tfa = TransferFunction([1, 2], [3, 4])
    >>> tfb = TransferFunction([5, 6], [7, 8])
    >>> tfadd(tfa, tfb)

    
    """

    if not isinstance(tfa, TransferFunction) or not isinstance(tfb, TransferFunction):
        raise ValueError("Los argumentos deben ser instancias de la clase TransferFunction.")

    # Calcula los coeficientes de la función de transferencia resultante
    numerador_resultante = np.polyadd(np.polymul(tfa.num, tfb.den), np.polymul(tfa.den, tfb.num))
    denominador_resultante = np.polymul(tfa.den, tfb.den)

    # Crea la función de transferencia resultante
    tfc = TransferFunction(numerador_resultante, denominador_resultante)
    
    return tfc

def sos2tf_analog(mySOS):
    """
    Convierte una matriz de secciones de segundo orden (SOS) en una función de transferencia analógica.

    Esta función toma una matriz que define las secciones de segundo orden (SOS) del sistema y devuelve la función
    de transferencia analógica resultante. Cada fila de la matriz SOS representa una sección de segundo orden
    con los coeficientes correspondientes.

    Los SOS siempre se definen como::
        
        mySOS= ( [ a1_1 a2_1 a3_1 b1_1 b2_1 b3_1 ] 
                 [ a1_2 a2_2 a3_2 b1_2 b2_2 b3_2 ] 
                 ... 
                 [ a1_N a2_N a3_N b1_N b2_N b3_N ] 
                ) 
    
    donde cada sección o línea de `mySOS` significa matemáticamente
    
    .. math:: T_i = (a_{1i} \, s^2 + a_{2i} \, s + a_{3i})/(b_{1i} \, s^2 + b_{2i} \, s + b_{3i})


    Parameters
    -----------
    mySOS : array_like
        Matriz que define las secciones de segundo orden (SOS) del sistema.


    Returns
    --------
    TransferFunction
        Función de transferencia analógica resultante.


    Raises
    ------
    ValueError
        Si 'mySOS' no es una matriz 2D o si las filas de la matriz no tienen exactamente 6 elementos.


    See Also
    -----------
    :func:`tf2sos_analog`
    :func:`pretty_print_SOS`
    :func:`pretty_print_lti`


    Examples
    --------
    >>> mySOS = np.array([[1, 0.5, 1, 1, 0.2, 1],
    ...                    [1, 1, 1, 1, 1, 1]])
    >>> tf_analog = sos2tf_analog(mySOS)
    >>> print(tf_analog)
    TransferFunctionContinuous(
    array([1., 1., 1.]),
    array([1. , 1.5, 1.7, 1.5, 1.2, 1. ]),
    dt: None
    )



    Nota
    ----
    Esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.
    """
    # Verificar si mySOS es una matriz 2D y si cada fila tiene exactamente 6 elementos
    if not isinstance(mySOS, np.ndarray) or mySOS.ndim != 2 or mySOS.shape[1] != 6:
        print(mySOS)
        raise ValueError("El argumento 'mySOS' debe ser una matriz 2D donde cada fila tiene 6 elementos.")

    # Obtener el número de secciones de segundo orden (SOS) y la forma de la matriz SOS
    SOSnumber, _ = mySOS.shape
    
    # Inicializar los numeradores y denominadores de la función de transferencia
    num = 1
    den = 1
    
    # Iterar sobre cada fila de la matriz SOS
    for ii in range(SOSnumber):
        # Obtener los coeficientes numéricos y denóminos de la sección de segundo orden actual
        sos_num, sos_den = _one_sos2tf(mySOS[ii,:])
        # Multiplicar los coeficientes numéricos y denóminos con los acumulados hasta el momento
        num = np.polymul(num, sos_num)
        den = np.polymul(den, sos_den)

    # Crear la función de transferencia a partir de los coeficientes acumulados
    tf = TransferFunction(num, den)
    
    # Devolver la función de transferencia resultante
    return tf

def tf2sos_analog(num, den=[]):
    """
    Convierte una función de transferencia en forma de coeficientes numéricos y denóminos en una matriz de secciones de segundo orden (SOS) para un sistema analógico.

    Esta función toma los coeficientes numéricos y denóminos de la función de transferencia y devuelve una matriz que
    define las secciones de segundo orden (SOS) del sistema analógico.

    Los SOS siempre se definen como::
        
        mySOS= ( [ a1_1 a2_1 a3_1 b1_1 b2_1 b3_1 ] 
                 [ a1_2 a2_2 a3_2 b1_2 b2_2 b3_2 ] 
                 ... 
                 [ a1_N a2_N a3_N b1_N b2_N b3_N ] 
                ) 
    
    donde cada sección o línea de `mySOS` significa matemáticamente
    
    .. math:: T_i = (a_{1i} \, s^2 + a_{2i} \, s + a_{3i})/(b_{1i} \, s^2 + b_{2i} \, s + b_{3i})


    Parameters
    -----------
    num : array_like, TransferFunction
        Coeficientes numéricos de la función de transferencia.
    den : array_like, opcional
        Coeficientes denóminos de la función de transferencia.


    Returns
    --------
    array_like
        Matriz que define las secciones de segundo orden (SOS) del sistema analógico.


    Raises
    ------
    ValueError
        Si 'num' o 'den' no son instancias de arrays de numpy.


    See Also
    -----------
    :func:`sos2tf_analog`
    :func:`pretty_print_SOS`
    :func:`pretty_print_lti`


    Examples
    --------
    >>> num = [1, 2, 3]
    >>> den = [4, 5, 6]
    >>> sos_analog = tf2sos_analog(num, den)
    >>> print(sos_analog)
    [[1. 2. 3. 4. 5. 6.]]


    Nota
    ----
    Esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.
    """
    
    if not isinstance(num, (list, np.ndarray, TransferFunction)):
        raise ValueError("El argumento 'num' debe ser una lista, instancias de arrays de numpy o un objeto TransferFunction.")
    
    if not isinstance(den, (list,np.ndarray)):
        raise ValueError("El argumento 'den' debe ser una lista o instancias de arrays de numpy.")

    # Convertir la función de transferencia en polos, ceros y ganancia
    if isinstance(num, TransferFunction):
        z, p, k = tf2zpk(num.num, num.den)
    else:
        z, p, k = tf2zpk(num, den)
    
    # Convertir los polos, ceros y ganancia en una matriz de secciones de segundo orden (SOS)
    sos = zpk2sos_analog(z, p, k)

    # Devolver la matriz de secciones de segundo orden (SOS)
    return sos

def zpk2sos_analog(zz, pp, kk):
    """
    Convierte los polos, ceros y ganancia de una función de transferencia en forma de matriz de secciones de segundo orden (SOS) para un sistema analógico.

    Esta función toma los polos, ceros y ganancia de una función de transferencia y devuelve una matriz que define las
    secciones de segundo orden (SOS) del sistema analógico.

    Los SOS siempre se definen como::
        
        mySOS= ( [ a1_1 a2_1 a3_1 b1_1 b2_1 b3_1 ] 
                 [ a1_2 a2_2 a3_2 b1_2 b2_2 b3_2 ] 
                 ... 
                 [ a1_N a2_N a3_N b1_N b2_N b3_N ] 
                ) 
    
    donde cada sección o línea de `mySOS` significa matemáticamente
    
    .. math:: T_i = (a_{1i} \, s^2 + a_{2i} \, s + a_{3i})/(b_{1i} \, s^2 + b_{2i} \, s + b_{3i})
            
    El algoritmo utilizado para convertir de ZPK a formato SOS sigue las sugerencias del libro "Design of Analog Filters" de R. Schaumann, Cap. 5:
        
    1. Asignar ceros a los polos más cercanos.
    2. Ordenar las secciones por Q creciente.
    3. Ordenar las ganancias para maximizar el rango dinámico. Ver [1] cap. 5.
    
    .. [1] Schaumann, Rolf, Haiqiao Xiao, and Van Valkenburg Mac. Design of analog filters 2nd. Edition. Oxford University Press, 2013. ISBN	0195373944, 9780195373943.
    

    Parameters
    -----------
    zz : array_like
        Zeros de la función de transferencia.
    pp : array_like
        Polos de la función de transferencia.
    kk : float
        Ganancia del sistema.

    Returns
    --------
    array_like
        Matriz que define las secciones de segundo orden (SOS) del sistema analógico.

    Raises
    ------
    AssertionError
        Si hay más ceros que polos.
    ValueError
        Si la factorización de la función de transferencia es incorrecta.


    See Also
    -----------
    :func:`sos2tf_analog`
    :func:`pretty_print_SOS`
    :func:`pretty_print_lti`
    :mod:`scipy.signal`


    Examples
    --------
    >>> zz = [1, 2, 3]
    >>> pp = [4, 5, 6]
    >>> kk = 2.5
    >>> sos_analog = zpk2sos_analog(zz, pp, kk)
    >>> print(sos_analog)
    [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]


    Nota
    ----
    Parte de esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.
    """
    if not isinstance(zz, (list, np.ndarray)) or not isinstance(pp, (list, np.ndarray)):
        raise ValueError("Los argumentos 'zz' y 'pp' deben ser listas o arrays.")
    if not isinstance(kk, (int, float)):
        raise ValueError("El argumento 'kk' debe ser un número.")
    
    # Verificar si el filtro está vacío
    if len(zz) == len(pp) == 0:
        return np.array([[0., 0., kk, 1., 0., 0.]])

    # Asegurar que hay más polos que ceros
    if len(zz) > len(pp):
        raise ValueError("El filtro debe tener igual o mayor cantidad de polos que ceros")

    # Calcular el número de secciones SOS
    n_sections = (len(pp) + 1) // 2
    sos = np.zeros((n_sections, 6))

    # Asegurar que los polos y ceros sean pares conjugados
    z = np.concatenate(_cplxreal(zz))
    p = np.concatenate(_cplxreal(pp))

    # Calcular omega_0 y Q para cada polo
    qq = 1 / (2 * np.cos(np.pi - np.angle(p)))

    # Inicializar matrices para polos y ceros de las secciones SOS
    p_sos = np.zeros((n_sections, 2), np.complex128)
    z_sos = np.zeros_like(p_sos)

    # Verificar si hay un cero por sección SOS
    one_z_per_section = n_sections == z.shape[0]

    # Iterar sobre las secciones SOS
    for si in range(n_sections):
        # Seleccionar el polo "peor"
        p1_idx = np.argmax(qq)
        p1 = p[p1_idx]
        p = np.delete(p, p1_idx)
        qq = np.delete(qq, p1_idx)

        # Emparejar ese polo con un cero
        if np.isreal(p1) and np.isreal(p).sum() == 0:
            # Caso especial para establecer una sección de primer orden
            if z.size == 0:
                z1 = np.nan
            else:
                z1_idx = _nearest_real_complex_idx(z, p1, 'real')
                z1 = z[z1_idx]
                z = np.delete(z, z1_idx)
            p2 = z2 = np.nan
        else:
            # Caso SOS
            if z.size == 0:
                z1 = np.nan
            else:
                z1_idx = np.argmin(np.abs(p1 - z))
                z1 = z[z1_idx]
                z = np.delete(z, z1_idx)

            # Determinar los polos y ceros restantes
            if np.isnan(z1):
                z2 = np.nan
                if np.isreal(p1):
                    idx = np.nonzero(np.isreal(p))[0]
                    assert len(idx) > 0
                    p2_idx = idx[np.argmax(qq)]
                    p2 = p[p2_idx]
                    p = np.delete(p, p2_idx)
                else:
                    p2 = p1.conj()
            else:
                if np.isreal(p1):
                    if np.isreal(z1):
                        idx = np.nonzero(np.isreal(p))[0]
                        assert len(idx) > 0
                        p2_idx = idx[np.argmin(np.abs(np.abs(p[idx]) - 1))]
                        p2 = p[p2_idx]
                        assert np.isreal(p2)
                        if one_z_per_section or len(z) == 0:
                            z2 = np.nan
                        else:
                            z2_idx = _nearest_real_complex_idx(z, p2, 'real')
                            z2 = z[z2_idx]
                            assert np.isreal(z2)
                            z = np.delete(z, z2_idx)
                    else:
                        z2 = z1.conj()
                        p2_idx = _nearest_real_complex_idx(p, z1, 'real')
                        p2 = p[p2_idx]
                        assert np.isreal(p2)
                    p = np.delete(p, p2_idx)
                else:
                    p2 = p1.conj()
                    if np.isreal(z1):
                        if one_z_per_section or len(z) == 0:
                            z2 = np.nan
                        else:
                            z2_idx = _nearest_real_complex_idx(z, p1, 'real')
                            z2 = z[z2_idx]
                            assert np.isreal(z2)
                            z = np.delete(z, z2_idx)
                    else:
                        z2 = z1.conj()

        p_sos[si] = [p1, p2]
        z_sos[si] = [z1, z2]

    assert len(p) == 0  # Se han utilizado todos los polos y ceros
    del p, z

    # Construir el sistema, invirtiendo el orden para que los "peores" estén al final
    p_sos = np.reshape(p_sos[::-1], (n_sections, 2))
    z_sos = np.reshape(z_sos[::-1], (n_sections, 2))

    # Asignar ganancias a cada sección SOS
    mmi = np.ones(n_sections)
    gains = np.ones(n_sections, np.array(kk).dtype)
    tf_j = TransferFunction(1.0, 1.0)

    # Calcular ganancias y construir la función de transferencia para cada sección SOS
    for si in range(n_sections):
        this_zz = z_sos[si, np.logical_not(np.isnan(z_sos[si]))]
        this_pp = p_sos[si, np.logical_not(np.isnan(p_sos[si]))]
        num, den = zpk2tf(this_zz, this_pp, 1)
        tf_j = tfcascade(tf_j, TransferFunction(num, den))
        this_zzpp = np.abs(np.concatenate([this_zz, this_pp]))
        this_zzpp = this_zzpp[this_zzpp > 0]
        _, mag, _ = tf_j.bode(np.logspace(np.floor(np.log10(np.min(this_zzpp))) - 2,
                                          np.ceil(np.log10(np.max(this_zzpp))) + 2, 100))
        mmi[si] = 10 ** (np.max(mag) / 20)

    # Calcular la primera ganancia para optimizar el rango dinámico
    gains[0] = kk * (mmi[-1] / mmi[0])

    # Calcular ganancias para cada sección SOS
    for si in range(n_sections):
        if si > 0:
            gains[si] = (mmi[si - 1] / mmi[si])
        num, den = zpk2tf(z_sos[si, np.logical_not(np.isnan(z_sos[si]))],
                           p_sos[si, np.logical_not(np.isnan(p_sos[si]))], gains[si])
        num = np.concatenate((np.zeros(np.max(3 - len(num), 0)), num))
        den = np.concatenate((np.zeros(np.max(3 - len(den), 0)), den))
        sos[si] = np.concatenate((num, den))

    # Verificar la factorización
    tf_verif = sos2tf_analog(sos)
    z_v, p_v, k_v = tf2zpk(tf_verif.num, tf_verif.den)
    num_t, den_t = zpk2tf(zz, pp, kk)

    if np.std(num_t - tf_verif.num) > 1e-10:
        raise ValueError('Factorización incorrecta: los ceros no coinciden')

    if np.std(den_t - tf_verif.den) > 1e-10:
        raise ValueError('Factorización incorrecta: los polos no coinciden')

    return sos

def pretty_print_lti(num, den=None, displaystr=True):
    """
    Genera una representación matemática de una función de transferencia lineal en función de sus coeficientes numéricos.

    Esta función toma los coeficientes del numerador y, opcionalmente, los del denominador de una función de transferencia lineal y genera una representación matemática en función de estos coeficientes. Los parámetros opcionales permiten especificar si se debe mostrar o devolver la cadena formateada.
    

    Parameters
    -----------
    num : array_like, lista o TransferFunction
        Coeficientes del numerador de la función de transferencia.
    den : array_like, opcional
        Coeficientes del denominador de la función de transferencia. Por defecto es None.
    displaystr : bool, opcional
        Indica si mostrar el resultado como salida o devolverlo como una cadena de texto. Por defecto es True.


    Returns
    --------
    None or str
        Si displaystr es True, muestra la cadena formateada, si no, devuelve la cadena.


    Raises
    ------
    ValueError
        Si los coeficientes numéricos no son de tipo array_like, lista, o un objeto TransferFunction.
        Si los coeficientes del denominador son proporcionados pero no son de tipo array_like.
        Si el argumento displaystr no es de tipo bool.


    See Also
    -----------
    :func:`pretty_print_bicuad_omegayq`
    :func:`_build_poly_str`

    Examples
    --------
    >>> num = [1, 2, 3]
    >>> den = [4, 5, 6]
    >>> pretty_print_lti(num, den)
    # Devuelve la cadena formateada de la función de transferencia.
    """
    
    if not isinstance(num, (list, np.ndarray, TransferFunction)):
        raise ValueError("El argumento 'num' debe ser una lista, un array o un objeto TransferFunction.")
    if den is not None and not isinstance(den, (list, np.ndarray)):
        raise ValueError("El argumento 'den' debe ser una lista o un array.")
    if not isinstance(displaystr, bool):
        raise ValueError("El argumento 'displaystr' debe ser de tipo booleano.")

    if den is None and isinstance(num, TransferFunction ):
        this_lti = num
    else:
        this_lti = TransferFunction(num, den)
    
    num_str_aux = _build_poly_str(this_lti.num)
    den_str_aux = _build_poly_str(this_lti.den)

    strout = r'\frac{' + num_str_aux + '}{' + den_str_aux + '}'

    if displaystr:
        display(Math(strout))
    else:
        return strout

def parametrize_sos(num, den):
    '''
    Parametriza una función de transferencia de segundo orden en función de sus coeficientes.


    Parameters
    -----------
    num : Poly
        Coeficientes del numerador.
    den : Poly
        Coeficientes del denominador.


    Returns
    --------
    tuple
        Una tupla que contiene los siguientes elementos:
            num : Poly
                Coeficientes del numerador parametrizado.
            den : Poly
                Coeficientes del denominador parametrizado.
            w_on : Rational
                Frecuencia natural de oscilación.
            Q_n : Rational
                Factor de calidad del numerador.
            w_od : Rational
                Frecuencia natural de oscilación del denominador.
            Q_d : Rational
                Factor de calidad del denominador.
            K : Rational
                Ganancia.

    Raises
    ------
    ValueError
        Si los coeficientes numéricos no son de tipo Poly.
        Si los coeficientes del denominador no son proporcionados o no son de tipo Poly.


    See Also
    -----------
    :func:`pretty_print_bicuad_omegayq`
    :func:`pretty_print_SOS`


    Examples
    --------
    >>> num = sp.Poly((a*s + b),s)
    >>> den = sp.Poly((c*s + d),s)
    >>> sos_bili, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)
    >>> 
    >>> num = sp.Poly((a*s),s)
    >>> sos_bili1, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)
    >>> 
    >>> num = sp.Poly((a),s)
    >>> sos_bili2, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)
    >>> 
    >>> num = sp.Poly((a*s**2 + b*s + c),s)
    >>> den = sp.Poly((d*s**2 + e*s + f),s)
    >>> sos_1, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)
    >>> 
    >>> num = sp.Poly((a*s**2 + c**2),s)
    >>> sos_2, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)
    >>> 
    >>> num = sp.Poly((a*s**2 + s*b),s)
    >>> sos_3, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)
    >>> 
    >>> num = sp.Poly(a,s)
    >>> sos_4, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)
    >>> 
    >>> num = sp.Poly(a*s**2 ,s)
    >>> sos_5, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)
    >>> 
    >>> num = sp.Poly((b*s),s)
    >>> sos_6, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)


    Nota
    ----
    Esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.
    '''    

    w_od = sp.Rational('0')
    Q_d = sp.Rational('0')
    w_on = sp.Rational('0')
    Q_n = sp.Rational('0')
    K = sp.Rational('0')
    
    if not isinstance(num, sp.Poly):
        raise ValueError("El argumento 'num' debe ser un Poly.")
    if not isinstance(den, sp.Poly):
        raise ValueError("El argumento 'den' debe ser un Poly.")

    den_coeffs = den.all_coeffs()
    num_coeffs = num.all_coeffs()

    if len(den_coeffs) == 3:
    # solo se permiten denominadores de segundo orden
        
        w_od = sp.sqrt(den_coeffs[2]/den_coeffs[0])
        
        omega_Q = den_coeffs[1]/den_coeffs[0]
        
        Q_d = sp.simplify(sp.expand(w_od / omega_Q))
        
        k_d = den_coeffs[0]
        
        # parametrización wo-Q
        den  = sp.poly( s**2 + s * sp.Mul(w_od, 1/Q_d, evaluate=False) + w_od**2, s)


        if num.is_monomial:
            
            if num.degree() == 2:
                # pasaaltos
                
                k_n = num_coeffs[0]
                
                num  = sp.poly( s**2, s)

            elif num.degree() == 1:
                # pasabanda
                
                k_n = num_coeffs[0] * Q_d / w_od
                
                # parametrización wo-Q
                num  = sp.poly( s * w_od / Q_d , s)

            else:
                # pasabajos
                
                k_n = num_coeffs[0] / w_od**2
                
                w_on = w_od
                
                Q_n = sp.oo
                    
                num  = sp.poly( w_od**2, s)

                
        else:
        # no monomial
        
            if num.degree() == 2:

                if num_coeffs[1].is_zero:
                    
                    # cero en w_on
                    w_on = sp.sqrt(num_coeffs[2]/num_coeffs[0])

                    k_n = num_coeffs[0]
                    
                    Q_n = sp.oo
                
                    num  = sp.poly( s**2 + w_on**2, s)

                elif num_coeffs[2].is_zero:
                
                    # cero en w=0 y en w_on
                    w_on = num_coeffs[1]/num_coeffs[0]

                    k_n = num_coeffs[0]

                    num = sp.poly( s*( s + w_on), s)
                
                else: 
                    # polinomio completo -> bicuad completo
                    
                    w_on = sp.sqrt(num_coeffs[2]/num_coeffs[0])
                
                    omega_Q = num_coeffs[1]/num_coeffs[0]
                    
                    Q_n = sp.simplify(sp.expand(w_on / omega_Q))
                    
                    k_n = num_coeffs[0]
                    
                    # parametrización wo-Q
                    num  = sp.poly( s**2 + s * sp.Mul(w_on, 1/Q_n, evaluate=False) + w_on**2, s)

            
            else:
                # solo de primer orden
                
                w_on = num_coeffs[1] / num_coeffs[0]
                
                k_n = num_coeffs[0]
                
                num  = sp.poly( s * w_on, s)

        
        K = sp.simplify(sp.expand(k_n / k_d))

    elif len(den_coeffs) == 2:
        # bilineal
        w_od = den_coeffs[1]/den_coeffs[0]
        
        k_d = den_coeffs[0]
        
        # parametrización wo-Q
        den  = sp.poly( s + w_od, s)        
    
        Q_d = sp.nan
        Q_n = sp.nan
    
        if num.is_monomial:
            
            if num.degree() == 1:
                
                k_n = num_coeffs[0]
                
                # parametrización wo-Q
                num = sp.poly( s, s)        

            else:
                                
                k_n = num_coeffs[0] / w_od
                
                w_on = w_od
                
                num  = sp.poly( w_od, s)

                
        else:
        # no monomial
        
            w_on = num_coeffs[1]/num_coeffs[0]
            
            k_n = num_coeffs[0]
            
            # parametrización wo-Q
            num = sp.poly( s + w_on, s)        
    
        K = sp.simplify(sp.expand(k_n / k_d))

    return (num, den, w_on, Q_n, w_od, Q_d, K )

def pretty_print_bicuad_omegayq(num, den=None, displaystr=True):
    
    """
    Genera una representación matemática de un sistema de segundo orden en función de su frecuencia natural (omega) y su factor de calidad (Q).

    Esta función toma los coeficientes del numerador y, opcionalmente, los del denominador de un sistema de segundo orden y genera una representación matemática en función de la frecuencia natural (omega) y el factor de calidad (Q). Los parámetros opcionales permiten especificar si se debe mostrar o devolver la cadena formateada.

    
    Parameters
    -----------
    num : array_like
        Los coeficientes del numerador del sistema de segundo orden.
    den : array_like, opcional
        Los coeficientes del denominador del sistema de segundo orden. Por defecto es None.
    displaystr : bool, opcional
        Indica si mostrar el resultado como salida o devolverlo como una cadena de texto. Por defecto es True.


    Returns
    --------
    None or str
        Si displaystr es True, muestra la cadena formateada, si no, devuelve la cadena.


    Raises
    ------
    ValueError
        Si los coeficientes numéricos no son proporcionados.
        Si los coeficientes numéricos no tienen una longitud de 3 elementos.


    See Also
    -----------
    :func:`pretty_print_SOS`
    :func:`_build_omegayq_str`


    Examples
    --------
    >>> pretty_print_bicuad_omegayq([1, 2, 1], [1, 1, 1])
    # Devuelve la cadena formateada del sistema de segundo orden.

    """

    if not isinstance(num, (list, np.ndarray)):
        raise ValueError("Los coeficientes numéricos deben ser proporcionados.")

    if den is None:
        if len(num) != 6:
            raise ValueError("Los coeficientes de una SOS deben tener una longitud de 6 elementos.")
        this_sos = num.reshape((1, 6))
    else:
        this_sos = np.hstack((
            np.pad(num, (3 - len(num), 0)),
            np.pad(den, (3 - len(den), 0)))
        ).reshape((1, 6))

    num = this_sos[0, :3]
    den = this_sos[0, 3:]

    if np.all(np.abs(num) > 0):
        # Segundo orden completo, parametrización omega y Q
        num_str_aux = _build_omegayq_str(num)
    elif np.all(num[[0, 2]] == 0) and num[1] > 0:
        # Estilo pasa banda: s . k = s . H . omega/Q
        num_str_aux = _build_omegayq_str(num, den=den)
    elif num[1] == 0 and np.all(num[[0, 2]] > 0):
        # Estilo ceros complejos conjugados: s² + omega²
        kk = num[0]
        if kk == 1.0:
            omega = np.sqrt(num[2])
            num_str_aux = r's^2 + {:3.4g}^2'.format(omega)
        else:
            omega = np.sqrt(num[2] / kk)
            num_str_aux = r'{:3.4g}(s^2 + {:3.4g}^2)'.format(kk, omega)
    else:
        num_str_aux = _build_poly_str(num)

    den_str_aux = _build_omegayq_str(den)

    strout = r'\frac{' + num_str_aux + '}{' + den_str_aux + '}'

    if displaystr:
        display(Math(strout))
    else:
        return strout
    
def pretty_print_SOS(mySOS, mode='default', displaystr=True):
    '''
    Imprime de forma "bonita" una expresión que define a un sistema de segundo orden (SOS)
    
    Esta función toma una matriz que define las secciones de segundo orden (SOS) y muestra la representación matemática de la cadena de sistemas de segundo orden. Los parámetros opcionales permiten especificar el modo de impresión y si se debe mostrar o devolver la cadena formateada.    
    
    Los SOS siempre deben definirse como::
        
        mySOS= ( [ a1_1 a2_1 a3_1 b1_1 b2_1 b3_1 ] 
                 [ a1_2 a2_2 a3_2 b1_2 b2_2 b3_2 ] 
                 ... 
                 [ a1_N a2_N a3_N b1_N b2_N b3_N ] 
                ) 
    
    donde cada sección o línea de `mySOS` significa matemáticamente
    
    .. math:: T_i = (a_{1i} \, s^2 + a_{2i} \, s + a_{3i})/(b_{1i} \, s^2 + b_{2i} \, s + b_{3i})
            
    
    Parameters
    -----------
    mySOS : matriz numpy
        La matriz que define los coeficientes de las secciones de segundo orden.
    mode : str, opcional
        El modo de impresión. Puede ser 'default' o 'omegayq'. Por defecto es 'default'.
    displaystr : bool, opcional
        Indica si mostrar el resultado como salida o devolverlo como una cadena de texto. Por defecto es True.


    Returns
    --------
    None or str
        Si displaystr es True, muestra la cadena formateada, si no, devuelve la cadena.


    Raises
    ------
    ValueError
        Si el modo de impresión no es válido.
        Si mySOS no es una matriz numpy.
        Si mySOS no tiene exactamente 6 columnas.
        Si displaystr no es un booleano.


    See Also
    -----------
    :func:`parametrize_sos`
    :func:`pretty_print_lti`
    :func:`pretty_print_bicuad_omegayq`


    Examples
    --------
    >>> mySOS = array([[1, 2, 1, 1, 1, 1], [1, 3, 1, 1, 4, 1]])
    >>> pretty_print_SOS(mySOS)
    [latex output]
    
    '''

    if not isinstance(mySOS, np.ndarray):
        raise ValueError('mySOS debe ser una matriz numpy')

    if mySOS.shape[1] != 6:
        raise ValueError('mySOS debe tener exactamente 6 columnas')

    if not isinstance(displaystr, bool):
        raise ValueError('displaystr debe ser un booleano')

    sos_str = '' 
    
    valid_modes = ['default', 'omegayq']
    if mode not in valid_modes:
        raise ValueError('mode debe ser uno de %s, no %s'
                         % (valid_modes, mode))
    SOSnumber, _ = mySOS.shape
    
    for ii in range(SOSnumber):
        
        if mode == "omegayq" and mySOS[ii, 3] > 0:
            sos_str += r' . ' + pretty_print_bicuad_omegayq(mySOS[ii, :], displaystr=False)
        else:
            num, den = _one_sos2tf(mySOS[ii, :])
            this_tf = TransferFunction(num, den)
            sos_str += r' . ' + pretty_print_lti(this_tf, displaystr=False)

    sos_str = sos_str[2:]

    if displaystr:
        display(Math(r' ' + sos_str))
    else:
        return sos_str

def analyze_sys(all_sys, sys_name=None, img_ext='none', same_figs=True, annotations=True, xaxis='omega', fs=2*np.pi):
    """
    Analiza el comportamiento de un sistema lineal en términos de:

          * Respuesta de magnitud y fase o gráfico de Bode
          * Mapa de polos y ceros
          * Retardo de grupo
          
        La función admite el sistema a analizar (*all_sys*) como:
            
            * uno o una lista de objetos TransferFunction
            * una matriz que define varias secciones de segundo orden (SOS).
            
        Si *all_sys* es una matriz SOS, la función muestra cada una de las SOS y el sistema resultante de la cascada de todas las SOS.
    
    Esta función toma un sistema lineal (ya sea una lista de objetos TransferFunction o una matriz que define una cascada de secciones de segundo orden) y realiza un análisis completo del comportamiento del sistema, incluyendo trazado de gráficos de Bode, mapa de polos y ceros, y gráfico de retardo de grupo. Los parámetros opcionales permiten personalizar el análisis según las necesidades del usuario.
    
    
    Parameters
    -----------
    all_sys : lista o matriz (Nx6)
        El sistema lineal a analizar. Ya sea una lista de objetos TransferFuncion de scipy.signal o una matriz que define una cascada de SOS.
    sys_name : string o lista, opcional
        Las etiquetas o descripción del sistema. Por defecto es None.
    img_ext : string ['none', 'png', 'svg'], opcional
        Cuando es diferente de 'none', la función guarda los resultados del gráfico en un archivo con la extensión indicada. Por defecto es 'none'.
    same_figs : booleano, opcional
        Usa siempre los mismos números de figura para trazar resultados. Cuando es False, cada llamada produce un nuevo grupo de figuras en un contenedor de gráficos separado. Por defecto es True.
    annotations : booleano, opcional
        Agrega anotaciones al gráfico del mapa de polos y ceros. Cuando es True, cada singularidad estará acompañada del valor de omega (es decir, la distancia radial al origen) y Q (es decir, una medida de proximidad al eje jw). Por defecto es True.
    xaxis : string, opcional ['omega', 'freq', 'norm']
        El significado del eje X: "omega" se mide en radianes/s y se prefiere para sistemas analógicos. "freq" se mide en Hz (1/s) y es válido tanto para sistemas digitales como analógicos. "norm" es una versión normalizada con la norma definida por fs. Por defecto es 'omega'.
    fs : valor real, opcional
        La frecuencia de muestreo del sistema digital o la norma para xaxis igual a "norm". Solo es válido si digital es True. Por defecto es None (definido en 1/dlti.dt).
    
    
    Raises
    ------
    ValueError
        Si la extensión de imagen no es válida.
        Si sys_name no es una lista o un string.
        Si all_sys no es una lista o una matriz.
        Si xaxis no es válido.


    Returns
    --------
    return_values : lista
        Lista con tres pares de manijas de figuras y ejes de cada gráfico mostrado.


    See Also
    -----------
    :func:`pretty_print_bicuad_omegayq`
    :func:`bodePlot`
    :func:`pzmap`


    Examples
    --------
    >>> # Analiza un sistema con w0 = 1 rad/s y Q = sqrt(2)/2
    >>> import numpy as np
    >>> from scipy import signal as sig
    >>> from pytc2.sistemas_lineales import analyze_sys, pretty_print_bicuad_omegayq
    >>> Q = np.sqrt(2)/2
    >>> w0 = 1
    >>> num = np.array([w0**2])
    >>> den = np.array([1., w0 / Q, w0**2])
    >>> H1 = sig.TransferFunction(num, den)
    >>> pretty_print_bicuad_omegayq(num, den)
    >>> analyze_sys(H1, sys_name='mi Example')
    >>> 
    >>> # Compara el sistema anterior con otros dos con valores diferentes de Q
    >>> Q = 5
    >>> w0 = 1
    >>> num = np.array([w0**2])
    >>> den = np.array([1., w0 / Q, w0**2])
    >>> H2 = sig.TransferFunction(num, den)
    >>> analyze_sys([H1, H2], sys_name=['H1', 'H2'])


    """

    valid_ext = ['none', 'png', 'svg']
    if img_ext not in valid_ext:
        raise ValueError('La extensión de imagen debe ser una de %s, no %s'
                         % (valid_ext, img_ext))

    if isinstance(all_sys, np.ndarray):
        if all_sys.shape[1] != 6:
            raise ValueError('La matriz all_sys debe tener 6 columnas')
        cant_sys = 1
        all_sys = [all_sys]
    elif isinstance(all_sys, list):
        cant_sys = len(all_sys)
    else:
        raise ValueError('all_sys debe ser una lista o una matriz numpy')

    if sys_name is None:
        sys_name = [str(ii) for ii in range(cant_sys)]
    elif not isinstance(sys_name, list) and not isinstance(sys_name, str):
        raise ValueError('sys_name debe ser una lista o un string')
    
    if isinstance(sys_name, str):
        sys_name = [sys_name]
        
    if isinstance(sys_name, list) and len(sys_name) != cant_sys:
        raise ValueError('sys_name debe tener igual cantidad de etiquetas que ')
        
        
    # Check valid values for xaxis
    valid_xaxis = ['omega', 'freq', 'norm']
    if xaxis not in valid_xaxis:
        raise ValueError('El valor de xaxis debe ser uno de %s, no %s'
                         % (valid_xaxis, xaxis))
    
    # Check valid type for same_figs
    if not isinstance(same_figs, bool):
        raise ValueError('same_figs debe ser un booleano')

    # Check valid type for annotations
    if not isinstance(annotations, bool):
        raise ValueError('annotations debe ser un booleano')

    # Check valid type for fs
    if fs is not None and not isinstance(fs, (int, float)):
        raise ValueError('fs debe ser None o un valor real')

    # Gráficos de BODE
    return_values = []
        
    if same_figs:
        fig_id = 1
    else:
        fig_id = 'none'
    axes_hdl = []

    for ii in range(cant_sys):
        
        if isinstance(all_sys, list) and isinstance(all_sys[ii], TransferFunction):
        
            if all_sys[ii].dt is None:
                this_digital = False
            else:
                this_digital = True

        else:
            this_digital = False

        fig_id, axes_hdl = bodePlot(all_sys[ii], fig_id, axes_hdl, filter_description=sys_name[ii], digital=this_digital, xaxis=xaxis, fs=fs)


    if img_ext != 'none':
        plt.savefig('_'.join(sys_name) + '_Bode.' + img_ext, format=img_ext)

    return_values += [[fig_id, axes_hdl]]
    
    # fig_id = 6
    # axes_hdl = ()

    # for ii in range(cant_sys):
    #     fig_id, axes_hdl = bodePlot(all_sys[ii], fig_id, axes_hdl, filter_description=sys_name[ii])

    # axes_hdl[0].set_ylim(bottom=-3)

    # if img_ext != 'none':
    #     plt.savefig('_'.join(sys_name) + '_Bode-3db.' + img_ext, format=img_ext)


    # Mapas de polos y ceros
    
    if same_figs:
        analog_fig_id = 2
        digital_fig_id = 3
    else:
        analog_fig_id = 'none'
        digital_fig_id = 'none'
    
    analog_axes_hdl = []
    digital_axes_hdl = []
    
    for ii in range(cant_sys):
    
        if isinstance(all_sys[ii], np.ndarray):
            
            thisFilter = sos2tf_analog(all_sys[ii])

            analog_fig_id, analog_axes_hdl = pzmap(thisFilter, filter_description=sys_name[ii], fig_id=analog_fig_id, axes_hdl=analog_axes_hdl, annotations=annotations, digital=this_digital, fs=fs)
            
        else:
                
            if all_sys[ii].dt is None:
                analog_fig_id, analog_axes_hdl = pzmap(all_sys[ii], filter_description=sys_name[ii], fig_id=analog_fig_id, axes_hdl=analog_axes_hdl, annotations=annotations)
                
            else:
                digital_fig_id, digital_axes_hdl = pzmap(all_sys[ii], filter_description=sys_name[ii], fig_id=digital_fig_id, axes_hdl=digital_axes_hdl, annotations=annotations)


    return_values += [[analog_fig_id, analog_axes_hdl]]
        
    return_values += [[digital_fig_id, digital_axes_hdl]]


    if isinstance(all_sys[ii], np.ndarray) or (isinstance(all_sys[ii], TransferFunction) and all_sys[ii].dt is None):
        analog_axes_hdl.legend()
        if img_ext != 'none':
            plt.figure(analog_fig_id)
            plt.savefig('_'.join(sys_name) + '_Analog_PZmap.' + img_ext, format=img_ext)
    else:
        digital_axes_hdl.legend()
        if img_ext != 'none':
            plt.figure(digital_fig_id)
            plt.savefig('_'.join(sys_name) + '_Digital_PZmap.' + img_ext, format=img_ext)

    
    # Gráficos de retardo de grupo
    
    if same_figs:
        fig_id = 4
    else:
        fig_id = 'none'
    
    for ii in range(cant_sys):
        
        if isinstance(all_sys, list) and isinstance(all_sys[ii], TransferFunction) and all_sys[ii].dt is None:
            this_digital = False
        else:
            this_digital = True
        
        fig_id, axes_hdl = GroupDelay(all_sys[ii], fig_id, filter_description=sys_name[ii], digital=this_digital, xaxis=xaxis, fs=fs)
    
    return_values += [[fig_id, axes_hdl]]
    
    # axes_hdl.legend(sys_name)

    # axes_hdl.set_ylim(bottom=0)

    if img_ext != 'none':
        plt.savefig('_'.join(sys_name) + '_GroupDelay.' + img_ext, format=img_ext)


    return return_values

def pzmap(myFilter, annotations=False, filter_description=None, fig_id='none', axes_hdl='none', digital=False, fs=2*np.pi):
    """
    Grafica el mapa de polos y ceros de un filtro dado.

    Parameters
    -----------
    myFilter : LTI object
        Objeto del filtro.
    annotations : bool, opcional
        Indica si se deben añadir anotaciones a los polos y ceros. 
        El valor predeterminado es False.
    filter_description : str, opcional
        Descripción del filtro. El valor predeterminado es None.
    fig_id : str or int, opcional
        Identificador de la figura. Si se establece en 'none', se creará una nueva figura.
        El valor predeterminado es 'none'.
    axes_hdl : str or axes handle, opcional
        Identificador o handle  del eje. Si se establece en 'none', se utilizará el eje actual.
        El valor predeterminado es 'none'.
    digital : bool, opcional
        Indica si el filtro es digital. El valor predeterminado es False.
    fs : float, opcional
        Frecuencia de muestreo. El valor predeterminado es 2*pi.


    Returns
    --------
    fig_id : int
        Identificador de la figura creada.
    axes_hdl : axes handle
        handle  del eje utilizado para el gráfico.


    Raises
    ------
    ValueError
        Si `fig_id` no es un string o un entero.
        Si `axes_hdl` no es un string o una handle  de eje válida.
        Si `digital` no es un booleano.
        Si `fs` no es un valor numérico.


    See Also
    -----------
    :func:`analyze_sys`
    :func:`bodePlot`
    :func:`pzmap`


    Examples
    --------
    >>> fig_id, ax_hdl = pzmap(myFilter, annotations=True, filter_description='Filtro Pasabajos')
    
    
    """
    
    # Chequeo de argumentos
    if not isinstance(fig_id, (str, int)):
        raise ValueError('fig_id debe ser un string o un entero.')

    if not isinstance(axes_hdl, (str, list, plt.Axes)):
        raise ValueError('axes_hdl debe ser un string o un handle  de eje válida.')

    if not isinstance(digital, bool):
        raise ValueError('digital debe ser un booleano.')

    if not isinstance(fs, (int, float)):
        raise ValueError('fs debe ser un valor numérico.')

    # Configuración de la figura y el eje
    if fig_id == 'none':
        fig_hdl = plt.figure()
        fig_id = fig_hdl.number
    else:
        if plt.fignum_exists(fig_id):
            fig_hdl = plt.figure(fig_id)
        else:
            fig_hdl = plt.figure(fig_id)
            fig_id = fig_hdl.number

    if not isinstance(axes_hdl, plt.Axes):
        axes_hdl = plt.gca()
        
    # Obtener los polos y ceros del filtro
    z, p, k = tf2zpk(myFilter.num, myFilter.den)

    # Añadir círculo unitario y ejes de cero
    unit_circle = patches.Circle((0, 0), radius=1, fill=False,
                                 color='gray', ls='dotted', lw=2)
    axes_hdl.add_patch(unit_circle)
    plt.axvline(0, color='0.7')
    plt.axhline(0, color='0.7')

    # Añadir líneas de círculo
    # maxRadius = np.abs(10*np.sqrt(p[0]))

    # Graficar los polos y configurar las propiedades del marcador
    if filter_description is None:
        poles = plt.plot(p.real, p.imag, 'x', markersize=9)
    else:
        poles = plt.plot(p.real, p.imag, 'x', markersize=9, label=filter_description)

    # Graficar los ceros y configurar las propiedades del marcador
    zeros = plt.plot(z.real, z.imag, 'o', markersize=9,
                     color='none',
                     markeredgecolor=poles[0].get_color(),  # mismo color que los polos
                     markerfacecolor='white'
                     )

    # agregar información a los polos y ceros
    # primero con polos
    w0, aux_idx = np.unique(np.abs(p), return_index=True)
    qq = 1 / (2 * np.cos(np.pi - np.angle(p[aux_idx])))

    # distancia de etiqueta al punto de datos
    lab_mod = 40

    # alternancia de signo para singularidades complejas conjugadas
    aux_sign = np.sign(np.random.uniform(-1, 1))

    for ii in range(len(w0)):

        rand_dir = np.random.uniform(0, 2*np.pi)

        xy_coorde = (lab_mod * np.cos(rand_dir), lab_mod * np.sin(rand_dir))

        if(xy_coorde[0] < 0.0):
            halign = 'left'
        else:
            halign = 'right'

        if(xy_coorde[1] < 0.0):
            valign = 'top'
        else:
            valign = 'bottom'

        if p[aux_idx[ii]].imag > 0.0:
            # anotar solo singularidades complejas conjugadas Q
            aux_sign = aux_sign * -1

            circle = patches.Circle((0, 0), radius=w0[ii], color=poles[0].get_color(), fill=False, ls=(0, (1, 10)), lw=0.7)
            axes_hdl.add_patch(circle)
            plt.axvline(0, color='0.7')
            plt.axhline(0, color='0.7')

            if annotations:
                axes_hdl.annotate('$\omega$ = {:3.3g} \n Q = {:3.3g}'.format(w0[ii], qq[ii]),
                                   xy=(p[aux_idx[ii]].real, p[aux_idx[ii]].imag * aux_sign), xycoords='data',
                                   xytext=xy_coorde, textcoords='offset points',
                                   arrowprops=dict(facecolor=poles[0].get_color(), shrink=0.15,
                                                   width=1, headwidth=5),
                                   horizontalalignment=halign, verticalalignment=valign,
                                   color=poles[0].get_color(),
                                   bbox=dict(edgecolor=poles[0].get_color(), facecolor=_complementaryColor(rgb2hex(poles[0].get_color())), alpha=0.4))

        else:
            # anotar con singularidades omega real
            if annotations:
                axes_hdl.annotate('$\omega$ = {:3.3g}'.format(w0[ii]),
                                  xy=(p[aux_idx[ii]].real, p[aux_idx[ii]].imag), xycoords='data',
                                  xytext=xy_coorde, textcoords='offset points',
                                  arrowprops=dict(facecolor=poles[0].get_color(), shrink=0.15,
                                                  width=1, headwidth=5),
                                  horizontalalignment=halign, verticalalignment=valign,
                                  color=poles[0].get_color(),
                                  bbox=dict(edgecolor=poles[0].get_color(), facecolor=_complementaryColor(rgb2hex(poles[0].get_color())), alpha=0.4))

    # y luego con los ceros
    w0, aux_idx = np.unique(np.abs(z), return_index=True)
    qq = 1 / (2 * np.cos(np.pi - np.angle(z[aux_idx])))

    # alternancia de signo para singularidades complejas conjugadas
    aux_sign = np.sign(np.random.uniform(-1, 1))

    for ii in range(len(w0)):

        aux_sign = aux_sign * -1

        rand_dir = np.random.uniform(0, 2*np.pi)

        xy_coorde = (lab_mod * np.cos(rand_dir), lab_mod * np.sin(rand_dir))

        if(xy_coorde[0] < 0.0):
            halign = 'left'
        else:
            halign = 'right'

        if(xy_coorde[1] < 0.0):
            valign = 'top'
        else:
            valign = 'bottom'

        if z[aux_idx[ii]].imag > 0.0:

            circle = patches.Circle((0, 0), radius=w0[ii], color=poles[0].get_color(), fill=False, ls=(0, (1, 10)), lw=0.7)
            axes_hdl.add_patch(circle)
            plt.axvline(0, color='0.7')
            plt.axhline(0, color='0.7')

            if annotations:
                axes_hdl.annotate('$\omega$ = {:3.3g} \n Q = {:3.3g}'.format(w0[ii], qq[ii]),
                                   xy=(z[aux_idx[ii]].real, z[aux_idx[ii]].imag * aux_sign), xycoords='data',
                                   xytext=xy_coorde, textcoords='offset points',
                                   arrowprops=dict(facecolor=poles[0].get_color(), shrink=0.15,
                                                   width=1, headwidth=5),
                                   horizontalalignment=halign, verticalalignment=valign,
                                   color=poles[0].get_color(),
                                   bbox=dict(edgecolor=poles[0].get_color(), facecolor=None, alpha=0.4))

        else:
            # anotar con singularidades omega real
            if annotations:
                axes_hdl.annotate('$\omega$ = {:3.3g}'.format(w0[ii]),
                                  xy=(z[aux_idx[ii]].real, z[aux_idx[ii]].imag), xycoords='data',
                                  xytext=xy_coorde, textcoords='offset points',
                                  arrowprops=dict(facecolor=poles[0].get_color(), shrink=0.15,
                                                  width=1, headwidth=5),
                                  horizontalalignment=halign, verticalalignment=valign,
                                  color=poles[0].get_color(),
                                  bbox=dict(edgecolor=poles[0].get_color(), facecolor=None, alpha=0.4))

    # Escalar ejes para que quepan
    r_old = axes_hdl.get_ylim()[1]

    r = 1.1 * np.amax(np.concatenate(([r_old/1.1], abs(z), abs(p), [1])))
    plt.axis('scaled')
    plt.axis([-r, r, -r, r])

    # Encontrar duplicados por mismas coordenadas de píxeles
    poles_xy = axes_hdl.transData.transform(np.vstack(poles[0].get_data()).T)
    zeros_xy = axes_hdl.transData.transform(np.vstack(zeros[0].get_data()).T)

    d = defaultdict(int)
    coords = defaultdict(tuple)
    for xy in poles_xy:
        key = tuple(np.rint(xy).astype('int'))
        d[key] += 1
        coords[key] = xy
    for key, value in d.items():
        if value > 1:
            x, y = axes_hdl.transData.inverted().transform(coords[key])
            plt.text(x, y,
                     r' ${}^{' + str(value) + '}$',
                     fontsize=13,
                     )

    d = defaultdict(int)
    coords = defaultdict(tuple)
    for xy in zeros_xy:
        key = tuple(np.rint(xy).astype('int'))
        d[key] += 1
        coords[key] = xy
    for key, value in d.items():
        if value > 1:
            x, y = axes_hdl.transData.inverted().transform(coords[key])
            plt.text(x, y,
                     r' ${}^{' + str(value) + '}$',
                     fontsize=13,
                     )

    # Configuraciones adicionales según si el filtro es digital o analógico
    if myFilter.dt is None:
        digital = False
    else:
        digital = True

    if digital:
        plt.xlabel(r'$\Re(z)$')
        plt.ylabel(r'$\Im(z)$')
    else:
        plt.xlabel(r'$\sigma$')
        plt.ylabel('j'+r'$\omega$')

    plt.grid(True, color='0.9', linestyle='-', which='both', axis='both')

    fig_hdl.suptitle('Mapa de Polos y Ceros')

    if not(filter_description is None):
        axes_hdl.legend()

    return fig_id, axes_hdl

def group_delay(freq, phase):
    """
    Calcula el retardo de grupo para una función de fase.


    Parameters
    -----------
    freq : array_like
        La grilla de frecuencia a la que se calcula la fase.
    phase : array_like
        La fase de la función para la cual se calcula el retardo de grupo.



    Returns
    --------
    gd : array_like
        Estimación del retardo de grupo, que es la derivada de la fase
        respecto a la frecuencia cambiada de signo.


    Raises
    ------
    ValueError
        Si `freq` y `phase` no tienen la misma longitud.
        Si `freq` o `phase` no son arreglos NumPy.


    Examples
    --------
    >>> freq = np.linspace(0, 10, 100)
    >>> phase = np.sin(freq)
    >>> group_delay(freq, phase)
    
    """

    # Chequeo de argumentos
    if not isinstance(freq, np.ndarray) or not isinstance(phase, np.ndarray):
        raise ValueError("freq y phase deben ser arreglos NumPy.")
    if len(freq) != len(phase):
        raise ValueError("freq y phase deben tener la misma longitud.")

    # Calcular la derivada de la fase respecto a la frecuencia
    groupDelay = -np.diff(phase) / np.diff(freq)

    # Agregar el último valor para que tenga la misma longitud que el arreglo original
    return np.append(groupDelay, groupDelay[-1])

def GroupDelay(myFilter, fig_id='none', filter_description=None, npoints=1000, digital=False, xaxis='omega', fs=2*np.pi):
    """
    Calcula y grafica el retardo de grupo de un filtro.


    Parameters
    -----------
    myFilter : array_like o scipy.signal.TransferFunction
        Coeficientes del filtro o objeto TransferFunction del filtro.
    fig_id : str o int, opcional
        Identificador de la figura. Si es 'none', crea una nueva figura. Por defecto es 'none'.
    filter_description : str, opcional
        Descripción del filtro. Por defecto es None.
    npoints : int, opcional
        Número de puntos para muestrear el eje de frecuencia. Por defecto es 1000.
    digital : bool, opcional
        Indicador de si el filtro es digital. Por defecto es False.
    xaxis : str, opcional
        Tipo de eje x ('omega', 'freq', 'norm'). Por defecto es 'omega'.
    fs : float, opcional
        Frecuencia de muestreo. Por defecto es 2*pi.


    Returns
    --------
    fig_id : int
        Identificador de la figura.
    axes_hdl : Axes
        Manejador de ejes de la figura.


    Raises
    ------
    ValueError
        Si myFilter no es un array NumPy ni un objeto TransferFunction.
        Si fig_id no es de tipo str, int o None.
        Si npoints no es un entero.
        Si digital no es un booleano.
        Si xaxis no es uno de los valores permitidos: 'omega', 'freq', 'norm'.
        Si fs no es un número.

    
    See Also
    -----------
    :func:`analyze_sys`
    :func:`bodePlot`
    :func:`pzmap`


    Example
    --------
    >>> fig_id, axes_hdl = GroupDelay(myFilter, fig_id=1, filter_description='Filtro pasa bajos', npoints=1000, digital=False, xaxis='omega', fs=2*np.pi)
    
    """

    # Verificar si myFilter es un array NumPy o un objeto TransferFunction
    if not isinstance(myFilter, (np.ndarray, TransferFunction)):
        raise ValueError("myFilter debe ser un array NumPy o un objeto TransferFunction.")

    # Verificar si fig_id es None, str o int
    if not isinstance(fig_id, (type(None), str, int)):
        raise ValueError("fig_id debe ser de tipo str, int o None.")

    # Verificar si npoints es un entero
    if not isinstance(npoints, int):
        raise ValueError("npoints debe ser un entero.")

    # Verificar si digital es un booleano
    if not isinstance(digital, bool):
        raise ValueError("digital debe ser un booleano.")

    # Verificar si xaxis es uno de los valores permitidos
    if xaxis not in ['omega', 'freq', 'norm']:
        raise ValueError("xaxis debe ser uno de los siguientes valores: 'omega', 'freq', 'norm'.")

    # Verificar si fs es un número
    if not isinstance(fs, (int, float)):
        raise ValueError("fs debe ser un número.")


    if isinstance(myFilter, np.ndarray):
        # Sección SOS
        cant_sos = myFilter.shape[0]
        phase = np.empty((npoints, cant_sos+1))
        sos_label = []

        # Calcular la respuesta de fase para cada sección SOS y el filtro completo
        for ii in range(cant_sos):
            num, den = _one_sos2tf(myFilter[ii, :])
            thisFilter = TransferFunction(num, den)

            this_zzpp = np.abs(np.concatenate([thisFilter.zeros, thisFilter.poles]))
            this_zzpp = this_zzpp[this_zzpp > 0]

            if digital:
                w, _, phase[:, ii] = thisFilter.bode(np.linspace(0, np.pi, npoints))
            else:
                w, _, phase[:, ii] = thisFilter.bode(np.logspace(np.floor(np.log10(np.min(this_zzpp))) - 1, np.ceil(np.log10(np.max(this_zzpp))) + 1, npoints))

            sos_label += [filter_description + ' - SOS {:d}'.format(ii)]

        # Filtro completo
        thisFilter = sos2tf_analog(myFilter)

        this_zzpp = np.abs(np.concatenate([thisFilter.zeros, thisFilter.poles]))
        this_zzpp = this_zzpp[this_zzpp > 0]

        if digital:
            w, _, phase[:, cant_sos] = thisFilter.bode(np.linspace(0, np.pi, npoints))
        else:
            w, _, phase[:, cant_sos] = thisFilter.bode(np.logspace(np.floor(np.log10(np.min(this_zzpp))) - 1, np.ceil(np.log10(np.max(this_zzpp))) + 1, npoints))

        sos_label += [filter_description]

        filter_description = sos_label

        phaseRad = phase * np.pi / 180.0

        phaseRad = phaseRad.reshape((npoints, 1+cant_sos))

        # Filtrar huecos y saltos en la respuesta de fase
        all_jump_x, all_jump_y = (np.abs(np.diff(phaseRad, axis=0)) > 4/5*np.pi).nonzero()

        for this_jump_x, this_jump_y in zip(all_jump_x, all_jump_y):
            phaseRad[this_jump_x+1:, this_jump_y] = phaseRad[this_jump_x+1:, this_jump_y] - np.pi

    else:
        # Objeto LTI
        cant_sos = 0

        this_zzpp = np.abs(np.concatenate([myFilter.zeros, myFilter.poles]))
        this_zzpp = this_zzpp[this_zzpp > 0]

        if digital:
            w, _, phase = myFilter.bode(np.linspace(0, np.pi, npoints))
        else:
            w, _, phase = myFilter.bode(np.logspace(np.floor(np.log10(np.min(this_zzpp))) - 1, np.ceil(np.log10(np.max(this_zzpp))) + 1, npoints))

        phaseRad = phase * np.pi / 180.0

        phaseRad = phaseRad.reshape((npoints, 1))

        # Filtrar huecos y saltos en la respuesta de fase
        all_jump = np.where(np.abs(np.diff(phaseRad, axis=0)) > 4/5*np.pi)[0]

        for this_jump_x in all_jump:
            phaseRad[this_jump_x+1:] = phaseRad[this_jump_x+1] - np.pi

    # Calcular el retardo de grupo
    groupDelay = -np.diff(phaseRad, axis=0) / np.diff(w).reshape((npoints-1, 1))

    # Convertir frecuencia a Hz si se solicita
    if xaxis == "freq":
        ww = w / 2 / np.pi
    elif xaxis == "norm":
        if fs is None:
            # Normalizar cada respuesta a su propio Nyquist
            wnorm = 2 * np.pi / myFilter.dt / 2
        else:
            # Normalizado a fs
            wnorm = 2 * np.pi * fs
        ww = w / wnorm
    else:
        ww = w

    # Crear o recuperar figura
    if fig_id == 'none':
        fig_hdl = plt.figure()
        fig_id = fig_hdl.number
    else:
        if plt.fignum_exists(fig_id):
            fig_hdl = plt.figure(fig_id)
        else:
            fig_hdl = plt.figure(fig_id)
            fig_id = fig_hdl.number

    # Graficar el retardo de grupo
    if digital:
        aux_hdl = plt.plot(ww[1:], groupDelay, label=filter_description)    # Gráfico de retardo de grupo
    else:
        aux_hdl = plt.semilogx(ww[1:], groupDelay, label=filter_description)    # Gráfico de retardo de grupo

    # Distinguir la respuesta SOS de la totalidad de la respuesta
    if cant_sos > 0:
        [aa.set_linestyle(':') for aa in aux_hdl[:-1]]
        aux_hdl[-1].set_linewidth(2)

    plt.grid(True)

    # Etiquetas y título del gráfico
    if xaxis == "freq":
        plt.xlabel('Frecuencia [Hz]')
    elif xaxis == "norm":
        plt.gca().set_xlim([0, 1])
        if fs is None:
            this_fs = 1 / myFilter.dt
        else:
            this_fs = fs
        plt.xlabel('Frecuencia normalizada a fs={:3.3f} [#]'.format(this_fs))
    else:
        plt.xlabel('Frecuencia angular [rad/seg]')

    plt.ylabel('Retardo de grupo [seg]')
    plt.title('Retardo de grupo')

    axes_hdl = plt.gca()

    # Mostrar la leyenda si hay descripción del filtro
    if not(filter_description is None):
        axes_hdl.legend()

    return fig_id, axes_hdl

def bodePlot(myFilter, fig_id='none', axes_hdl='none', filter_description=None, npoints=1000, digital=False, xaxis='omega', fs=2*np.pi):
    """
    Grafica el diagrama de Bode (magnitud y fase) de un filtro.


    Parameters
    -----------
    myFilter : array_like o scipy.signal.TransferFunction
        Coeficientes del filtro o objeto TransferFunction del filtro.
    fig_id : str o int, opcional
        Identificador de la figura. Si es 'none', crea una nueva figura. Por defecto es 'none'.
    axes_hdl : str o array_like de Axes, opcional
        Manejador de ejes de la figura. Si es 'none', crea nuevos ejes. Por defecto es 'none'.
    filter_description : str, opcional
        Descripción del filtro. Por defecto es None.
    npoints : int, opcional
        Número de puntos para muestrear el eje de frecuencia. Por defecto es 1000.
    digital : bool, opcional
        Indicador de si el filtro es digital. Por defecto es False.
    xaxis : str, opcional
        Tipo de eje x ('omega', 'freq', 'norm'). Por defecto es 'omega'.
    fs : float, opcional
        Frecuencia de muestreo. Por defecto es 2*pi.


    Returns
    --------
    fig_id : int
        Identificador de la figura.
    axes_hdl : array_like de Axes
        Manejadores de ejes de la figura.


    Raises
    ------
    ValueError
        Si myFilter no es un array NumPy ni un objeto TransferFunction.
        Si los argumentos fig_id o axes_hdl no son válidos.
        Si xaxis no es uno de los valores permitidos: 'omega', 'freq', 'norm'.


    See Also
    -----------
    :func:`analyze_sys`
    :func:`GroupDelay`
    :func:`pzmap`


    Example
    --------
    >>> fig_id, axes_hdl = bodePlot(myFilter, fig_id=1, axes_hdl='none', filter_description='Filtro pasa bajos', npoints=1000, digital=False, xaxis='omega', fs=2*np.pi)

    
    """

    # Verificar si fig_id es válido
    if not isinstance(fig_id, (str, int)):
        raise ValueError("fig_id debe ser una cadena de caracteres o un número entero.")

    # Verificar si axes_hdl es válido
    if not isinstance(axes_hdl, (str, list, np.ndarray)):
        raise ValueError("axes_hdl debe ser una cadena de caracteres, una lista o un handle de Axes.")

    # Verificar si xaxis es uno de los valores permitidos
    if xaxis not in ['omega', 'freq', 'norm']:
        raise ValueError("xaxis debe ser uno de los siguientes valores: 'omega', 'freq', 'norm'.")

    # Convertir myFilter a un objeto TransferFunction si es un array NumPy
    if isinstance(myFilter, np.ndarray):
        # Convertir sección SOS a una TransferFunction completa
        wholeFilter = sos2tf_analog(myFilter)

        # Obtener todas las singularidades
        this_zzpp = np.abs(np.concatenate([wholeFilter.zeros, wholeFilter.poles]))
        this_zzpp = this_zzpp[this_zzpp > 0]

        # Calcular el eje de frecuencia según las singularidades del filtro completo
        if digital:
            ww = np.linspace(0, np.pi, npoints)
        else:
            ww = np.logspace(np.floor(np.log10(np.min(this_zzpp))) - 1, np.ceil(np.log10(np.max(this_zzpp))) + 1, npoints)

        cant_sos = myFilter.shape[0]
        mag = np.empty((npoints, cant_sos + 1))
        phase = np.empty_like(mag)
        sos_label = []

        # Calcular la respuesta de magnitud y fase para cada sección SOS y el filtro completo
        for ii in range(cant_sos):

            num, den = _one_sos2tf(myFilter[ii, :])
            thisFilter = TransferFunction(num, den)

            this_zzpp = np.abs(np.concatenate([thisFilter.zeros, thisFilter.poles]))
            this_zzpp = this_zzpp[this_zzpp > 0]

            _, mag[:, ii], phase[:, ii] = thisFilter.bode(ww)

            sos_label += [filter_description + ' - SOS {:d}'.format(ii)]

        _, mag[:, cant_sos], phase[:, cant_sos] = wholeFilter.bode(ww)

        sos_label += [filter_description]

        filter_description = sos_label

    else:
        # Si myFilter es un objeto TransferFunction
        cant_sos = 0

        this_zzpp = np.abs(np.concatenate([myFilter.zeros, myFilter.poles]))
        this_zzpp = this_zzpp[this_zzpp > 0]

        if digital:
            ww, mag, phase = myFilter.bode(n=np.linspace(0, np.pi, npoints))
        else:
            ww, mag, phase = myFilter.bode(np.logspace(np.floor(np.log10(np.min(this_zzpp))) - 1, np.ceil(np.log10(np.max(this_zzpp))) + 1, npoints))

    # Convertir frecuencia a Hz si se solicita
    if xaxis == "freq":
        ww = ww / 2 / np.pi
    elif xaxis == "norm":
        if fs is None:
            # Normalizar cada respuesta a su propio Nyquist
            wnorm = 2 * np.pi / myFilter.dt / 2
        else:
            # Normalizado a fs
            wnorm = 2 * np.pi * fs
        ww = ww / wnorm

    # Crear o recuperar figura y ejes
    if fig_id == 'none':
        fig_hdl, axes_hdl = plt.subplots(2, 1, sharex='col')
        fig_id = fig_hdl.number
    else:
        if plt.fignum_exists(fig_id):
            fig_hdl = plt.figure(fig_id)
            axes_hdl = fig_hdl.get_axes()
        else:
            fig_hdl = plt.figure(fig_id)
            axes_hdl = fig_hdl.subplots(2, 1, sharex='col')
            fig_id = fig_hdl.number

    (mag_ax_hdl, phase_ax_hdl) = axes_hdl

    # Graficar respuesta de magnitud
    plt.sca(mag_ax_hdl)
    if digital:
        if filter_description is None:
            aux_hdl = plt.plot(ww, mag)    # Bode magnitude plot
        else:
            aux_hdl = plt.plot(ww, mag, label=filter_description)    # Bode magnitude plot
    else:
        if filter_description is None:
            aux_hdl = plt.semilogx(ww, mag)    # Bode magnitude plot
        else:
            aux_hdl = plt.semilogx(ww, mag, label=filter_description)    # Bode magnitude plot

    if cant_sos > 0:
        # Distinguir respuesta SOS de la respuesta total
        [aa.set_linestyle(':') for aa in aux_hdl[:-1]]
        aux_hdl[-1].set_linewidth(2)

    plt.grid(True)
    plt.ylabel('Magnitud [dB]')
    plt.title('Respuesta de Magnitud')

    if not(filter_description is None):
        mag_ax_hdl.legend()

    # Graficar respuesta de fase
    plt.sca(phase_ax_hdl)
    if digital:
        if filter_description is None:
            aux_hdl = plt.plot(ww, np.pi / 180 * phase)    # Bode phase plot
        else:
            aux_hdl = plt.plot(ww, np.pi / 180 * phase, label=filter_description)    # Bode phase plot
    else:
        if filter_description is None:
            aux_hdl = plt.semilogx(ww, np.pi / 180 * phase)    # Bode phase plot
        else:
            aux_hdl = plt.semilogx(ww, np.pi / 180 * phase, label=filter_description)    # Bode phase plot

    # Escalar los ejes para ajustar
    ylim = plt.gca().get_ylim()
    ticks = np.linspace(start=np.round(ylim[0] / np.pi) * np.pi, stop=np.round(ylim[1] / np.pi) * np.pi, num=5, endpoint=True)
    ylabs = []
    for aa in ticks:
        if aa == 0:
            ylabs += ['0']
        else:
            bb = Fraction(aa / np.pi).limit_denominator(1000000)
            if np.abs(bb.numerator) != 1:
                if np.abs(bb.denominator) != 1:
                    str_aux = r'$\frac{{{:d}}}{{{:d}}} \pi$'.format(bb.numerator, bb.denominator)
                else:
                    str_aux = r'${:d}\pi$'.format(bb.numerator)
            else:
                if np.abs(bb.denominator) == 1:
                    if np.sign(bb.numerator) == -1:
                        str_aux = r'$-\pi$'
                    else:
                        str_aux = r'$\pi$'
                else:
                    if np.sign(bb.numerator) == -1:
                        str_aux = r'$-\frac{{\pi}}{{{:d}}}$'.format(bb.denominator)
                    else:
                        str_aux = r'$\frac{{\pi}}{{{:d}}}$'.format(bb.denominator)
            ylabs += [str_aux]
    plt.yticks(ticks, labels=ylabs)

    if cant_sos > 0:
        # Distinguir respuesta SOS de la respuesta total
        [aa.set_linestyle(':') for aa in aux_hdl[:-1]]
        aux_hdl[-1].set_linewidth(2)

    plt.grid(True)

    if xaxis == "freq":
        plt.xlabel('Frecuencia [Hz]')
    elif xaxis == "norm":
        plt.gca().set_xlim([0, 1])
        if fs is None:
            # Normalizar cada respuesta a su propio Nyquist
            this_fs = 1 / myFilter.dt
        else:
            # Normalizado a fs
            this_fs = fs
        plt.xlabel('Frecuencia normalizada a fs={:3.3f} [#]'.format(this_fs))
    else:
        plt.xlabel('Frecuencia angular [rad/seg]')

    plt.ylabel('Fase [rad]')
    plt.title('Respuesta de Fase')

    if not(filter_description is None):
        phase_ax_hdl.legend()

    return fig_id, axes_hdl

def plot_plantilla(filter_type='', fpass=0.25, ripple=0.5, fstop=0.6, attenuation=40, fs=2):
    """
    Plotea una plantilla de diseño de filtro digital.

    Parameters
    -----------
    filter_type : str, opcional
        Tipo de filtro ('lowpass', 'highpass', 'bandpass', 'bandstop'). Por defecto es 'lowpass'.
    fpass : float o tupla, opcional
        Frecuencia de paso o tupla de frecuencias de paso para los filtros 'bandpass' o 'bandstop'.
    ripple : float, opcional
        Máxima ondulación en la banda de paso (en dB). Por defecto es 0.5 dB.
    fstop : float o tupla, opcional
        Frecuencia de detención o tupla de frecuencias de detención para los filtros 'bandpass' o 'bandstop'.
    attenuation : float, opcional
        Atenuación mínima en la banda de detención (en dB). Por defecto es 40 dB.
    fs : float, opcional
        Frecuencia de muestreo. Por defecto es 2.
        
    Returns
    --------
    None


    Raises
    ------
    ValueError
        Si los argumentos no son del tipo o valor correcto.
    

    See Also
    -----------
    :func:`analyze_sys`

    
    Example
    --------
    >>> plot_plantilla(filter_type='bandpass', fpass=(0.2, 0.4), ripple=0.3, fstop=(0.1, 0.5), attenuation=50, fs=2)

    """

    if not isinstance(fpass, (tuple, np.ndarray, Integral, Real)):
        raise ValueError("fpass debe ser un float o una tupla de frecuencias de paso.")
    
    if not isinstance(fstop, (tuple, np.ndarray, Integral, Real)):
        raise ValueError("fstop debe ser un float o una tupla de frecuencias de detención.")
    
    if not isinstance(attenuation, (tuple, np.ndarray, Integral, Real)):
        try:
            attenuation = np.float64(attenuation)
        except ValueError:
            raise ValueError("attenuation debe ser un valor numérico o convertible a float.")    
    
    if not isinstance(ripple, (tuple, np.ndarray, Integral, Real)):
        try:
            ripple = np.float64(ripple)
        except ValueError:
            raise ValueError("attenuation debe ser un valor numérico o convertible a float.")    
    
    if not isinstance(fs, (Integral, Real)):
        try:
            fs = np.float64(fs)
        except ValueError:
            raise ValueError("fs debe ser un valor numérico o convertible a float.")    

    # Obtener los límites actuales de los ejes
    xmin, xmax, ymin, ymax = plt.axis()

    # Banda de paso digital
    plt.fill([xmin, xmin, fs / 2, fs / 2], [ymin, ymax, ymax, ymin], 'lightgreen', alpha=0.2, lw=1, label='Banda de paso digital')

    # analizar los valores de la plantilla para ver qué tipo de plantilla es
    tipos_permitidos = ['lowpass', 'highpass', 'bandpass', 'bandstop']
    if filter_type not in tipos_permitidos:
            
        if isinstance(fpass, (tuple, np.ndarray)) and isinstance(fstop, (tuple, np.ndarray)):
            if fstop[0] < fpass[0]:
                filter_type = 'bandpass'
            else:                
                filter_type = 'bandstop'
        else:
            if fstop < fpass:
                filter_type = 'highpass'
            else:                
                filter_type = 'lowpass'

    if filter_type == 'lowpass':
        # Definir regiones de banda de detención y banda de paso para el filtro pasa bajos
        fstop_start = fstop
        fstop_end = xmax
        fpass_start = xmin
        fpass_end = fpass
    elif filter_type == 'highpass':
        # Definir regiones de banda de detención y banda de paso para el filtro pasa altos
        fstop_start = xmin
        fstop_end = fstop
        fpass_start = fpass
        fpass_end = xmax
    elif filter_type == 'bandpass':
        if len(fpass) == 2 and len(fstop) == 2:
            fstop_start = xmin
            fstop_end = fstop[0]
            fpass_start = fpass[0]
            fpass_end = fpass[1]
            fstop2_start = fstop[1]
            fstop2_end = xmax
        else:
            raise ValueError("En modo bandpass, fpass y fstop deben ser tuplas con 2 valores.")
    elif filter_type == 'bandstop':
        if len(fpass) == 2 and len(fstop) == 2:
            fpass_start = xmin
            fpass_end = fpass[0]
            fstop_start = fstop[0]
            fstop_end = fstop[1]
            fpass2_start = fpass[1]
            fpass2_end = xmax
        else:
            raise ValueError("En modo bandstop, fpass y fstop deben ser tuplas con 2 valores.")
    else:
        raise ValueError("filtro_type debe ser 'lowpass', 'highpass', 'bandpass', o 'bandstop'.")

    # Plotea regiones de banda de detención y banda de paso
    plt.fill([fstop_start, fstop_end, fstop_end, fstop_start], [-attenuation, -attenuation, ymax, ymax], 'lightgrey', alpha=0.4, hatch='x', lw=1, ls='--', ec='k')
    plt.fill([fpass_start, fpass_start, fpass_end, fpass_end], [ymin, -ripple, -ripple, ymin], 'lightgrey', alpha=0.4, hatch='x', lw=1, ls='--', ec='k', label='Plantilla')
    
    # Plotea región adicional de banda de detención para filtro pasa banda
    if filter_type == 'bandpass':
        plt.fill([fstop2_start, fstop2_end, fstop2_end, fstop2_start], [-attenuation, -attenuation, ymax, ymax], 'lightgrey', alpha=0.4, hatch='x', lw=1, ls='--', ec='k')
    
    # Plotea región adicional de banda de paso para filtro rechaza banda
    if filter_type == 'bandstop':
        plt.fill([fpass2_start, fpass2_start, fpass2_end, fpass2_end], [ymin, -ripple, -ripple, ymin], 'lightgrey', alpha=0.4, hatch='x', lw=1, ls='--', ec='k')
    
    # Establece los límites de los ejes
    plt.axis([xmin, xmax, np.max([ymin, -100]), np.max([ymax, 1])])
    
   ########################
  ## Funciones internas #
 ########################
#%%

def _nearest_real_complex_idx(fro, to, which):
    '''
    Obtiene el índice del siguiente elemento real o complejo más cercano basado en la distancia.

    Parameters
    -----------
    fro : array_like
        Arreglo de partida que contiene los elementos a comparar.
    to : array_like
        Valor de referencia para encontrar el elemento más cercano.
    which : {'real', 'complex'}
        Especifica si se busca el elemento real o complejo más cercano.


    Returns
    --------
    int
        Índice del elemento más cercano en el arreglo de partida.


    Raises
    ------
    AssertionError
        Si el argumento 'which' no es 'real' o 'complex'.


    See Also
    -----------
    :func:`zpk2sos_analog`

    
    Example
    --------
    >>> fro = np.array([1, 2, 3, 4])
    >>> to = 2.5
    >>> nearest_idx = _nearest_real_complex_idx(fro, to, 'real')
    >>> print("El índice del elemento real más cercano a", to, "es:", nearest_idx)


    Esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.

    '''
    # Verificar que 'which' sea 'real' o 'complex'
    assert which in ('real', 'complex')
    # Ordenar los índices según la distancia al valor de referencia
    order = np.argsort(np.abs(fro - to))
    # Crear una máscara para seleccionar elementos reales o complejos
    mask = np.isreal(fro[order])
    if which == 'complex':
        mask = ~mask
    # Devolver el índice del primer elemento que cumple con la condición
    return order[np.nonzero(mask)[0][0]]

def _cplxreal(z, tol=None):
    """
    Separa en partes complejas y reales, combinando pares conjugados.

    El vector de entrada unidimensional `z` se divide en sus elementos complejos (`zc`) y reales (`zr`).
    Cada elemento complejo debe ser parte de un par conjugado complejo, que se combinan en un solo número
    (con parte imaginaria positiva) en la salida. Dos números complejos se consideran un par conjugado si sus partes
    real e imaginaria difieren en magnitud por menos de ``tol * abs(z)``.

    Parameters
    -----------
    z : array_like
        Vector de números complejos para ordenar y dividir.
    tol : float, opcional
        Tolerancia relativa para probar la realidad e igualdad conjugada.
        El valor predeterminado es ``100 * espaciado(1)`` del tipo de datos de `z`
        (es decir, 2e-14 para float64).

    Returns
    --------
    zc : ndarray
        Elementos complejos de `z`, donde cada par se representa por un solo valor
        con parte imaginaria positiva, ordenada primero por parte real y luego
        por magnitud de la parte imaginaria. Los pares se promedian cuando se combinan
        para reducir el error.
    zr : ndarray
        Elementos reales de `z` (aquellos que tienen parte imaginaria menor que
        `tol` veces su magnitud), ordenados por valor.


    Raises
    ------
    ValueError
        Si hay números complejos en `z` para los cuales no se puede encontrar un conjugado.
    

    See Also
    ---------
    :func:`zpk2sos_analog`


    Exampless
    ---------
    >>> a = [4, 3, 1, 2-2j, 2+2j, 2-1j, 2+1j, 2-1j, 2+1j, 1+1j, 1-1j]
    >>> zc, zr = _cplxreal(a)
    >>> print(zc)
    [ 1.+1.j  2.+1.j  2.+1.j  2.+2.j]
    >>> print(zr)
    [ 1.  3.  4.]

    Esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.
    """

    z = np.atleast_1d(z)
    if z.size == 0:
        return z, z
    elif z.ndim != 1:
        raise ValueError('_cplxreal solo acepta entradas 1-D')

    if tol is None:
        # Obtener la tolerancia del dtype de la entrada
        tol = 100 * np.finfo((1.0 * z).dtype).eps

    # Ordenar por parte real, magnitud de la parte imaginaria (acelerar más la clasificación)
    z = z[np.lexsort((abs(z.imag), z.real))]

    # Separar reales de pares conjugados
    real_indices = abs(z.imag) <= tol * abs(z)
    zr = z[real_indices].real

    if len(zr) == len(z):
        # La entrada es completamente real
        return np.array([]), zr

    # Separar mitades positivas y negativas de conjugados
    z = z[~real_indices]
    zp = z[z.imag > 0]
    zn = z[z.imag < 0]

    if len(zp) != len(zn):
        raise ValueError('El array contiene un valor complejo sin su conjugado correspondiente.')

    # Encontrar carreras de partes reales (aproximadamente) iguales
    same_real = np.diff(zp.real) <= tol * abs(zp[:-1])
    diffs = np.diff(np.concatenate(([0], same_real, [0])))
    run_starts = np.nonzero(diffs > 0)[0]
    run_stops = np.nonzero(diffs < 0)[0]

    # Ordenar cada carrera por sus partes imaginarias
    for i in range(len(run_starts)):
        start = run_starts[i]
        stop = run_stops[i] + 1
        for chunk in (zp[start:stop], zn[start:stop]):
            chunk[...] = chunk[np.lexsort([abs(chunk.imag)])]

    # Verificar que los negativos coincidan con los positivos
    if any(abs(zp - zn.conj()) > tol * abs(zn)):
        raise ValueError('El array contiene un valor complejo sin su conjugado correspondiente.')

    # Promediar la inexactitud numérica en partes reales vs imaginarias de pares
    zc = (zp + zn.conj()) / 2

    return zc, zr

def _one_sos2tf(mySOS):
    """
    Convierte una sección de segundo orden (SOS) en coeficientes de función de transferencia.

    Parameters
    -----------
    mySOS : array_like
        Vector que define una sección de segundo orden (SOS) del sistema.


    Returns
    --------
    num : ndarray
        Coeficientes del numerador de la función de transferencia.
    den : ndarray
        Coeficientes del denominador de la función de transferencia.


    Raises
    ------
    ValueError
        Si la entrada no es un vector con al menos 6 elementos.


    See Also
    -----------
    :func:`sos2tf_analog`


    Examples
    --------
    >>> mySOS = [1, -1.9, 1, 1, -1.6, 0.64]
    >>> num, den = _one_sos2tf(mySOS)
    >>> print(num)
    [1, -1.9, 1]
    >>> print(den)
    [1, -1.6, 0.64]

    Esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.
    """
    if not isinstance(mySOS, (list, tuple, np.ndarray)):
        raise ValueError("El argumento 'mySOS' debe ser un vector.")
    if len(mySOS) < 6:
        raise ValueError("El argumento 'mySOS' debe tener al menos 6 elementos.")

    # Verificar ceros en los coeficientes de orden superior
    if mySOS[0] == 0 and mySOS[1] == 0:
        num = mySOS[2]
    elif mySOS[0] == 0:
        num = mySOS[1:3]
    else:
        num = mySOS[:3]

    if mySOS[3] == 0 and mySOS[4] == 0:
        den = mySOS[-1]
    elif mySOS[3] == 0:
        den = mySOS[4:]
    else:
        den = mySOS[3:]

    return num, den

def _build_poly_str(this_poly):
    """
    Construye una cadena de caracteres que representa un polinomio.

    Parameters
    -----------
    this_poly : ndarray
        Coeficientes del polinomio.


    Returns
    --------
    str
        Cadena de caracteres que representa el polinomio.


    Raises
    ------
    ValueError
        Si `this_poly` no es un array de numpy.


    See Also
    -----------
    :func:`pretty_print_lti`
    :func:`pretty_print_bicuad_omegayq`


    Examples
    --------
    >>> this_poly = np.array([1, -2, 3])
    >>> poly_str = _build_poly_str(this_poly)
    >>> print(poly_str)
    's^2 - 2 s + 3'


    Esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.
    """
    if not isinstance(this_poly, np.ndarray):
        raise ValueError("El argumento 'this_poly' debe ser un array de numpy.")

    poly_str = ''

    for ii in range(this_poly.shape[0]):

        if this_poly[ii] != 0.0:

            if (this_poly.shape[0] - 2) == ii:
                poly_str += '+ s '

            elif (this_poly.shape[0] - 1) != ii:
                poly_str += '+ s^{:d} '.format(this_poly.shape[0] - ii - 1)

            if (this_poly.shape[0] - 1) == ii:
                poly_str += '+ {:3.4g} '.format(this_poly[ii])
            else:
                if this_poly[ii] != 1.0:
                    poly_str += '\,\, {:3.4g} '.format(this_poly[ii])

    return poly_str[2:]

def _build_omegayq_str(this_quad_poly, den=np.array([])):
    """
    Construye una cadena de caracteres que representa omega y Q.


    Parameters
    -----------
    this_quad_poly : ndarray
        Coeficientes del polinomio cuadrático.
    den : ndarray, opcional
        Coeficientes del denominador. El valor predeterminado es np.array([]).


    Raises
    ------
    ValueError
        Si `this_quad_poly` no es un array de numpy.
        Si `den` no es un array de numpy.

        
    Returns
    --------
    str
        Cadena de caracteres que representa omega y Q.


    See Also
    -----------
    :func:`pretty_print_bicuad_omegayq`


    Examples
    --------
    >>> this_quad_poly = np.array([1, 2, 3])
    >>> den = np.array([4, 5, 6])
    >>> omegaq_str = _build_omegayq_str(this_quad_poly, den)
    >>> print(omegaq_str)
    's\,0.08333\,\frac{2}{0.1667}'


    Esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.
    """
    if not isinstance(this_quad_poly, np.ndarray):
        raise ValueError("El argumento 'this_quad_poly' debe ser un array de numpy.")

    if not isinstance(den, np.ndarray):
        raise ValueError("El argumento 'den' debe ser un array de numpy.")

    if den.shape[0] > 0:
        # Estilo de numerador para banda pasante s. hh . omega/ Q
        
        omega = np.sqrt(den[2]) # del denominador
        Q = omega / den[1] # del denominador
        
        hh = this_quad_poly[1] * Q / omega
        
        poly_str = r's\,{:3.4g}\,\frac{{{:3.4g}}}{{{:3.4g}}}'.format(hh, omega, Q )
    
    else:
        # Todos los demás polinomios cuadráticos completos
        omega = np.sqrt(this_quad_poly[2])
        
        if this_quad_poly[1] == 0:
            poly_str = r's^2 + {:3.4g}^2'.format(omega)
        else:
            Q = omega / this_quad_poly[1]
            poly_str = r's^2 + s \frac{{{:3.4g}}}{{{:3.4g}}} + {:3.4g}^2'.format(omega, Q, omega)
        
                
    return poly_str

def _complementaryColor(my_hex):
    """
    Returns el color RGB complementario.

    Parameters
    -----------
    my_hex : str
        Código hexadecimal del color.


    Raises
    ------
    ValueError
        Si `my_hex` no es una cadena de caracteres válida o no tiene la longitud correcta.


    Returns
    --------
    str
        Código hexadecimal del color complementario.


    See Also
    -----------
    :func:`pzmap`


    Examples
    --------
    >>> _complementaryColor('FFFFFF')
    '000000'


    Esta documentación ha sido generada por ChatGPT, una IA desarrollada por OpenAI.
    """
    
    if not isinstance(my_hex, str) or len(my_hex) != 7: # '#' + 6 RGB
        raise ValueError("El argumento 'my_hex' debe ser una cadena de 6 caracteres representando un código hexadecimal válido.")
    
    if my_hex[0] == '#':
        my_hex = my_hex[1:]
    
    rgb = (my_hex[0:2], my_hex[2:4], my_hex[4:6])
    comp = ['%02X' % (255 - int(a, 16)) for a in rgb]
    return '#' + ''.join(comp)
