#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:14:13 2023

@author: mariano
"""

import numpy as np

import sympy as sp

from .general import print_latex


#%%
   #############################################################################
  ## (Simbólicas) Funciones de conversión de matrices de cuadripolos lineales #
 #############################################################################
#%%

def S2Ts_s(Spar):
    """
    Convierte una matriz de parámetros de dispersión (S) simbólica 
    en el modelo de parámetros de transferencia de dispersión (Ts).

    Esta función toma una matriz simbólica que representa los parámetros de dispersión (S) de un sistema y calcula la matriz de parámetros de transferencia de dispersión (Ts) correspondiente.

    
    Parameters
    ----------
    Spar : Symbolic Matrix
        Matriz de parámetros de dispersión S.


    Returns
    -------
    Ts : Symbolic Matrix
        Matriz de parámetros de transferencia de dispersión Ts.


    Raises
    ------
    ValueError
        Si Spar no es una instancia de Symbolic Matrix.
        Si Spar no tiene el formato correcto [ [Spar_11, Spar_12], [Spar_21, Spar_22] ].
        Si Spar_12 es nulo.


    See Also
    --------
    :func:`Ts2S_s`
    :func:`S2Tabcd_s`
    :func:`Model_conversion`


    Examples
    --------
    >>> Spar = sp.Matrix([[sp.symbols('Z11'), sp.symbols('Z12')],
    ...                   [sp.symbols('Z21'), sp.symbols('Z22')]])
    >>> Ts = S2Ts_s(Spar)
    >>> print(Ts)
    Matrix([[1, -Z12], [Z11, Z22/Z21]])


    Notes
    -----
    - La matriz Spar debe tener la forma [ [Spar_11, Spar_12], [Spar_21, Spar_22] ].
    - Spar_12 no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    """
    # Verificar si Spar es una instancia de Symbolic Matrix
    if not isinstance(Spar, sp.Matrix):
        raise ValueError("Spar debe ser una instancia de Symbolic Matrix")
    
    # Verificar que Spar tenga el formato correcto
    if Spar.shape != (2, 2):
        raise ValueError("Spar debe tener el formato [ [Spar_11, Spar_12], [Spar_21, Spar_22] ]")
    
    # Verificar que Spar_12 no sea nulo
    if Spar[0, 1] == 0:
        raise ValueError("Spar_12 no puede ser nulo")
    
    # Inicialización de la matriz de parámetros de transferencia
    Ts = sp.Matrix([[0, 0], [0, 0]])
    
    # Cálculo de los elementos de la matriz Ts
    Ts[0, 0] = sp.Rational('1')                               # Ts11 = 1
    Ts[0, 1] = -Spar[1, 1]                                    # Ts12 = -S22
    Ts[1, 0] = Spar[0, 0]                                     # Ts21 = S11
    Ts[1, 1] = -sp.simplify(sp.expand(sp.Determinant(Spar)))  # Ts22 = -det(Spar)

    return sp.simplify(sp.expand(1 / Spar[1, 0] * Ts))  # Tsxx = 1/S21 * Tsxx 

def Ts2S_s(Ts):
    """
    Convierte una matriz de transferencia de scattering (Ts) simbólica 
    al modelo de parámetros scattering (S).

    Parameters
    ----------
    Ts : Symbolic Matrix
        Matriz de parámetros S.

    Returns
    -------
    Spar : Symbolic Matrix
        Matriz de parámetros de scattering.

    Raises
    ------
    ValueError
        Si Ts no es una instancia de Symbolic Matrix.
        Si Ts no tiene el formato correcto [ [Ts_11, Ts_12], [Ts_21, Ts_22] ].
        Si Ts_11 es nulo.

    See Also
    --------
    :func:`S2Ts_s`
    :func:`S2Tabcd_s`
    :func:`Model_conversion`

    Examples
    --------
    >>> Ts = sp.Matrix([[sp.symbols('A'), sp.symbols('B')],
    ...                 [sp.symbols('C'), sp.symbols('D')]])
    >>> Spar = Ts2S_s(Ts)
    >>> print(Spar)
    Matrix([[C, A], [B, -D]])

    Notes
    -----
    - La matriz Ts debe tener la forma [ [Ts_11, Ts_12], [Ts_21, Ts_22] ].
    - Ts_11 no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    """
    # Verificar si Ts es una instancia de Symbolic Matrix
    if not isinstance(Ts, sp.Matrix):
        raise ValueError("Ts debe ser una instancia de Symbolic Matrix")
    
    # Verificar que Ts tenga el formato correcto
    if Ts.shape != (2, 2):
        raise ValueError("Ts debe tener el formato [ [Ts_11, Ts_12], [Ts_21, Ts_22] ]")
    
    # Verificar que Ts_11 no sea nulo
    if Ts[0, 0] == 0:
        raise ValueError("Ts_11 no puede ser nulo")
    
    # Inicialización de la matriz de parámetros de scattering
    Spar = sp.Matrix([[0, 0], [0, 0]])
    
    # Cálculo de los elementos de la matriz Spar
    Spar[0, 0] = Ts[1, 0]                                    # S11 = TS21
    Spar[0, 1] = sp.simplify(sp.expand(sp.Determinant(Ts)))  # S12 = -det(Ts)
    Spar[1, 0] = sp.Rational('1')                            # S21 = 1
    Spar[1, 1] = -Ts[0, 1]                                   # S22 = -Ts12

    return sp.simplify(sp.expand(1 / Ts[0, 0] * Spar)) # Sxx = 1/Ts11 * Sxx 

def Ts2Tabcd_s(Ts, Z0=sp.Rational('1')):
    """Converts a symbolic scattering parameter matrix (Ts) to the symbolic ABCD or Tabcd model.

    This function converts a symbolic scattering parameter matrix (Ts) to the symbolic ABCD or Tabcd model.

    Parameters
    ----------
    Ts : Symbolic Matrix
        The Ts parameter matrix.
    Z0 : sp.Expr, optional
        The reference impedance, defaults to 1.

    Returns
    -------
    Tabcd : Symbolic Matrix
        The ABCD parameter matrix.

    Raises
    ------
    ValueError
        If Ts is not an instance of sp.Matrix.
        If Z0 is not an instance of sp.Expr.


    See Also
    --------
    :func:`Ts2S_s`
    :func:`S2Tabcd_s`
    :func:`Model_conversion`


    Examples
    --------
    >>> Ts = sp.MatrixSymbol('Ts', 2, 2)
    >>> Ts2Tabcd_s(Ts)
    Tabcd

    """
    # Check if Ts is an instance of sp.Matrix
    if not isinstance(Ts, sp.Matrix):
        raise ValueError("Ts must be an instance of sp.Matrix.")
    
    # Check if Z0 is an instance of sp.Expr
    if not isinstance(Z0, sp.Expr):
        raise ValueError("Z0 must be an instance of sp.Expr.")

    # Convert Ts to S-parameter matrix and then to ABCD matrix
    return sp.simplify(sp.expand(S2Tabcd_s(Ts2S_s(Ts), Z0)))

def Tabcd2S_s(Tabcd, Z0=sp.Rational('1')):
    '''
    Convierte una matriz de parámetros ABCD (Tabcd) simbólica 
    al modelo de parámetros scattering (S).

    Parameters
    ----------
    Tabcd : Symbolic Matrix
        Matriz de parámetros ABCD.
    Z0 : sympy expression, optional
        Impedancia característica del medio. Por defecto es 1.

    Returns
    -------
    Spar : Symbolic Matrix
        Matriz de parámetros de scattering.

    Raises
    ------
    ValueError
        Si Tabcd no es una instancia de Symbolic Matrix.
        Si Tabcd no tiene el formato correcto [ [A, B], [C, D] ].
        Si la matriz Tabcd no es invertible.


    See Also
    --------
    :func:`Ts2S_s`
    :func:`S2Tabcd_s`
    :func:`Model_conversion`


    Examples
    --------
    >>> Tabcd = sp.Matrix([[sp.symbols('A'), sp.symbols('B')],
    ...                    [sp.symbols('C'), sp.symbols('D')]])
    >>> Spar = Tabcd2S_s(Tabcd)
    >>> print(Spar)
    Matrix([[C/A, 1/B], [A/D - B*C/(A*D), -1/D]])

    Notes
    -----
    - La matriz Tabcd debe tener el formato [ [A, B], [C, D] ].
    - La matriz Tabcd debe ser invertible para realizar la conversión correctamente.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si Tabcd es una instancia de Symbolic Matrix
    if not isinstance(Tabcd, sp.Matrix):
        raise ValueError("Tabcd debe ser una instancia de Symbolic Matrix")

    # Verificar que Tabcd tenga el formato correcto
    if Tabcd.shape != (2, 2):
        raise ValueError("Tabcd debe tener el formato [ [A, B], [C, D] ]")

    # Verificar si la matriz Tabcd es invertible
    if Tabcd.det() == 0:
        raise ValueError("La matriz Tabcd no es invertible")

    # Inicialización de la matriz de parámetros de scattering
    Spar = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Spar
    common = Tabcd[0, 0] + Tabcd[0, 1]/Z0 + Tabcd[1, 0]*Z0 + Tabcd[1, 1]
    Spar[0, 0] = Tabcd[0, 0] + Tabcd[0, 1]/Z0 - Tabcd[1, 0]*Z0 - Tabcd[1, 1]
    Spar[0, 1] = 2 * sp.simplify(sp.expand(sp.Determinant(Tabcd)))
    Spar[1, 0] = sp.Rational('2')
    Spar[1, 1] = -Tabcd[0, 0] + Tabcd[0, 1]/Z0 - Tabcd[1, 0]*Z0 + Tabcd[1, 1]

    return sp.simplify(sp.expand(1 / common * Spar))

def S2Tabcd_s(Spar, Z0=sp.Rational('1')):
    '''
    Convierte una matriz de parámetros scattering (S) simbólica 
    al modelo de parámetros ABCD (Tabcd).

    Parameters
    ----------
    Spar : Symbolic Matrix
        Matriz de parámetros S.
    Z0 : sympy expression, optional
        Impedancia característica del medio. Por defecto es 1.

    Returns
    -------
    Tabcd : Symbolic Matrix
        Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si Spar no es una instancia de Symbolic Matrix.
        Si Spar no tiene el formato correcto [ [S11, S12], [S21, S22] ].
        Si Spar[1, 0] es nulo.


    See Also
    --------
    :func:`Tabcd2S_s`
    :func:`Y2Tabcd_s`
    :func:`Model_conversion`


    Examples
    --------
    >>> Spar = sp.Matrix([[sp.symbols('S11'), sp.symbols('S12')],
    ...                   [sp.symbols('S21'), sp.symbols('S22')]])
    >>> Tabcd = S2Tabcd_s(Spar)
    >>> print(Tabcd)
    Matrix([[(1 - S11)*(1 + S22) + S12*S21, Z0*((1 + S11)*(1 + S22) - S12*S21)], [1/Z0*((1 - S11)*(1 - S22) - S12*S21), (1 - S11)*(1 + S22) + S12*S21]])

    Notes
    -----
    - La matriz Spar debe tener el formato [ [S11, S12], [S21, S22] ].
    - Spar[1, 0] no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si Spar es una instancia de Symbolic Matrix
    if not isinstance(Spar, sp.Matrix):
        raise ValueError("Spar debe ser una instancia de Symbolic Matrix")

    # Verificar que Spar tenga el formato correcto
    if Spar.shape != (2, 2):
        raise ValueError("Spar debe tener el formato [ [S11, S12], [S21, S22] ]")

    # Verificar que Spar[1, 0] no sea nulo
    if Spar[1, 0] == 0:
        raise ValueError("Spar[1, 0] no puede ser nulo")

    # Inicialización de la matriz de parámetros ABCD
    Tabcd = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Tabcd
    common = 2 * Spar[1, 0]
    Tabcd[0, 0] = (1 - Spar[0, 0]) * (1 + Spar[1, 1]) + Spar[1, 0] * Spar[0, 1]
    Tabcd[0, 1] = Z0 * ((1 + Spar[0, 0]) * (1 + Spar[1, 1]) - Spar[1, 0] * Spar[0, 1])
    Tabcd[1, 0] = 1 / Z0 * ((1 - Spar[0, 0]) * (1 - Spar[1, 1]) - Spar[1, 0] * Spar[0, 1])
    Tabcd[1, 1] = (1 - Spar[0, 0]) * (1 + Spar[1, 1]) + Spar[1, 0] * Spar[0, 1]

    return sp.simplify(sp.expand(1 / common * Tabcd))

def Y2Tabcd_s(YY):
    """
    Convierte una matriz de admitancia de dos puertos (YY) simbólica 
    al modelo de parámetros ABCD (Tabcd).

    Parameters
    ----------
    YY : Symbolic Matrix
        Matriz de admitancia de dos puertos.

    Returns
    -------
    TT : Symbolic Matrix
        Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si YY no es una instancia de Symbolic Matrix.
        Si YY no tiene el formato correcto [ [Y11, Y12], [Y21, Y22] ].
        Si Y21 es nulo.


    See Also
    --------
    :func:`Ts2S_s`
    :func:`Tabcd2Y_s`
    :func:`Model_conversion`


    Examples
    --------
    >>> YY = sp.Matrix([[sp.symbols('Y11'), sp.symbols('Y12')],
    ...                 [sp.symbols('Y21'), sp.symbols('Y22')]])
    >>> TT = Y2Tabcd_s(YY)
    >>> print(TT)
    Matrix([[-Y22/Y21, -1/Y21], [-D/Y21, -Y22/Y21]])

    Notes
    -----
    - La matriz YY debe tener el formato [ [Y11, Y12], [Y21, Y22] ].
    - YY[1, 0] no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    """
    # Verificar si YY es una instancia de Symbolic Matrix
    if not isinstance(YY, sp.Matrix):
        raise ValueError("YY debe ser una instancia de Symbolic Matrix")

    # Verificar que YY tenga el formato correcto
    if YY.shape != (2, 2):
        raise ValueError("YY debe tener el formato [ [Y11, Y12], [Y21, Y22] ]")

    # Verificar que YY[1, 0] no sea nulo
    if YY[1, 0] == 0:
        raise ValueError("Y21 no puede ser nulo")

    # Inicialización de la matriz de parámetros ABCD
    TT = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz TT
    TT[0, 0] = -YY[1, 1] 
    TT[0, 1] = -sp.Rational('1')
    TT[1, 0] = -sp.expand(sp.Determinant(YY))
    TT[1, 1] = -YY[1, 1]

    return sp.simplify(sp.expand(1/YY[1, 0] * TT))

def Z2Tabcd_s(ZZ):
    '''
    Convierte la matriz de impedancia (ZZ) simbólica 
    al modelo de parámetros ABCD (Tabcd).

    Parameters
    ----------
    ZZ : Symbolic Matrix
        Matriz de impedancia.

    Returns
    -------
    TT : Symbolic Matrix
        Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si ZZ no es una instancia de Symbolic Matrix.
        Si ZZ no tiene el formato correcto [ [Z11, Z12], [Z21, Z22] ].
        Si Z21 es nulo.


    See Also
    --------
    :func:`Tabcd2Z_s`
    :func:`Tabcd2Y_s`
    :func:`Model_conversion`


    Examples
    --------
    >>> ZZ = sp.Matrix([[sp.symbols('Z11'), sp.symbols('Z12')],
    ...                 [sp.symbols('Z21'), sp.symbols('Z22')]])
    >>> TT = Z2Tabcd_s(ZZ)
    >>> print(TT)
    Matrix([[Z11/Z21, Z12*Z22/Z21 - Z11], [1/Z21, Z22/Z21]])

    Notes
    -----
    - La matriz ZZ debe tener el formato [ [Z11, Z12], [Z21, Z22] ].
    - Z21 no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si ZZ es una instancia de Symbolic Matrix
    if not isinstance(ZZ, sp.Matrix):
        raise ValueError("ZZ debe ser una instancia de Symbolic Matrix")

    # Verificar que ZZ tenga el formato correcto
    if ZZ.shape != (2, 2):
        raise ValueError("ZZ debe tener el formato [ [Z11, Z12], [Z21, Z22] ]")

    # Verificar que ZZ[1, 0] no sea nulo
    if ZZ[1, 0] == 0:
        raise ValueError("Z21 no puede ser nulo")

    # Inicialización de la matriz de parámetros ABCD
    TT = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz TT
    TT[0, 0] = ZZ[0, 0] / ZZ[1, 0]
    TT[0, 1] = sp.expand(sp.Determinant(ZZ)) / ZZ[1, 0]
    TT[1, 0] = 1 / ZZ[1, 0]
    TT[1, 1] = ZZ[1, 1] / ZZ[1, 0]

    return sp.simplify(sp.expand(TT))

def Tabcd2Z_s(TT):
    '''
    Convierte una matriz de parámetros ABCD (TT) simbólica 
    al modelo de impedancia de dos puertos (ZZ).

    Parameters
    ----------
    TT : Symbolic Matrix
        Matriz de parámetros ABCD.

    Returns
    -------
    ZZ : Symbolic Matrix
        Matriz de impedancia de dos puertos.

    Raises
    ------
    ValueError
        Si TT no es una instancia de Symbolic Matrix.
        Si TT no tiene el formato correcto [ [A, B], [C, D] ].
        Si C es nulo.


    See Also
    --------
    :func:`Z2Tabcd_s`
    :func:`Tabcd2Y_s`
    :func:`Model_conversion`


    Examples
    --------
    >>> TT = sp.Matrix([[sp.symbols('T11'), sp.symbols('T12')],
    ...                 [sp.symbols('T21'), sp.symbols('T22')]])
    >>> ZZ = Tabcd2Z_s(TT)
    >>> print(ZZ)
    Matrix([[Z11/Z21, DT/Z21], [1/Z21, Z22/Z21]])

    Notes
    -----
    - La matriz TT debe tener el formato [ [A, B], [C, D] ].
    - C no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si TT es una instancia de Symbolic Matrix
    if not isinstance(TT, sp.Matrix):
        raise ValueError("TT debe ser una instancia de Symbolic Matrix")

    # Verificar que TT tenga el formato correcto
    if TT.shape != (2, 2):
        raise ValueError("TT debe tener el formato [ [A, B], [C, D] ]")

    # Verificar que TT[1, 0] no sea nulo
    if TT[1, 0] == 0:
        raise ValueError("C no puede ser nulo")

    # Inicialización de la matriz de impedancia de dos puertos
    ZZ = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz ZZ
    ZZ[0, 0] = sp.simplify(sp.expand(TT[0, 0] / TT[1, 0]))
    ZZ[0, 1] = sp.simplify(sp.expand(sp.Determinant(TT) / TT[1, 0]))
    ZZ[1, 0] = sp.simplify(sp.expand(1 / TT[1, 0]))
    ZZ[1, 1] = sp.simplify(sp.expand(TT[1, 1] / TT[1, 0]))

    return ZZ

def Tabcd2Y_s(TT):
    '''
    Convierte una matriz de parámetros ABCD (TT) simbólica 
    al modelo de admitancia de dos puertos (YY).

    Parameters
    ----------
    TT : Symbolic Matrix
        Matriz de parámetros ABCD.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz de admitancia de dos puertos.

    Raises
    ------
    ValueError
        Si TT no es una instancia de Symbolic Matrix.
        Si TT no tiene el formato correcto [ [A, B], [C, D] ].
        Si B es nulo.


    See Also
    --------
    :func:`Y2Tabcd_s`
    :func:`Tabcd2Z_s`
    :func:`Model_conversion`


    Examples
    --------
    >>> TT = sp.Matrix([[sp.symbols('T11'), sp.symbols('T12')],
    ...                 [sp.symbols('T21'), sp.symbols('T22')]])
    >>> YY = Tabcd2Y_s(TT)
    >>> print(YY)
    Matrix([[D/B, -DT/B], [-1/B, A/B]])

    Notes
    -----
    - La matriz TT debe tener el formato [ [A, B], [C, D] ].
    - B no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si TT es una instancia de Symbolic Matrix
    if not isinstance(TT, sp.Matrix):
        raise ValueError("TT debe ser una instancia de Symbolic Matrix")

    # Verificar que TT tenga el formato correcto
    if TT.shape != (2, 2):
        raise ValueError("TT debe tener el formato [ [A, B], [C, D] ]")

    # Verificar que TT[0, 1] no sea nulo
    if TT[0, 1] == 0:
        raise ValueError("B no puede ser nulo")

    # Inicialización de la matriz de admitancia de dos puertos
    YY = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz YY
    YY[0, 0] = sp.simplify(sp.expand(TT[1, 1] / TT[0, 1]))
    YY[0, 1] = sp.simplify(sp.expand(-sp.Determinant(TT) / TT[0, 1]))
    YY[1, 0] = sp.simplify(sp.expand(-1 / TT[0, 1]))
    YY[1, 1] = sp.simplify(sp.expand(TT[0, 0] / TT[0, 1]))

    return YY

def I2Tabcd_s(gamma, z01, z02=None):
    '''
    Convierte una ganancia compleja expresada en neppers (gamma) 
    y la impedancia de referencia (z01,2) en una matriz de parámetros ABCD (TT).

    Parameters
    ----------
    gamma : Symbol
        Ganancia compleja expresada en neppers (Re{gamma}) y radianes (Im{gamma}).
        
    z01 : Symbol
        Impedancia de referencia del puerto 1.

    z02 : Symbol, opcional
        Impedancia de referencia del puerto 2. Si no se proporciona, se asume z02 = z01.

    Returns
    -------
    TT : Symbolic Matrix
        Matriz ABCD en función de los parámetros imagen.

    Raises
    ------
    ValueError
        Si z01 no es un símbolo o no es un número real positivo.
        Si z02 no es un símbolo o no es un número real positivo.
        Si gamma no es un número complejo.


    See Also
    --------
    :func:`Y2Tabcd_s`
    :func:`Tabcd2Z_s`
    :func:`Model_conversion`


    Examples
    --------
    >>> gamma = sp.symbols('gamma')
    >>> z01 = sp.symbols('z01')
    >>> z02 = sp.symbols('z02')
    >>> TT = I2Tabcd_s(gamma, z01, z02)
    >>> print(TT)
    Matrix([[cosh(gamma)*sqrt(z01/z02), sinh(gamma)*sqrt(z01*z02)], [sinh(gamma)/sqrt(z01*z02), cosh(gamma)*sqrt(z02/z01)]])

    Notes
    -----
    - Esta función está diseñada para trabajar con expresiones simbólicas utilizando el módulo SymPy.

    '''
    # Verificar que gamma sea un número complejo
    if not isinstance(gamma, sp.Expr) or not gamma.is_complex(): 
        raise ValueError("gamma debe ser un número complejo")

    # Verificar que z01 sea un símbolo o un número real positivo
    if not (isinstance(z01, sp.Symbol) or (isinstance(z01, int) and z01 > 0) or (isinstance(z01, float) and z01 > 0)):
        raise ValueError("z01 debe ser un símbolo o un número real positivo")

    # Verificar si z02 es proporcionado y, de serlo, que sea un símbolo o un número real positivo
    if z02 is not None:
        if not (isinstance(z02, sp.Symbol) or (isinstance(z02, int) and z02 > 0) or (isinstance(z02, float) and z02 > 0)):
            raise ValueError("z02 debe ser un símbolo o un número real positivo")

    # Si z02 no es proporcionado, se asume z02 = z01
    if z02 is None:
        z02 = z01

    # Construcción de la matriz de parámetros ABCD
    TT = sp.Matrix([[sp.cosh(gamma) * sp.sqrt(z01 / z02), sp.sinh(gamma) * sp.sqrt(z01 * z02)],
                    [sp.sinh(gamma) / sp.sqrt(z01 * z02), sp.cosh(gamma) * sp.sqrt(z02 / z01)]])

    return TT

def Model_conversion(src_model, dst_model):
    '''
    Convierte modelos de cuadripolos lineales de un formato a otro.

    Parameters
    ----------
    src_model : dict
        Diccionario que describe el modelo de origen.
        Debe tener las claves: 
            - 'model_name': nombre del modelo ('Z', 'T', etc.).
            - 'matrix': matriz de parámetros del modelo.
            - 'dep_var': variables dependientes del modelo.
            - 'indep_var': variables independientes del modelo.
            - (opcional) 'proxy_matrix': matriz de parámetros proxy en caso de ser necesario.
            - (opcional) 'neg_i2_current': indicador booleano si la corriente i2 se define con signo negativo.

    dst_model : dict
        Diccionario que describe el modelo de salida.
        Debe tener las mismas claves que src_model.

    Returns
    -------
    dict
        Diccionario que contiene la matriz convertida y el nombre del modelo resultante.

    Raises
    ------
    ValueError
        Si los modelos de origen y destino son iguales.
        Si falta alguna clave en src_model o dst_model.
        Si la variable independiente no es un símbolo o no es un número real positivo.


    See Also
    --------
    :func:`Y2Tabcd_s`
    :func:`Tabcd2Z_s`
    :func:`S2Ts_s`


    Example
    -------
    >>> import sympy as sp
    >>> # Parámetros Z (impedancia - circuito abierto)
    >>> ZZ = sp.Matrix([[z11, z12], [z21, z22]])
    >>> # Variables dependientes
    >>> vv = sp.Matrix([[v1], [v2]])
    >>> # Variables independientes
    >>> ii = sp.Matrix([[i1], [i2]])
    >>> # Parámetros Tdcba (transmisión inversa, DCBA)
    >>> TTi = sp.Matrix([[Ai, Bi], [-Ci, -Di]])
    >>> # Variables dependientes
    >>> ti_dep = sp.Matrix([[v2], [i2]])
    >>> # Variables independientes. (Signo negativo de corriente)
    >>> ti_ind = sp.Matrix([[v1], [i1]])
    >>> # Diccionario con la definición de cada modelo
    >>> src_model = {'model_name': 'Z', 'matrix': ZZ, 'dep_var': vv, 'indep_var': ii}
    >>> dst_model = {'model_name': 'T', 'matrix': TTi, 'dep_var': ti_dep, 'indep_var': ti_ind, 'neg_i2_current': True}
    >>> T_z = Model_conversion(src_model, dst_model)
    >>> print(T_z['matrix'])
    Matrix([[cosh(gamma)*sqrt(z01/z02), sinh(gamma)*sqrt(z01*z02)], [sinh(gamma)/sqrt(z01*z02), cosh(gamma)*sqrt(z02/z01)]])

    Notes
    -----
    - Esta función está diseñada para trabajar con expresiones simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si los modelos de origen y destino son iguales
    if src_model['model_name'] == dst_model['model_name']:
        raise ValueError("Los modelos de origen y destino son iguales")

    # Verificar que src_model tenga las claves necesarias
    required_keys = ['model_name', 'matrix', 'dep_var', 'indep_var']
    for key in required_keys:
        if key not in src_model:
            raise ValueError(f"Falta la clave '{key}' en src_model")

    # Verificar que dst_model tenga las claves necesarias
    for key in required_keys:
        if key not in dst_model:
            raise ValueError(f"Falta la clave '{key}' en dst_model")

    # Verificar que las variables independientes sean símbolos o números reales positivos
    for var in src_model['indep_var']:
        if not (isinstance(var, sp.Symbol) or (isinstance(var, int) and var > 0) or (isinstance(var, float) and var > 0)):
            raise ValueError("La variable independiente debe ser un símbolo o un número real positivo")

    # Si 'proxy_matrix' está presente en src_model, usarla como src_matrix; de lo contrario, usar src_model['matrix']
    src_matrix = src_model['proxy_matrix'] if 'proxy_matrix' in src_model else src_model['matrix']

    # Resolver para las variables dependientes de dst_model
    aa = sp.solve([src_matrix * src_model['indep_var'] - src_model['dep_var']], dst_model['dep_var'])

    # Reemplazar el determinante por Delta
    det_src_matrix = sp.det(src_matrix)
    if 'neg_i2_current' in src_model:
        det_src_matrix = -det_src_matrix

    dd = sp.Symbol('\Delta')
    QQ = sp.Matrix([[0, 0], [0, 0]])

    for jj, dep_var in enumerate(dst_model['dep_var'], start=1):
        yyy = sp.collect(sp.expand(aa[dep_var]), dst_model['indep_var'][0])
        yyy = sp.collect(yyy, dst_model['indep_var'][1])
        if dep_var.name == 'i2' and 'neg_i2_current' in dst_model:
            yyy = -yyy
        for kk, indep_var in enumerate(dst_model['indep_var'], start=1):
            bb = sp.cancel(yyy.coeff(indep_var, 1))
            if indep_var.name == 'i2' and 'neg_i2_current' in dst_model:
                bb = -bb
            QQ[jj - 1, kk - 1] = bb.subs(det_src_matrix, dd)

    return {'matrix': QQ, 'name': f"{dst_model['model_name']}_{src_model['model_name']}"}

#%%
   ############################################################################
  ## (NUMERICAS) Funciones de conversión de matrices de cuadripolos lineales #
 ############################################################################
#%%

def Y2Tabcd(YY):
    """
    Convierte una matriz de admitancia de dos puertos (YY) a la matriz de parámetros ABCD (TT).

    Parameters
    ----------
    YY : numpy.ndarray
        Matriz de admitancia de dos puertos.

    Returns
    -------
    TT : numpy.ndarray
        Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si YY no es una matriz de 2x2.
        Si Y21 es cero.


    See Also
    --------
    :func:`Z2Tabcd`
    :func:`Tabcd2Y`
    :func:`y2mai`


    Example
    -------
    >>> YY = np.array([[Y11, Y12], [Y21, Y22]])
    >>> TT = Y2Tabcd(YY)
    >>> print(TT)
    [[-Y22/Y21 -1/Y21]
     [-D/Y21 -Y22/Y21]]

    Notes
    -----
    - Esta función asume que YY tiene el formato [ [Y11, Y12], [Y21, Y22] ].
    - YY[1, 0] no puede ser cero para evitar una división por cero.

    """
    # Verificar que YY sea una matriz de 2x2
    if YY.shape != (2, 2):
        raise ValueError("YY debe ser una matriz de 2x2")

    # Verificar que Y21 no sea cero para evitar división por cero
    if YY[1, 0] == 0:
        raise ValueError("Y21 no puede ser cero")

    # Inicializar la matriz de parámetros ABCD
    TT = np.zeros_like(YY)

    # Calcular los elementos de la matriz TT
    TT[0, 0] = -YY[1, 1] / YY[1, 0]
    TT[0, 1] = -1 / YY[1, 0]
    TT[1, 0] = -np.linalg.det(YY) / YY[1, 0]
    TT[1, 1] = -YY[1, 1] / YY[1, 0]

    return TT

def Z2Tabcd(ZZ):
    '''
    Convierte una matriz de impedancia de dos puertos (ZZ) a la matriz de parámetros ABCD (TT).

    Parameters
    ----------
    ZZ : numpy.ndarray
        Matriz de impedancia de dos puertos.

    Returns
    -------
    TT : numpy.ndarray
        Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si ZZ no es una matriz de 2x2.
        Si Z21 es cero.


    See Also
    --------
    :func:`Y2Tabcd`
    :func:`Tabcd2Z`
    :func:`may2y`


    Example
    -------
    >>> ZZ = np.array([[Z11, Z12], [Z21, Z22]])
    >>> TT = Z2Tabcd(ZZ)
    >>> print(TT)
    [[Z11/Z21 DT/Z21]
     [1/Z21 Z22/Z21]]

    Notes
    -----
    - Esta función asume que ZZ tiene el formato [ [Z11, Z12], [Z21, Z22] ].
    - ZZ[1, 0] no puede ser cero para evitar una división por cero.

    '''
    # Verificar que ZZ sea una matriz de 2x2
    if ZZ.shape != (2, 2):
        raise ValueError("ZZ debe ser una matriz de 2x2")

    # Verificar que Z21 no sea cero para evitar división por cero
    if ZZ[1, 0] == 0:
        raise ValueError("Z21 no puede ser cero")

    # Inicializar la matriz de parámetros ABCD
    TT = np.zeros_like(ZZ)

    # Calcular los elementos de la matriz TT
    TT[0, 0] = ZZ[0, 0] / ZZ[1, 0]
    TT[0, 1] = np.linalg.det(ZZ) / ZZ[1, 0]
    TT[1, 0] = 1 / ZZ[1, 0]
    TT[1, 1] = ZZ[1, 1] / ZZ[1, 0]

    return TT

def Tabcd2Y(TT):
    '''
    Convierte una matriz de parámetros ABCD (TT) a la matriz de admitancia de dos puertos (YY).

    Parameters
    ----------
    TT : numpy.ndarray
        Matriz de parámetros ABCD.

    Returns
    -------
    YY : numpy.ndarray
        Matriz de admitancia de dos puertos.

    Raises
    ------
    ValueError
        Si TT no es una matriz de 2x2.
        Si B es cero.


    See Also
    --------
    :func:`Y2Tabcd`
    :func:`Tabcd2Z`
    :func:`may2y`


    Example
    -------
    >>> TT = np.array([[A, B], [C, D]])
    >>> YY = Tabcd2Y(TT)
    >>> print(YY)
    [[D/B -DT/B]
     [-1/B A/B]]

    Notes
    -----
    - Esta función asume que TT tiene el formato [ [A, B], [C, D] ].
    - B no puede ser cero para evitar una división por cero.

    '''
    # Verificar que TT sea una matriz de 2x2
    if TT.shape != (2, 2):
        raise ValueError("TT debe ser una matriz de 2x2")

    # Verificar que B no sea cero para evitar división por cero
    if TT[0, 1] == 0:
        raise ValueError("B no puede ser cero")

    # Inicializar la matriz de admitancia YY
    YY = np.zeros_like(TT)

    # Calcular los elementos de la matriz YY
    YY[0, 0] = TT[1, 1] / TT[0, 1]
    YY[0, 1] = -np.linalg.det(TT) / TT[0, 1]
    YY[1, 0] = -1 / TT[0, 1]
    YY[1, 1] = TT[0, 0] / TT[0, 1]

    return YY

def y2mai(YY):
    '''
    Convierte una matriz de admitancia indefinida (YY) a una matriz admitancia indefinida (Ymai).

    Parameters
    ----------
    YY : sympy.Matrix
        Matriz admitancia indefinida.

    Returns
    -------
    Ymai : sympy.Matrix
        Matriz admitancia indefinida.

    Raises
    ------
    ValueError
        Si YY no es una instancia de sympy.Matrix.


    See Also
    --------
    :func:`may2y`
    :func:`Y2Tabcd`
    :func:`I2Tabcd`


    Example
    -------
    >>> YY = sp.Matrix([[Y11, Y12], [Y21, Y22]])
    >>> Ymai = y2mai(YY)
    >>> print(Ymai)
    [[Y11 Y12 -Y11 - Y12]
     [Y21 Y22 -Y21 - Y22]
     [-Y11 - Y21 Y11 + Y21]
     [-Y12 - Y22 Y12 + Y22]]

    Notes
    -----
    - Esta función suma las corrientes de entrada y salida para obtener la matriz admitancia indefinida.
    - Se espera que YY sea una instancia de sympy.Matrix.

    '''
    # Verificar si YY es una instancia de sympy.Matrix
    if not isinstance(YY, sp.Matrix):
        raise ValueError("YY debe ser una instancia de sympy.Matrix")

    # Insertar filas y columnas para sumar las corrientes de entrada y salida
    Ymai = YY.row_insert(YY.shape[0], sp.Matrix([-sum(YY[:, ii]) for ii in range(YY.shape[1])]).transpose())
    Ymai = Ymai.col_insert(Ymai.shape[1], sp.Matrix([-sum(Ymai[ii, :]) for ii in range(Ymai.shape[0])]))
    Ymai[-1] = sum(YY)

    return Ymai

def may2y(Ymai, nodes2del):
    '''
    Convierte una matriz admitancia indefinida (Ymai) a una matriz admitancia (YY) luego de eliminar filas y columnas indicadas en nodes2del.

    Parameters
    ----------
    Ymai : sympy.Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Índices de las filas y columnas que se eliminarán.

    Returns
    -------
    YY : sympy.Matrix
        Matriz admitancia.

    Raises
    ------
    ValueError
        Si Ymai no es una instancia de sympy.Matrix.
        Si nodes2del no es una lista o un entero.
        Si los elementos de nodes2del no son enteros o están fuera del rango de índices de Ymai.


    See Also
    --------
    :func:`y2mai`
    :func:`Y2Tabcd`
    :func:`I2Tabcd`


    Example
    -------
    >>> Ymai = sp.Matrix([[Y11, Y12, Y13], [Y21, Y22, Y23], [Y31, Y32, Y33]])
    >>> nodes2del = [0, 2]
    >>> YY = may2y(Ymai, nodes2del)
    >>> print(YY)
    [[Y22]]

    Notes
    -----
    - Esta función elimina las filas y columnas indicadas en nodes2del de Ymai para obtener la matriz admitancia YY.
    - Se espera que Ymai sea una instancia de sympy.Matrix.
    - nodes2del puede ser una lista de índices o un solo entero.
    - Los índices en nodes2del deben ser enteros y estar dentro del rango de índices de Ymai.

    '''
    # Verificar si Ymai es una instancia de sympy.Matrix
    if not isinstance(Ymai, sp.Matrix):
        raise ValueError("Ymai debe ser una instancia de sympy.Matrix")

    # Verificar si nodes2del es una lista o un entero
    if not isinstance(nodes2del, list) and not isinstance(nodes2del, int):
        raise ValueError("nodes2del debe ser una lista o un entero")

    # Convertir nodes2del a lista si es un entero
    if isinstance(nodes2del, int):
        nodes2del = [nodes2del]

    # Verificar si los elementos de nodes2del son enteros
    if not all(isinstance(node, int) for node in nodes2del):
        raise ValueError("Los elementos de nodes2del deben ser enteros")

    # Verificar si los elementos de nodes2del están dentro del rango de índices de Ymai
    if not all(0 <= node < Ymai.rows for node in nodes2del):
        raise ValueError("Los elementos de nodes2del están fuera del rango de índices de Ymai")

    # Eliminar las filas y columnas indicadas en nodes2del
    YY = Ymai.copy()
    for node in sorted(nodes2del, reverse=True):
        YY.row_del(node)
        YY.col_del(node)

    return YY

def I2Tabcd(gamma, z01, z02=None):
    '''
    Convierte una ganancia compleja expresada en neppers (gamma) 
    y la impedancia de referencia (z01,2) en una matriz de parámetros ABCD (TT).

    Parameters
    ----------
    gamma : float or complex
        Ganancia compleja expresada en neppers (Re{gamma}) y radianes (Im{gamma}).
        
    z01 : float
        Impedancia de referencia del puerto 1.

    z02 : float, opcional
        Impedancia de referencia del puerto 2. Si no se proporciona, se asume z02 = z01.

    Returns
    -------
    TT : numpy.ndarray
        Matriz ABCD en función de los parámetros imagen.

    Raises
    ------
    ValueError
        Si z01 no es un número real positivo.
        Si z02 no es un número real positivo.


    See Also
    --------
    :func:`y2mai`
    :func:`Tabcd2Y`
    :func:`Y2Tabcd`


    Examples
    --------
    >>> gamma = 0.5 + 1j
    >>> z01 = 50
    >>> z02 = 75
    >>> TT = I2Tabcd(gamma, z01, z02)
    >>> print(TT)
    [[ 1.2039728+0.82171471j  0.82171471+0.30545599j]
     [ 0.30545599+0.82171471j  1.2039728+0.82171471j]]

    Notes
    -----
    - Esta función calcula la matriz de parámetros ABCD en función de una ganancia compleja gamma y las impedancias de referencia z01 y z02.
    - Si z02 no se proporciona, se asume que z02 = z01.
    - Se espera que z01 y z02 sean números reales positivos.

    '''
    # Verificar si z01 es un número real positivo
    if not isinstance(z01, (int, float)) or z01 <= 0:
        raise ValueError("z01 debe ser un número real positivo")

    # Verificar si z02 es proporcionado y, de serlo, que sea un número real positivo
    if z02 is not None:
        if not isinstance(z02, (int, float)) or z02 <= 0:
            raise ValueError("z02 debe ser un número real positivo")

    # Si z02 no es proporcionado, se asume z02 = z01
    if z02 is None:
        z02 = z01

    # Construcción de la matriz de parámetros ABCD
    TT = np.array([[np.cosh(gamma) * np.sqrt(z01 / z02), np.sinh(gamma) * np.sqrt(z01 * z02)],
                   [np.sinh(gamma) / np.sqrt(z01 * z02), np.cosh(gamma) * np.sqrt(z02 / z01)]])

    return TT

#%%
#%%
   ##############################################################
  ## (Simbólicas) Parámetros de cuadripolo para redes estandar #
 ##############################################################
#%%

def SparZ_s(Zexc, Z01=sp.Rational(1), Z02=sp.Rational(1)):
    '''
    Convierte una matriz de transferencia de scattering (Ts) simbólica 
    al modelo de parámetros scattering (S).

    Parameters
    ----------
    Zexc : sympy.Symbol
           Función de excitación de la impedancia a representar.
    Z01 : sympy.Symbol, optional
          Impedancia de referencia en el plano 1. Por defecto es 1.
    Z02 : sympy.Symbol, optional
          Impedancia de referencia en el plano 2. Por defecto es 1.

    Returns
    -------
    Spar : sympy.Matrix
           Matriz de parámetros de scattering de Z.

    Raises
    ------
    ValueError
        Si Zexc no es una instancia de Symbolic.
        Si Z01 no es una instancia de Symbolic.
        Si Z02 no es una instancia de Symbolic.


    See Also
    --------
    :func:`SparY_s`
    :func:`TabcdLYZ_s`
    :func:`TabcdZ_s`


    Examples
    --------
    >>> Zexc = sp.symbols('Zexc')
    >>> Z01 = sp.symbols('Z01')
    >>> Z02 = sp.symbols('Z02')
    >>> Spar = SparZ_s(Zexc, Z01, Z02)
    >>> print(Spar)
    Matrix([[Zexc, 2*Z01], [2*Z01, Zexc]])

    Notes
    -----
    - Esta función está diseñada para trabajar con impedancias simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si Zexc, Z01 y Z02 son instancias de Symbolic
    if not isinstance(Zexc, sp.Symbol):
        raise ValueError("Zexc debe ser una instancia de Symbolic")
    if not isinstance(Z01, sp.Symbol):
        raise ValueError("Z01 debe ser una instancia de Symbolic")
    if not isinstance(Z02, sp.Symbol):
        raise ValueError("Z02 debe ser una instancia de Symbolic")

    # Inicialización de la matriz de parámetros de scattering
    Spar = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Spar
    common = Zexc + sp.Rational(2) * Z01
    Spar[0, 0] = Zexc
    Spar[0, 1] = sp.Rational(2) * Z01
    Spar[1, 0] = sp.Rational(2) * Z01
    Spar[1, 1] = Zexc

    return sp.simplify(sp.expand(1 / common * Spar))

def SparY_s(Yexc, Y01=sp.Rational('1'), Y02=sp.Rational('1')):
    '''
    Convierte una matriz de transferencia de scattering (Ts) simbólica 
    al modelo de parámetros scattering (S).

    Parameters
    ----------
    Yexc : Symbolic impedance
           Función de excitación de la admitancia a representar.

    Y01 : Symbolic impedance, optional
          Admitancia de referencia en el plano 1. Por defecto es 1.

    Y02 : Symbolic impedance, optional
          Admitancia de referencia en el plano 2. Por defecto es 1.

    Returns
    -------
    Spar : Symbolic Matrix
           Matriz de parámetros de scattering de Y.

    Raises
    ------
    ValueError
        Si Yexc no es una instancia de Symbolic.
        Si Y01 no es una instancia de Symbolic.
        Si Y02 no es una instancia de Symbolic.


    See Also
    --------
    :func:`SparZ_s`
    :func:`TabcdLYZ_s`
    :func:`TabcdLZY`


    Examples
    --------
    >>> Zexc = sp.symbols('Yexc')
    >>> Y01 = sp.symbols('Y01')
    >>> Y02 = sp.symbols('Y02')
    >>> SparY = SparZ_s(Yexc, Y01, Y02)
    >>> print(SparY)
    Matrix([[Zexc, 2*Z01], [2*Z01, Zexc]])

    Notes
    -----
    - Esta función está diseñada para trabajar con admitancias simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si Yexc, Y01 y Y02 son instancias de Symbolic
    if not isinstance(Yexc, sp.Symbol):
        raise ValueError("Yexc debe ser una instancia de Symbolic")
    if not isinstance(Y01, sp.Symbol):
        raise ValueError("Y01 debe ser una instancia de Symbolic")
    if not isinstance(Y02, sp.Symbol):
        raise ValueError("Y02 debe ser una instancia de Symbolic")

    # Inicialización de la matriz de parámetros de scattering
    Spar = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Spar
    Spar[0, 0] = -Yexc
    Spar[0, 1] = sp.Rational(2) * Y01
    Spar[1, 0] = sp.Rational(2) * Y01
    Spar[1, 1] = -Yexc
    
    # Normalización por el término común
    common = Yexc + sp.Rational(2) * Y01
    return sp.simplify(sp.expand(1 / common * Spar))

def TabcdLYZ_s(Yexc, Zexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido por una Y en derivación seguida  por 
    una Z en serie.

    Parameters
    ----------
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.
    
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.

    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si Yexc no es una instancia de Symbolic.
        Si Zexc no es una instancia de Symbolic.


    See Also
    --------
    :func:`SparZ_s`
    :func:`TabcdZ`
    :func:`TabcdLZY`


    '''
    # Verificar si Yexc y Zexc son instancias de Symbolic
    if not isinstance(Yexc, sp.Symbol):
        raise ValueError("Yexc debe ser una instancia de Symbolic")
    if not isinstance(Zexc, sp.Symbol):
        raise ValueError("Zexc debe ser una instancia de Symbolic")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = sp.Rational('1')  
    Tpar[0, 1] = Zexc
    Tpar[1, 0] = Yexc
    Tpar[1, 1] = sp.Rational('1') + sp.simplify(sp.expand(Zexc * Yexc))
    
    return Tpar

def TabcdLZY_s(Zexc, Yexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido por una Z en serie seguida de una Y en 
    derivación.

    Parameters
    ----------
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.
    
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.

    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si Zexc no es una instancia de Symbolic.
        Si Yexc no es una instancia de Symbolic.


    See Also
    --------
    :func:`SparZ_s`
    :func:`TabcdLYZ_s`
    :func:`TabcdY_s`

    '''
    # Verificar si Zexc y Yexc son instancias de Symbolic
    if not isinstance(Zexc, sp.Symbol):
        raise ValueError("Zexc debe ser una instancia de Symbolic")
    if not isinstance(Yexc, sp.Symbol):
        raise ValueError("Yexc debe ser una instancia de Symbolic")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = sp.Rational('1') + sp.simplify(sp.expand(Zexc * Yexc))
    Tpar[0, 1] = Zexc
    Tpar[1, 0] = Yexc
    Tpar[1, 1] = sp.Rational('1')  
    
    return Tpar

def TabcdZ_s(Zexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido únicamente por una Z en serie.

    Parameters
    ----------
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.


    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si Zexc no es una instancia de Symbolic.


    See Also
    --------
    :func:`SparZ_s`
    :func:`TabcdLYZ_s`
    :func:`TabcdY_s`

    '''
    # Verificar si Zexc es una instancia de Symbolic
    if not isinstance(Zexc, sp.Symbol):
        raise ValueError("Zexc debe ser una instancia de Symbolic")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = sp.Rational('1') 
    Tpar[0, 1] = Zexc
    Tpar[1, 0] = sp.Rational('0') 
    Tpar[1, 1] = sp.Rational('1')  
    
    return Tpar

def TabcdY_s(Yexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido únicamente por una Y en derivación.

    Parameters
    ----------
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.

    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si Yexc no es una instancia de Symbolic.

    '''
    # Verificar si Yexc es una instancia de Symbolic
    if not isinstance(Yexc, sp.Symbol):
        raise ValueError("Yexc debe ser una instancia de Symbolic")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = sp.Rational('1') 
    Tpar[0, 1] = sp.Rational('0')
    Tpar[1, 0] = Yexc
    Tpar[1, 1] = sp.Rational('1')  
    
    return Tpar

#%%
   #############################################################
  ## (Numéricas) Parámetros de cuadripolo para redes estandar #
 #############################################################
#%%

def TabcdLYZ(Yexc, Zexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido por una Y en derivación seguida  por 
    una Z en serie.

    Parameters
    ----------
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.
    
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.

    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si Yexc no es una instancia de Symbolic.
        Si Zexc no es una instancia de Symbolic.

    '''
    # Verificar si Yexc y Zexc son instancias de Symbolic
    if not isinstance(Yexc, sp.Symbol):
        raise ValueError("Yexc debe ser una instancia de Symbolic")
    if not isinstance(Zexc, sp.Symbol):
        raise ValueError("Zexc debe ser una instancia de Symbolic")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = np.array([[0, 0], [0, 0]])
    
    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = 1 
    Tpar[0, 1] = Zexc
    Tpar[1, 0] = Yexc
    Tpar[1, 1] = 1 + Zexc * Yexc
    
    return Tpar

def TabcdLZY(Zexc, Yexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido por una Z en serie seguida una Y en 
    derivación.

    Parameters
    ----------
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.
    
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.

    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si Zexc no es una instancia de Symbolic.
        Si Yexc no es una instancia de Symbolic.

    '''
    # Verificar si Zexc y Yexc son instancias de Symbolic
    if not isinstance(Zexc, sp.Symbol):
        raise ValueError("Zexc debe ser una instancia de Symbolic")
    if not isinstance(Yexc, sp.Symbol):
        raise ValueError("Yexc debe ser una instancia de Symbolic")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = np.array([[0.0, 0], [0, 0]])
    
    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = 1 + Zexc * Yexc 
    Tpar[0, 1] = Zexc
    Tpar[1, 0] = Yexc
    Tpar[1, 1] = 1
    
    return Tpar

def TabcdZ(Zexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido únicamente por una Z en serie.

    Parameters
    ----------
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.


    Returns
    -------
    Tabcd : np.ndarray
           Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si Zexc no es una instancia de Symbolic.

    '''
    # Verificar si Zexc es una instancia de Symbolic
    if not isinstance(Zexc, sp.Symbol):
        raise ValueError("Zexc debe ser una instancia de Symbolic")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = np.array([[0.0, 0.0], [0, 0]])
    
    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = 1 
    Tpar[0, 1] = Zexc
    Tpar[1, 0] = 0
    Tpar[1, 1] = 1

    return Tpar

def TabcdY(Yexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido únicamente por una Y en derivación.

    Parameters
    ----------
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.


    Returns
    -------
    Tabcd : np.ndarray
           Matriz de parámetros ABCD.

    Raises
    ------
    ValueError
        Si Yexc no es una instancia de Symbolic.

    '''
    # Verificar si Yexc es una instancia de Symbolic
    if not isinstance(Yexc, sp.Symbol):
        raise ValueError("Yexc debe ser una instancia de Symbolic")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = np.array([[0.0, 0], [0, 0]])
    
    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = 1 
    Tpar[0, 1] = 0
    Tpar[1, 0] = Yexc
    Tpar[1, 1] = 1
    
    return Tpar

#%%
   ###############################################################
  ## Cálculo de transferencias e inmitancias a partir de la MAI #
 ###############################################################
#%%

def calc_MAI_ztransf_ij_mn(Ymai, ii=2, jj=3, mm=0, nn=1, verbose=False):
    """Calculates the impedance transfer V_ij / I_mn.

    This function calculates the impedance transfer V_ij / I_mn of a given
    multiport network represented by its admittance matrix.

    Parameters
    ----------
    Ymai : sp.Matrix
        The indefinite admittance matrix.
    ii : int, optional
        The index i of the output element, defaults to 2.
    jj : int, optional
        The index j of the output element, defaults to 3.
    mm : int, optional
        The index m of the input element, defaults to 0.
    nn : int, optional
        The index n of the input element, defaults to 1.
    verbose : bool, optional
        If True, prints intermediate calculations, defaults to False.

    Returns
    -------
    Tz : sp.Expr
        The impedance transfer.

    Raises
    ------
    ValueError
        If any of the indices is not an integer.
        If Ymai is not an instance of sp.Matrix.

    Examples
    --------
    >>> Ymai = sp.MatrixSymbol('Y', 4, 4)
    >>> calc_MAI_ztransf_ij_mn(Ymai)
    Tz

    """
    # Check if Ymai is an instance of sp.Matrix
    if not isinstance(Ymai, sp.Matrix):
        raise ValueError("Ymai must be an instance of sp.Matrix.")

    # Check if the indices are integers
    if not all(isinstance(val, int) for val in [ii, jj, mm, nn]):
        raise ValueError("Indices must be integers.")

    # Calculate cofactors
    num = Ymai.minor_submatrix(max(ii, jj), max(mm, nn)).minor_submatrix(min(ii, jj), min(mm, nn))
    den = Ymai.minor_submatrix(min(mm, nn), min(mm, nn))

    # Calculate determinants of the cofactors
    num_det = sp.simplify(num.det())
    den_det = sp.simplify(den.det())

    # Sign correction
    sign_correction = mm + nn + ii + jj
    Tz = sp.simplify(-1 ** sign_correction * num_det / den_det)

    # Print intermediate calculations if verbose is True
    if verbose:
        print("Intermediate calculations:")
        print(f"num: {num}, den: {den}, num_det: {num_det}, den_det: {den_det}")
        print(f"Tz: {Tz}")

    return Tz

def calc_MAI_vtransf_ij_mn(Ymai, ii=2, jj=3, mm=0, nn=1, verbose=False):
    """Calculates the voltage transfer V_ij / V_mn.

    This function calculates the voltage transfer V_ij / V_mn of a given
    multiport network represented by its admittance matrix.

    Parameters
    ----------
    Ymai : sp.Matrix
        The indefinite admittance matrix.
    ii : int, optional
        The index i of the output element, defaults to 2.
    jj : int, optional
        The index j of the output element, defaults to 3.
    mm : int, optional
        The index m of the input element, defaults to 0.
    nn : int, optional
        The index n of the input element, defaults to 1.
    verbose : bool, optional
        If True, prints intermediate calculations, defaults to False.

    Returns
    -------
    Av : sp.Expr
        The voltage transfer.

    Raises
    ------
    ValueError
        If any of the indices is not an integer.
        If Ymai is not an instance of sp.Matrix.

    Examples
    --------
    >>> Ymai = sp.MatrixSymbol('Y', 4, 4)
    >>> calc_MAI_vtransf_ij_mn(Ymai)
    Av

    """
    # Check if Ymai is an instance of sp.Matrix
    if not isinstance(Ymai, sp.Matrix):
        raise ValueError("Ymai must be an instance of sp.Matrix.")

    # Check if the indices are integers
    if not all(isinstance(val, int) for val in [ii, jj, mm, nn]):
        raise ValueError("Indices must be integers.")

    # Calculate cofactors
    num = Ymai.minor_submatrix(max(ii, jj), max(mm, nn)).minor_submatrix(min(ii, jj), min(mm, nn))
    den = Ymai.minor_submatrix(max(mm, nn), max(mm, nn)).minor_submatrix(min(mm, nn), min(mm, nn))

    # Calculate determinants of the cofactors
    num_det = sp.simplify(num.det())
    den_det = sp.simplify(den.det())

    # Sign correction
    sign_correction = mm + nn + ii + jj
    Av = sp.simplify(-1 ** sign_correction * num_det / den_det)

    # Print intermediate calculations if verbose is True
    if verbose:
        print("Intermediate calculations:")
        print(f"num: {num}, den: {den}, num_det: {num_det}, den_det: {den_det}")
        print(f"Av: {Av}")

    return Av

def calc_MAI_impedance_ij(Ymai, ii=0, jj=1, verbose=False):
    """Calculates the impedance transfer V_ij / V_mn.

    This function calculates the impedance transfer V_ij / V_mn of a given
    multiport network represented by its admittance matrix.

    Parameters
    ----------
    Ymai : sp.Matrix
        The indefinite admittance matrix.
    ii : int, optional
        The index i of the output element, defaults to 0.
    jj : int, optional
        The index j of the output element, defaults to 1.
    verbose : bool, optional
        If True, prints intermediate calculations, defaults to False.

    Returns
    -------
    ZZ : sp.Expr
        The impedance transfer.

    Raises
    ------
    ValueError
        If ii or jj is not an integer.
        If Ymai is not an instance of sp.Matrix.

    Examples
    --------
    >>> Ymai = sp.MatrixSymbol('Y', 4, 4)
    >>> calc_MAI_impedance_ij(Ymai)
    ZZ

    """
    # Check if Ymai is an instance of sp.Matrix
    if not isinstance(Ymai, sp.Matrix):
        raise ValueError("Ymai must be an instance of sp.Matrix.")

    # Check if ii and jj are integers
    if not isinstance(ii, int) or not isinstance(jj, int):
        raise ValueError("ii and jj must be integers.")

    # Calculate cofactor of second order
    num = Ymai.minor_submatrix(max(ii, jj), max(ii, jj)).minor_submatrix(min(ii, jj), min(ii, jj))
    # Any cofactor of first order
    den = Ymai.minor_submatrix(min(ii, jj), min(ii, jj))

    # Determinants of the cofactors
    num_det = sp.simplify(num.det())
    den_det = sp.simplify(den.det())

    # Calculate impedance transfer
    ZZ = sp.simplify(num_det / den_det)

    # Print intermediate calculations if verbose is True
    if verbose:
        print("Intermediate calculations:")
        print(f"num: {num}, den: {den}, num_det: {num_det}, den_det: {den_det}")
        print(f"ZZ: {ZZ}")

    return ZZ

