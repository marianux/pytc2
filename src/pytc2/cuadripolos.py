#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:14:13 2023

@author: mariano
"""

import numpy as np
import sympy as sp
from numbers import Real, Complex

import pandas as pd
import re as re

from pytc2.general import s, print_console_alert

import platform
import subprocess
import os


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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import S2Ts_s
    >>> Spar = sp.Matrix([[sp.symbols('S11'), sp.symbols('S12')],
    ...                   [sp.symbols('S21'), sp.symbols('S22')]])
    >>> Ts = S2Ts_s(Spar)
    >>> print(Ts)
    Matrix([[1/S21, -S22/S21], [S11/S21, -S11*S22/S21 + S12]])

    Notes
    -----
    - La matriz Spar debe tener la forma [ [Spar_11, Spar_12], [Spar_21, Spar_22] ].
    - Spar_12 no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    """
    # Verificar si Spar es una instancia de Symbolic Matrix
    if not isinstance(Spar, sp.MatrixBase):
        raise ValueError("Spar debe ser una instancia de Symbolic Matrix")
    
    # Verificar que Spar tenga el formato correcto
    if Spar.shape != (2, 2):
        raise ValueError("Spar debe tener el formato [ [Spar_11, Spar_12], [Spar_21, Spar_22] ]")
    
    # Verificar que Spar_12 no sea nulo
    if Spar[1, 0] == 0:
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import Ts2S_s
    >>> Ts = sp.Matrix([[sp.symbols('Ts11'), sp.symbols('Ts12')],
    ...                 [sp.symbols('Ts21'), sp.symbols('Ts22')]])
    >>> Spar = Ts2S_s(Ts)
    >>> print(Spar)
    Matrix([[Ts21/Ts11, Ts22 - Ts12*Ts21/Ts11], [1/Ts11, -Ts12/Ts11]])


    Notes
    -----
    - La matriz Ts debe tener la forma [ [Ts_11, Ts_12], [Ts_21, Ts_22] ].
    - Ts_11 no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    """
    # Verificar si Ts es una instancia de Symbolic Matrix
    if not isinstance(Ts, sp.MatrixBase):
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

def Ts2Tabcd_s(Ts, Z01=sp.Rational('1'), Z02=sp.Rational('1')):
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import Ts2Tabcd_s
    >>> Z0 = sp.symbols('Z0')
    >>> Ts = sp.Matrix([[sp.symbols('Ts_11'), sp.symbols('Ts_12')],
    ...                 [sp.symbols('Ts_21'), sp.symbols('Ts_22')]])
    >>> Tabcd = Ts2Tabcd_s(Ts, Z0)
    >>> print(Tabcd)
    Matrix([[Ts_11/2 - Ts_12/2 - Ts_21/2 + Ts_22/2, Z0*(Ts_11 - Ts_12 + Ts_21 - Ts_22)/2], [(Ts_11 + Ts_12 - Ts_21 - Ts_22)/(2*Z0), Ts_11/2 - Ts_12/2 - Ts_21/2 + Ts_22/2]])

    """
    # Check if Ts is an instance of sp.Matrix
    if not isinstance(Ts, sp.MatrixBase):
        raise ValueError("Ts must be an instance of sp.Matrix.")
    
    # Check if Z0 is an instance of sp.Expr
    if not isinstance(Z01, sp.Expr):
        raise ValueError("Z01 debe ser un número complejo (expresión simbolólica de SymPy)")

    if not isinstance(Z02, sp.Expr):
        raise ValueError("Z02 debe ser un número complejo (expresión simbolólica de SymPy)")

    # Convert Ts to S-parameter matrix and then to ABCD matrix
    return sp.simplify(sp.expand(S2Tabcd_s(Ts2S_s(Ts), Z01=Z01, Z02=Z02)))

def Tabcd2S_s(Tabcd, Z01=sp.Rational('1'), Z02=sp.Rational('1')):
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import Tabcd2S_s
    >>> Tabcd = sp.Matrix([[sp.symbols('A'), sp.symbols('B')],
    ...                    [sp.symbols('C'), sp.symbols('D')]])
    >>> Spar = Tabcd2S_s(Tabcd)
    >>> print(Spar)
    Matrix([[(A + B - C - D)/(A + B + C + D), 2*(A*D - B*C)/(A + B + C + D)], [2/(A + B + C + D), (-A + B - C + D)/(A + B + C + D)]])


    Notes
    -----
    - La matriz Tabcd debe tener el formato [ [A, B], [C, D] ].
    - La matriz Tabcd debe ser invertible para realizar la conversión correctamente.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si Tabcd es una instancia de Symbolic Matrix
    if not isinstance(Tabcd, sp.MatrixBase):
        raise ValueError("Tabcd debe ser una instancia de Symbolic Matrix")

    # Verificar que Tabcd tenga el formato correcto
    if Tabcd.shape != (2, 2):
        raise ValueError("Tabcd debe tener el formato [ [A, B], [C, D] ]")

    # # Verificar si la matriz Tabcd es invertible
    # if Tabcd.det() == 0:
    #     raise ValueError("La matriz Tabcd no es invertible")

    if not isinstance(Z01, sp.Expr):
        raise ValueError("Z01 debe ser un número complejo (expresión simbolólica de SymPy)")

    if not isinstance(Z02, sp.Expr):
        raise ValueError("Z02 debe ser un número complejo (expresión simbolólica de SymPy)")

    # Inicialización de la matriz de parámetros de scattering
    Spar = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Spar
    common = Tabcd[0, 0]*Z02 + Tabcd[0, 1] + Tabcd[1, 0]*Z01*Z02 + Tabcd[1, 1]*Z01
    Spar[0, 0] = Tabcd[0, 0] *Z02 + Tabcd[0, 1] - Tabcd[1, 0]*sp.conjugate(Z01)*Z02 - Tabcd[1, 1]*sp.conjugate(Z01)
    Spar[0, 1] = sp.Rational('2') * sp.sqrt(sp.re(Z01)*sp.re(Z02)) * sp.simplify(sp.expand(sp.Determinant(Tabcd))) 
    Spar[1, 0] = sp.Rational('2') * sp.sqrt(sp.re(Z01)*sp.re(Z02)) 
    Spar[1, 1] = -Tabcd[0, 0] *sp.conjugate(Z02) + Tabcd[0, 1] - Tabcd[1, 0]*sp.conjugate(Z02)*Z01 + Tabcd[1, 1]*Z01

    return sp.simplify(sp.expand(1 / common * Spar))

def S2Tabcd_s(Spar, Z01=sp.Rational('1'), Z02=sp.Rational('1')):
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import S2Tabcd_s
    >>> Spar = sp.Matrix([[sp.symbols('S11'), sp.symbols('S12')],
    ...                   [sp.symbols('S21'), sp.symbols('S22')]])
    >>> Tabcd = S2Tabcd_s(Spar)
    >>> print(Tabcd)
    Matrix([[(-S11*S22 - S11 + S12*S21 + S22 + 1)/(2*S21), (S11*S22 + S11 - S12*S21 + S22 + 1)/(2*S21)], [(S11*S22 - S11 - S12*S21 - S22 + 1)/(2*S21), (-S11*S22 - S11 + S12*S21 + S22 + 1)/(2*S21)]])


    Notes
    -----
    - La matriz Spar debe tener el formato [ [S11, S12], [S21, S22] ].
    - Spar[1, 0] no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si Spar es una instancia de Symbolic Matrix
    if not isinstance(Spar, sp.MatrixBase):
        raise ValueError("Spar debe ser una instancia de Symbolic Matrix")

    # Verificar que Spar tenga el formato correcto
    if Spar.shape != (2, 2):
        raise ValueError("Spar debe tener el formato [ [S11, S12], [S21, S22] ]")

    # Verificar que Spar[1, 0] no sea nulo
    if Spar[1, 0] == 0:
        raise ValueError("Spar[1, 0] no puede ser nulo")

    if not isinstance(Z01, sp.Expr):
        raise ValueError("Z01 debe ser un número complejo (expresión simbolólica de SymPy)")

    if not isinstance(Z02, sp.Expr):
        raise ValueError("Z02 debe ser un número complejo (expresión simbolólica de SymPy)")

    # Inicialización de la matriz de parámetros ABCD
    Tabcd = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Tabcd
    common = 2 * Spar[1, 0] * sp.sqrt(sp.re(Z01)*sp.re(Z02))
    Tabcd[0, 0] = (sp.conjugate(Z01) + Spar[0, 0] * Z01) * (1 - Spar[1, 1]) + Spar[1, 0] * Spar[0, 1] * Z01
    Tabcd[0, 1] = (sp.conjugate(Z01) + Spar[0, 0] * Z01) * (sp.conjugate(Z02) + Spar[1, 1] * Z02) - Spar[1, 0] * Spar[0, 1] * Z01 * Z02
    Tabcd[1, 0] = (1 - Spar[0, 0]) * (1 - Spar[1, 1]) - Spar[1, 0] * Spar[0, 1]
    Tabcd[1, 1] = (1 - Spar[0, 0]) * (sp.conjugate(Z02) + Spar[1, 1] * Z02) + Spar[1, 0] * Spar[0, 1] * Z02

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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import Y2Tabcd_s
    >>> YY = sp.Matrix([[sp.symbols('Y11'), sp.symbols('Y12')],
    ...                 [sp.symbols('Y21'), sp.symbols('Y22')]])
    >>> TT = Y2Tabcd_s(YY)
    >>> print(TT)
    Matrix([[-Y22/Y21, -1/Y21], [-(Y11*Y22 - Y12*Y21)/Y21, -Y22/Y21]])


    Notes
    -----
    - La matriz YY debe tener el formato [ [Y11, Y12], [Y21, Y22] ].
    - YY[1, 0] no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    """
    # Verificar si YY es una instancia de Symbolic Matrix
    if not isinstance(YY, sp.MatrixBase):
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import Z2Tabcd_s
    >>> ZZ = sp.Matrix([[sp.symbols('Z11'), sp.symbols('Z12')],
    ...                 [sp.symbols('Z21'), sp.symbols('Z22')]])
    >>> TT = Z2Tabcd_s(ZZ)
    >>> print(TT)
    Matrix([[Z11/Z21, (Z11*Z22 - Z12*Z21)/Z21], [1/Z21, Z22/Z21]])


    Notes
    -----
    - La matriz ZZ debe tener el formato [ [Z11, Z12], [Z21, Z22] ].
    - Z21 no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si ZZ es una instancia de Symbolic Matrix
    if not isinstance(ZZ, sp.MatrixBase):
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import Tabcd2Z_s
    >>> TT = sp.Matrix([[sp.symbols('A'), sp.symbols('B')],
    ...                 [sp.symbols('C'), sp.symbols('D')]])
    >>> ZZ = Tabcd2Z_s(TT)
    >>> print(ZZ)
    Matrix([[A/C, (A*D - B*C)/C], [1/C, D/C]])


    Notes
    -----
    - La matriz TT debe tener el formato [ [A, B], [C, D] ].
    - C no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si TT es una instancia de Symbolic Matrix
    if not isinstance(TT, sp.MatrixBase):
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import Tabcd2Y_s
    >>> TT = sp.Matrix([[sp.symbols('A'), sp.symbols('B')],
    ...                 [sp.symbols('C'), sp.symbols('D')]])
    >>> YY = Tabcd2Y_s(TT)
    >>> print(YY)
    Matrix([[D/B, -(A*D - B*C)/B], [-1/B, A/B]])

    Notes
    -----
    - La matriz TT debe tener el formato [ [A, B], [C, D] ].
    - B no puede ser nulo.
    - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si TT es una instancia de Symbolic Matrix
    if not isinstance(TT, sp.MatrixBase):
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import I2Tabcd_s
    >>> gamma = sp.symbols('gamma')
    >>> z01 = sp.symbols('z01')
    >>> z02 = sp.symbols('z02')
    >>> TT = I2Tabcd_s(gamma, z01, z02)
    >>> print(TT)
    Matrix([[sqrt(z01/z02)*cosh(gamma), sqrt(z01*z02)*sinh(gamma)], [sinh(gamma)/sqrt(z01*z02), sqrt(z02/z01)*cosh(gamma)]])


    Notes
    -----
    - Esta función está diseñada para trabajar con expresiones simbólicas utilizando el módulo SymPy.

    '''
    # Verificar que gamma sea un número complejo
    if not isinstance(gamma, sp.Expr): 
        raise ValueError("gamma debe ser un número complejo (expresión simbolólica de SymPy)")

    # Verificar que z01 sea un símbolo o un número real positivo
    if not isinstance(z01, sp.Expr):
        raise ValueError("z01 debe ser un real (expresión simbolólica de SymPy)")

    # Verificar si z02 es proporcionado y, de serlo, que sea un símbolo o un número real positivo
    if not isinstance(z02, (sp.Expr, type(None))):
        raise ValueError("z02 debe ser un real (expresión simbolólica de SymPy)")

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
        - 'proxy_matrix': (opcional) matriz de parámetros auxiliar. Por ejemplo para 
        relacionar modelos que no tengan variables en común (S->Z).
        Se necesitará una conversión intermedia, en PyTC2 se
        adopta :math:`T_{ABCD}` como modelo intermedio.                          
        - 'neg_i2_current': (opcional) indicador booleano si la corriente i2 se define con signo negativo.

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
    >>> from pytc2.cuadripolos import Model_conversion
    >>> v1, v2, i1, i2 = sp.symbols('v1, v2, i1, i2', complex=True)
    >>> z11, z12, z21, z22 = sp.symbols('z11, z12, z21, z22', complex=True)
    >>> Ai, Bi, Ci, Di = sp.symbols('Ai, Bi, Ci, Di', complex=True)
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
    Matrix([[z22/z12, -\Delta/z12], [-1/z12, z11/z12]])


    Notes
    -----
    - Esta función está diseñada para trabajar con expresiones simbólicas utilizando el módulo SymPy.

    '''
    # Verificar que src_model tenga las claves necesarias
    required_keys = ['model_name', 'matrix', 'dep_var', 'indep_var']
    for key in required_keys:
        if key not in src_model:
            raise ValueError(f"Falta la clave '{key}' en src_model")

    # Verificar que dst_model tenga las claves necesarias
    for key in required_keys:
        if key not in dst_model:
            raise ValueError(f"Falta la clave '{key}' en dst_model")

    # Verificar si los modelos de origen y destino son iguales
    if src_model['model_name'] == dst_model['model_name']:
        return {'matrix': sp.Matrix([[1,1],[1,1]]), 'name': f"{dst_model['model_name']}_{src_model['model_name']}"}


    # Verificar que las variables independientes sean símbolos ""o números reales positivos""
    for var in src_model['indep_var']:
        # if not (isinstance(var, sp.Expr) or (isinstance(var, int) and var > 0) or (isinstance(var, float) and var > 0)):
        if not isinstance(var, sp.Expr) :
            raise ValueError("La variable independiente debe ser un símbolo")

    for var in dst_model['indep_var']:
        # if not (isinstance(var, sp.Expr) or (isinstance(var, int) and var > 0) or (isinstance(var, float) and var > 0)):
        if not isinstance(var, sp.Expr) :
            raise ValueError("La variable independiente debe ser un símbolo")

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


def y2mai(YY):
    '''
    Convierte una matriz de admitancia definida (YY) a una matriz admitancia indefinida (Ymai).

    Parameters
    ----------
    YY : sympy.Matrix
        Matriz admitancia definida.

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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import y2mai
    >>> YY = sp.Matrix([[sp.symbols('Y11'), sp.symbols('Y12')],
    ...                 [sp.symbols('Y21'), sp.symbols('Y22')]])
    >>> Ymai = y2mai(YY)
    >>> print(Ymai)
    Matrix([[Y11, Y12, -Y11 - Y12], [Y21, Y22, -Y21 - Y22], [-Y11 - Y21, -Y12 - Y22, Y11 + Y12 + Y21 + Y22]])


    Notes
    -----
    - Esta función suma las corrientes de entrada y salida para obtener la matriz admitancia indefinida.
    - Se espera que YY sea una instancia de sympy.Matrix.

    '''
    # Verificar si YY es una instancia de sympy.Matrix
    if not isinstance(YY, sp.MatrixBase):
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import may2y
    >>> Ymai = sp.Matrix([[sp.symbols('Y11'), sp.symbols('Y12'), sp.symbols('Y13')],
    ...                 [sp.symbols('Y21'), sp.symbols('Y22'), sp.symbols('Y23')],
    ...                 [sp.symbols('Y31'), sp.symbols('Y32'), sp.symbols('Y33')]])
    >>> nodes2del = [0, 2]
    >>> YY = may2y(Ymai, nodes2del)
    >>> print(YY)
    Matrix([[Y22]])


    Notes
    -----
    - Esta función elimina las filas y columnas indicadas en nodes2del de Ymai para obtener la matriz admitancia YY.
    - Se espera que Ymai sea una instancia de sympy.Matrix.
    - nodes2del puede ser una lista de índices o un solo entero.
    - Los índices en nodes2del deben ser enteros y estar dentro del rango de índices de Ymai.

    '''
    # Verificar si Ymai es una instancia de sympy.Matrix
    if not isinstance(Ymai, sp.MatrixBase):
        raise ValueError("Ymai debe ser una instancia de sympy.Matrix")

    # Verificar si nodes2del es una lista o un entero
    if not isinstance(nodes2del, (list, int)) :
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
    >>> import numpy as np
    >>> from pytc2.cuadripolos import Y2Tabcd
    >>> YY = np.array([[6.0, -3.0], [-3.0, 5.0]])
    >>> TT = Y2Tabcd(YY)
    >>> print(TT)
    [[1.66666667 0.33333333]
     [7.         2. ]]
    
    >>> # Recordar la conversión entre modelos:
    [[-Y22/Y21 -1/Y21]
     [-D/Y21 -Y11/Y21]]
    
    
    Notes
    -----
    - Esta función asume que YY tiene el formato [ [Y11, Y12], [Y21, Y22] ].
    - YY[1, 0] no puede ser cero para evitar una división por cero.

    """
    
    if not isinstance(YY, np.ndarray):
        raise ValueError("YY debe ser una instancia de np.ndarray")
    
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
    TT[1, 1] = -YY[0, 0] / YY[1, 0]

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
    >>> import numpy as np
    >>> from pytc2.cuadripolos import Z2Tabcd
    >>> ZZ = np.array([[6., 3.], [3., 5.]])
    >>> TT = Z2Tabcd(ZZ)
    >>> print(TT)
    [[2.         7.        ]
     [0.33333333 1.66666667]]
    
    >>> # Recordar la conversión entre modelos:
    [[Z11/Z21 DT/Z21]
     [1/Z21 Z22/Z21]]
    

    Notes
    -----
    - Esta función asume que ZZ tiene el formato [ [Z11, Z12], [Z21, Z22] ].
    - ZZ[1, 0] no puede ser cero para evitar una división por cero.

    '''
    if not isinstance(ZZ, np.ndarray):
        raise ValueError("ZZ debe ser una instancia de np.ndarray")
    
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

def Tabcd2Z(TT):
    '''
    Convierte una matriz de parámetros ABCD (TT) a la matriz de impedancia de dos puertos (ZZ).

    Parameters
    ----------
    TT : numpy.ndarray
        Matriz de parámetros ABCD.

    Returns
    -------
    ZZ : numpy.ndarray
        Matriz de impedancia de dos puertos.

    Raises
    ------
    ValueError
        Si TT no es una matriz de 2x2.
        Si B es cero.


    See Also
    --------
    :func:`Z2Tabcd`
    :func:`Tabcd2Y`
    :func:`may2y`


    Example
    -------
    >>> import numpy as np
    >>> from pytc2.cuadripolos import Tabcd2Z
    >>> TT = np.array([[2, 7], [1./3., 5./3.]])
    >>> ZZ = Tabcd2Z(TT)
    >>> print(ZZ)
    [[ 6.,  3.]
     [3.,  5.]]    
    
    >>> # Recordar la conversión entre modelos:
     [[A/C DT/C]
     [1/C D/C]]

    Notes
    -----
    - Esta función asume que TT tiene el formato [ [A, B], [C, D] ].
    - C no puede ser cero para evitar una división por cero.

    '''
    if not isinstance(TT, np.ndarray):
        raise ValueError("TT debe ser una instancia de np.ndarray")
    
    # Verificar que TT sea una matriz de 2x2
    if TT.shape != (2, 2):
        raise ValueError("TT debe ser una matriz de 2x2")

    # Verificar que B no sea cero para evitar división por cero
    if TT[1, 0] == 0:
        raise ValueError("C no puede ser cero")

    # Inicializar la matriz de admitancia YY
    ZZ = np.zeros_like(TT)

    # Calcular los elementos de la matriz YY
    ZZ[0, 0] = TT[0, 0] 
    ZZ[0, 1] = np.linalg.det(TT) 
    ZZ[1, 0] = 1 
    ZZ[1, 1] = TT[1, 1]

    return 1/TT[1, 0] * ZZ

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
    >>> import numpy as np
    >>> from pytc2.cuadripolos import Tabcd2Y
    >>> TT = np.array([[5./3., 1./3.], [7., 2.]])
    >>> YY = Tabcd2Y(TT)
    >>> print(YY)
    [[ 6. -3.]
     [-3.  5.]]    
    
    >>> # Recordar la conversión entre modelos:
    [[D/B -DT/B]
     [-1/B A/B]]

    Notes
    -----
    - Esta función asume que TT tiene el formato [ [A, B], [C, D] ].
    - B no puede ser cero para evitar una división por cero.

    '''
    if not isinstance(TT, np.ndarray):
        raise ValueError("TT debe ser una instancia de np.ndarray")
    
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
    >>> import numpy as np
    >>> from pytc2.cuadripolos import I2Tabcd
    >>> gamma = 0.5 + 1.j
    >>> z01 = 50.
    >>> z02 = 75.
    >>> TT = I2Tabcd(gamma, z01, z02)
    >>> print(TT)
    [[4.97457816e-01+3.58022793e-01j 1.72412844e+01+5.81058484e+01j]
     [4.59767584e-03+1.54948929e-02j 7.46186724e-01+5.37034190e-01j]]

    
    >>> # Recordar la conversión entre modelos:
    TT = np.array([[np.cosh(gamma) * np.sqrt(z01 / z02), np.sinh(gamma) * np.sqrt(z01 * z02)],
                   [np.sinh(gamma) / np.sqrt(z01 * z02), np.cosh(gamma) * np.sqrt(z02 / z01)]])


    Notes
    -----
    - Esta función calcula la matriz de parámetros ABCD en función de una ganancia compleja gamma y las impedancias de referencia z01 y z02.
    - Si z02 no se proporciona, se asume que z02 = z01.
    - Se espera que z01 y z02 sean números reales positivos.

    '''
    
    # Verificar si gamma es un número complejo
    if not isinstance(gamma, Complex):
        raise ValueError("Gamma debe ser un número complejo")

    # Verificar si z01 es un número real positivo
    if not isinstance(z01, Real) or z01 <= 0:
        raise ValueError("z01 debe ser un número real positivo")


    # Verificar si z02 es proporcionado y, de serlo, que sea un número real positivo
    if not isinstance(z02, (Real, type(None))) or z02 <= 0:
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

def SparZ_s(Zexc, Z01=sp.Rational(1), Z02=None):
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import SparZ_s
    >>> Zexc = sp.symbols('Z')
    >>> Z01 = sp.symbols('Z01')
    >>> Z02 = sp.symbols('Z02')
    >>> Spar = SparZ_s(Zexc, Z01, Z01)
    >>> print(Spar)
    Matrix([[Z/(Z + 2*Z01), 2*Z01/(Z + 2*Z01)], [2*Z01/(Z + 2*Z01), Z/(Z + 2*Z01)]])
    
    >>> # Recordar la definición de los parámetros S de una Z en serie:
    1/(Z + 2*Z01) * [[Z,     2*Z01], 
                     [2*Z01, Z]])

    Notes
    -----
    - Esta función está diseñada para trabajar con impedancias simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si Zexc, Z01 y Z02 son instancias de Symbolic
    if not isinstance(Zexc, sp.Expr):
        raise ValueError("Zexc debe ser una instancia de Symbolic")
    if not isinstance(Z01, sp.Expr):
        raise ValueError("Z01 debe ser una instancia de Symbolic")
    if not isinstance(Z02, (sp.Expr, type(None))):
        raise ValueError("Z02 debe ser una instancia de Symbolic")

    if Z02 is None:
        Z02 = Z01

    # Inicialización de la matriz de parámetros de scattering
    Spar = sp.Matrix([[0, 0], [0, 0]])

    # Cálculo de los elementos de la matriz Spar
    common = Zexc + Z02 + Z01
    Spar[0, 0] = Zexc + Z02 - Z01
    Spar[0, 1] = sp.Rational(2) * Z01 * sp.sqrt(Z02/Z01)
    Spar[1, 0] = sp.Rational(2) * Z02 * sp.sqrt(Z01/Z02)
    Spar[1, 1] = Zexc + Z01 - Z02

    return sp.simplify(sp.expand(1 / common * Spar))

def SparY_s(Yexc, Y01=sp.Rational('1'), Y02=None):
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
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import SparY_s
    >>> Yexc = sp.symbols('Yexc')
    >>> Y01 = sp.symbols('Y01')
    >>> Y02 = sp.symbols('Y02')
    >>> SparY = SparY_s(Yexc, Y01)
    >>> print(SparY)
    Matrix([[-Yexc/(2*Y01 + Yexc), 2*Y01/(2*Y01 + Yexc)], [2*Y01/(2*Y01 + Yexc), -Yexc/(2*Y01 + Yexc)]])
    
    >>> # Recordar la definición de los parámetros S de una Y en derivación:
    1/(Y + 2*Y01) * [[-Y,     2*Y01], 
                     [2*Y01,  -Y]])

    Notes
    -----
    - Esta función está diseñada para trabajar con admitancias simbólicas utilizando el módulo SymPy.

    '''
    # Verificar si Yexc, Y01 y Y02 son instancias de Symbolic
    if not isinstance(Yexc, sp.Expr):
        raise ValueError("Yexc debe ser una instancia de Symbolic")
    if not isinstance(Y01, sp.Expr):
        raise ValueError("Y01 debe ser una instancia de Symbolic")
    if not isinstance(Y02, (sp.Expr, type(None))):
        raise ValueError("Y02 debe ser una instancia de Symbolic")

    if Y02 is None:
        Y02 = Y01

    # Inicialización de la matriz de parámetros de scattering
    Spar = sp.Matrix([[0, 0], [0, 0]])

    # Normalización por el término común
    common = Yexc + Y01 + Y02

    # Cálculo de los elementos de la matriz Spar
    Spar[0, 0] = Y01 - Yexc - Y02
    Spar[0, 1] = sp.Rational(2) * Y01 * sp.sqrt(Y01/Y02)
    Spar[1, 0] = sp.Rational(2) * Y02 * sp.sqrt(Y02/Y01)
    Spar[1, 1] = Y02 - Yexc - Y01
    
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


    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import TabcdLYZ_s
    >>> Y = sp.symbols('Y')
    >>> Z = sp.symbols('Z')
    >>> TT = TabcdLYZ_s(Y, Z)
    >>> print(TT)
    Matrix([[1, Z], [Y, Y*Z + 1]])

    '''
    # Verificar si Yexc y Zexc son instancias de Symbolic
    if not isinstance(Yexc, sp.Expr):
        raise ValueError("Yexc debe ser una instancia de Symbolic")
    if not isinstance(Zexc, sp.Expr):
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
    
    
    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import TabcdLZY_s
    >>> Y = sp.symbols('Y')
    >>> Z = sp.symbols('Z')
    >>> TT = TabcdLZY_s(Z, Y)
    >>> print(TT)
    Matrix([[Y*Z + 1, Z], [Y, 1]])
    

    '''
    # Verificar si Zexc y Yexc son instancias de Symbolic
    if not isinstance(Zexc, sp.Expr):
        raise ValueError("Zexc debe ser una instancia de Symbolic")
    if not isinstance(Yexc, sp.Expr):
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
    
    
    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import TabcdZ_s
    >>> Z = sp.symbols('Z')
    >>> TT = TabcdZ_s(Z)
    >>> print(TT)
    Matrix([[1, Z], [0, 1]])
    

    '''
    # Verificar si Zexc es una instancia de Symbolic
    if not isinstance(Zexc, sp.Expr):
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


    See Also
    --------
    :func:`SparZ_s`
    :func:`TabcdLYZ_s`
    :func:`TabcdY_s`
    
    
    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.cuadripolos import TabcdY_s
    >>> Y = sp.symbols('Y')
    >>> TT = TabcdY_s(Y)
    >>> print(TT)
    Matrix([[1, 0], [Y, 1]])
  
    '''
    # Verificar si Yexc es una instancia de Symbolic
    if not isinstance(Yexc, sp.Expr):
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


    Examples
    --------
    >>> from pytc2.cuadripolos import TabcdLYZ
    >>> TT = TabcdLYZ(Yexc=2., Zexc=3.)
    >>> print(TT)
    [[1 3]
     [2 7]]
    
    >>> # Recordar la definición de la matriz como:
    ([[1, Z], [Y, Y*Z + 1]])

    '''
    # Verificar si Zexc, Yexc es un número complejo
    if not isinstance(Zexc, Complex):
        raise ValueError("Zexc debe ser un número complejo")

    if not isinstance(Yexc, Complex):
        raise ValueError("Yexc debe ser un número complejo")


    # Inicialización de la matriz de parámetros ABCD
    Tpar = np.array([[0., 0.], [0., 0.]])
    
    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = 1. 
    Tpar[0, 1] = Zexc
    Tpar[1, 0] = Yexc
    Tpar[1, 1] = 1. + Zexc * Yexc
    
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


    Examples
    --------
    >>> from pytc2.cuadripolos import TabcdLZY
    >>> TT = TabcdLZY(Yexc=2., Zexc=3.)
    >>> print(TT)
    [[7. 3.]
     [2. 1.]]

    >>> # Recordar la definición de la matriz como:
    [[Y*Z + 1, Z], [Y, 1]]

    '''
    # Verificar si Zexc, Yexc es un número complejo
    if not isinstance(Zexc, Complex):
        raise ValueError("Zexc debe ser un número complejo")

    if not isinstance(Yexc, Complex):
        raise ValueError("Yexc debe ser un número complejo")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = np.array([[0., 0.], [0., 0.]])
    
    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = 1. + Zexc * Yexc 
    Tpar[0, 1] = Zexc
    Tpar[1, 0] = Yexc
    Tpar[1, 1] = 1.
    
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


    Examples
    --------
    >>> from pytc2.cuadripolos import TabcdZ
    >>> TT = TabcdZ(Zexc=3.)
    >>> print(TT)
    [[1. 3.]
     [0. 1.]]
    
    >>> # Recordar la definición de la matriz como:
    [[1, Z], 
     [0, 1]]

    '''
    # Verificar si Zexc, Yexc es un número complejo
    if not isinstance(Zexc, Complex):
        raise ValueError("Zexc debe ser un número complejo")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = np.array([[0.0, 0.0], [0., 0.]])
    
    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = 1. 
    Tpar[0, 1] = Zexc
    Tpar[1, 0] = 0.
    Tpar[1, 1] = 1.

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


    Examples
    --------
    >>> from pytc2.cuadripolos import TabcdY
    >>> TT = TabcdY(Yexc=2.)
    >>> print(TT)
    [[1. 0.]
     [2. 1.]]
    
    >>> # Recordar la definición de la matriz como:
    [[1, 0], 
     [Y, 1]]


    '''
    # Verificar si Zexc, Yexc es un número complejo
    if not isinstance(Yexc, Complex):
        raise ValueError("Yexc debe ser un número complejo")

    # Inicialización de la matriz de parámetros ABCD
    Tpar = np.array([[0., 0.], [0., 0.]])
    
    # Cálculo de los elementos de la matriz Tpar
    Tpar[0, 0] = 1.
    Tpar[0, 1] = 0.
    Tpar[1, 0] = Yexc
    Tpar[1, 1] = 1.
    
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
    >>> # Para la siguiente red eléctrica:
    >>> # Numeramos los polos de 0 a n=3
    >>> # 
    >>> #     0-------+--Y1----2---Y3--3---
    >>> #                      |           /
    >>> #                     Y2           / R
    >>> #                      |           /
    >>> #     1----------------+-------1----
    >>> #  
    >>> from pytc2.general import print_latex, a_equal_b_latex_s
    >>> from pytc2.cuadripolos import calc_MAI_ztransf_ij_mn
    >>> import sympy as sp
    >>> input_port = [0, 1]
    >>> output_port = [3, 1]
    >>> Y1, Y2, Y3 = sp.symbols('Y1 Y2 Y3', complex=True)
    >>> G = sp.symbols('G', real=True, positive=True)
    >>> #      Nodos: 0      1        2        3
    >>> Ymai = sp.Matrix([  
    >>>                  [ Y1,    0,      -Y1,      0],
    >>>                  [ 0,    Y2+G,    -Y2,     -G],
    >>>                  [ -Y1,  -Y2,    Y1+Y2+Y3, -Y3],
    >>>                  [ 0,    -G,      -Y3,      Y3+G ]
    >>>                  ])
    >>> s = sp.symbols('s ', complex=True)
    >>> # Butter de 3er orden doblemente cargado
    >>> Ymai = Ymai.subs(Y1, 1/s/sp.Rational('1'))
    >>> Ymai = Ymai.subs(Y3, 1/s/sp.Rational('1'))
    >>> Ymai = Ymai.subs(Y2, s*sp.Rational('2'))
    >>> # con_detalles = False
    >>> con_detalles = True
    >>> # Calculo la Z en el puerto de entrada a partir de la MAI
    >>> Zmai = calc_MAI_ztransf_ij_mn(Ymai, output_port[0], output_port[1], input_port[0], input_port[1], verbose=con_detalles)
    >>> print_latex(a_equal_b_latex_s('Z(s)', Zmai))
    Zmai = -1/(2*G*s**2 + G + 2*s)

    """
    # Check if Ymai is an instance of sp.Matrix
    if not isinstance(Ymai, sp.MatrixBase):
        raise ValueError("Ymai must be an instance of sp.Matrix.")

    # Check if the indices are integers
    if not all(isinstance(val, int) for val in [ii, jj, mm, nn]):
        raise ValueError("Indices must be integers.")
    # Check if verbose is an instance of bool
    if not isinstance(verbose, bool):
        raise ValueError("verbose must be an instance of bool")


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
    >>> # Para la siguiente red eléctrica:
    >>> # Numeramos los polos de 0 a n=3
    >>> # 
    >>> #     0-------+--Y1----2---Y3--3---
    >>> #                      |           /
    >>> #                     Y2           / R
    >>> #                      |           /
    >>> #     1----------------+-------1----
    >>> #  
    >>> from pytc2.general import print_latex, a_equal_b_latex_s
    >>> from pytc2.cuadripolos import calc_MAI_vtransf_ij_mn
    >>> import sympy as sp
    >>> input_port = [0, 1]
    >>> output_port = [3, 1]
    >>> Y1, Y2, Y3 = sp.symbols('Y1 Y2 Y3', complex=True)
    >>> G = sp.symbols('G', real=True, positive=True)
    >>> #      Nodos: 0      1        2        3
    >>> Ymai = sp.Matrix([  
    >>>                  [ Y1,    0,      -Y1,      0],
    >>>                  [ 0,    Y2+G,    -Y2,     -G],
    >>>                  [ -Y1,  -Y2,    Y1+Y2+Y3, -Y3],
    >>>                  [ 0,    -G,      -Y3,      Y3+G ]
    >>>                  ])
    >>> s = sp.symbols('s ', complex=True)
    >>> # Butter de 3er orden doblemente cargado
    >>> Ymai = Ymai.subs(Y1, 1/s/sp.Rational('1'))
    >>> Ymai = Ymai.subs(Y3, 1/s/sp.Rational('1'))
    >>> Ymai = Ymai.subs(Y2, s*sp.Rational('2'))
    >>> # con_detalles = False
    >>> con_detalles = True
    >>> # Calculo la Z en el puerto de entrada a partir de la MAI
    >>> Vmai = calc_MAI_vtransf_ij_mn(Ymai, output_port[0], output_port[1], input_port[0], input_port[1], verbose=con_detalles)
    >>> print_latex(a_equal_b_latex_s('T(s)', Vmai ))
    Vmai = -1/(2*G*s + 2*s**2*(G*s + 1) + 1)

    """
    # Check if Ymai is an instance of sp.Matrix
    if not isinstance(Ymai, sp.MatrixBase):
        raise ValueError("Ymai must be an instance of sp.Matrix.")

    # Check if the indices are integers
    if not all(isinstance(val, int) for val in [ii, jj, mm, nn]):
        raise ValueError("Indices must be integers.")

    # Check if verbose is an instance of bool
    if not isinstance(verbose, bool):
        raise ValueError("verbose must be an instance of bool")

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
    >>> # Para la siguiente red eléctrica:
    >>> # Numeramos los polos de 0 a n=3
    >>> # 
    >>> #     0-------+--Y1----2---Y3--3---
    >>> #                      |           /
    >>> #                     Y2           / R
    >>> #                      |           /
    >>> #     1----------------+-------1----
    >>> #  
    >>> from pytc2.general import print_latex, a_equal_b_latex_s
    >>> from pytc2.cuadripolos import calc_MAI_impedance_ij
    >>> import sympy as sp
    >>> input_port = [0, 1]
    >>> output_port = [3, 1]
    >>> Y1, Y2, Y3 = sp.symbols('Y1 Y2 Y3', complex=True)
    >>> G = sp.symbols('G', real=True, positive=True)
    >>> #      Nodos: 0      1        2        3
    >>> Ymai = sp.Matrix([  
    >>>                  [ Y1,    0,      -Y1,      0],
    >>>                  [ 0,    Y2+G,    -Y2,     -G],
    >>>                  [ -Y1,  -Y2,    Y1+Y2+Y3, -Y3],
    >>>                  [ 0,    -G,      -Y3,      Y3+G ]
    >>>                  ])
    >>> s = sp.symbols('s ', complex=True)
    >>> # Butter de 3er orden doblemente cargado
    >>> Ymai = Ymai.subs(Y1, 1/s/sp.Rational('1'))
    >>> Ymai = Ymai.subs(Y3, 1/s/sp.Rational('1'))
    >>> Ymai = Ymai.subs(Y2, s*sp.Rational('2'))
    >>> # con_detalles = False
    >>> con_detalles = True
    >>> # Calculo la Z en el puerto de entrada a partir de la MAI
    >>> Zmai = calc_MAI_impedance_ij(Ymai, input_port[0], input_port[1], verbose=con_detalles)
    >>> print_latex(a_equal_b_latex_s('Z(s)', Zmai  ))
    Zmai  = (2*G*s + 2*s**2*(G*s + 1) + 1)/(2*G*s**2 + G + 2*s)

    """
    # Check if Ymai is an instance of sp.Matrix
    if not isinstance(Ymai, sp.MatrixBase):
        raise ValueError("Ymai must be an instance of sp.Matrix.")

    # Check if ii and jj are integers
    if not isinstance(ii, int) or not isinstance(jj, int):
        raise ValueError("ii and jj must be integers.")

    # Check if verbose is an instance of bool
    if not isinstance(verbose, bool):
        raise ValueError("verbose must be an instance of bool")

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

#%%
   ##################################
  ## Análisis Nodal Modificado MNA #
 ##################################
#%%

# Modelos usados para los OpAmps. Ver 2.1 Schaumann, R. Design of Analog Filters.
opamp_models_str = ['OA_ideal','OA_1polo', 'OA_integrador']

# parámetros típicos en los opamps de LTspice.
parametros_opamp = ( 'aop', 'gbw', 'aol' )

# posibles nombres para puertos de entrada/salida en el esquemático.
posibles_entradas = ( 'v_v1', 'v_vi', 'v_vin' )
posibles_salidas = ( 'v_v2', 'v_vo', 'v_vout' )

def smna(file_schematic, opamp_model = 'OA_ideal', bAplicarValoresComponentes = True, bAplicarParametros = True):
    """Realiza el análisis nodal modificado (simbólico) al circuito definido en 
    el archivo *file_schematic*. Los formatos aceptados son LTspice 
    y Netlist (ver ejemplos).
    
    El código se basa en el [publicado](https://tiburonboy.github.io/Symbolic-Modified-Nodal-Analysis-using-Python/Introduction.html) por Tony a.k.a @tiburonboy.

    Parameters
    ----------
    file_schematic : nombre de archivo del circuito
        Los formatos de esquemáticos aceptados son de LTspice y Netlist.
    opamp_model : string
        Cómo se van a tratar a los OpAmps que haya en el esquemático. Las opciones son:
            
        * 'OA_ideal': El opamp tiene ganancia y ancho de banda infinito.
        * 'OA_1polo': El opamp tiene una ganancia con un solo polo situado en :math:`:math:`\\frac{G_{BW}}{s+\\frac{G_{BW}}{A_{OL}}}``.
        * 'OA_integrador': El opamp se comporta como un integrador: :math:`\\frac{G_{BW}}{s}`.
        
    bAplicarValoresComponentes : bool, optional
        Se aplicarán los valores de cada componente a la ecuación MNA.
    bAplicarParametros : bool, optional
        Se aplicarán los valores de cada parámerto hallado en el esquemático a la ecuación MNA.

    Returns
    -------
    report: text string
            The net list report.
    df_netlist: pandas dataframe
            circuit net list info loaded into a dataframe
    df_netlist_unknown_currents: pandas dataframe
            branches with unknown currents
    A: sympy matrix
            The A matrix is (m+n) by (m+n) and is the combination of 4 smaller matrices, G, B, C, and D.
            The G matrix is n by n, where n is the number of nodes. The matrix is formed by the interconnections
            between the resistors, capacitors and VCCS type elements. In the original paper G is called Yr, 
            where Yr is a reduced form of the nodal matrix excluding the contributions due to voltage 
            sources, current controlling elements, etc. In python row and columns are: G[row, column]
            The B matrix is an n by m matrix with only 0, 1 and -1 elements, where n = number of nodes
            and m is the number of current unknowns, i_unk. There is one column for each unknown current.
            The code loop through all the branches and process elements that have stamps for the B matrix: 
            The C matrix is an m by n matrix with only 0, 1 and -1 elements (except for controlled sources).
            The code is similar to the B matrix code, except the indices are swapped. The code loops through 
            all the branches and process elements that have stamps for the C matrix: 
            The D matrix is an m by m matrix, where m is the number of unknown currents. 
    X: list
            The X matrix is an (n+m) by 1 vector that holds the unknown quantities (node voltages 
            and the currents through the independent voltage sources). The top n elements are the n node 
            voltages. The bottom m elements represent the currents through the m independent voltage 
            sources in the circuit. The V matrix is n by 1 and holds the unknown voltages. The J matrix 
            is m by 1 and holds the unknown currents through the voltage sources
    Z: list
            The Z matrix holds the independent voltage and current sources and is the combination
            of 2 smaller matrices I and Ev. The Z matrix is (m+n) by 1, n is the number of nodes, 
            and m is the number of independent voltage sources. The I matrix is n by 1 and contains 
            the sum of the currents through the passive elements into the corresponding node (either 
            zero, or the sum of independent current sources). The Ev matrix is m by 1 and holds the 
            values of the independent voltage sources.

    Raises
    ------
    ValueError
        Si el archivo no existe o no es un nombre válido.

    Examples
    --------
    >>> # Para la siguiente red eléctrica:
    >>> # Numeramos los polos de 0 a n=3
    >>> # 
    >>> #     0-------+--Y1----2---Y3--3---
    >>> #                      |           /
    >>> #                     Y2           / R
    >>> #                      |           /
    >>> #     1----------------+-------1----
    >>> #  
    >>> from pytc2.general import print_latex, a_equal_b_latex_s
    >>> from pytc2.cuadripolos import calc_MAI_impedance_ij
    >>> import sympy as sp
    >>> input_port = [0, 1]
    >>> output_port = [3, 1]
    >>> Y1, Y2, Y3 = sp.symbols('Y1 Y2 Y3', complex=True)
    >>> G = sp.symbols('G', real=True, positive=True)
    >>> #      Nodos: 0      1        2        3
    >>> Ymai = sp.Matrix([  
    >>>                  [ Y1,    0,      -Y1,      0],
    >>>                  [ 0,    Y2+G,    -Y2,     -G],
    >>>                  [ -Y1,  -Y2,    Y1+Y2+Y3, -Y3],
    >>>                  [ 0,    -G,      -Y3,      Y3+G ]
    >>>                  ])
    >>> s = sp.symbols('s ', complex=True)
    >>> # Butter de 3er orden doblemente cargado
    >>> Ymai = Ymai.subs(Y1, 1/s/sp.Rational('1'))
    >>> Ymai = Ymai.subs(Y3, 1/s/sp.Rational('1'))
    >>> Ymai = Ymai.subs(Y2, s*sp.Rational('2'))
    >>> # con_detalles = False
    >>> con_detalles = True
    >>> # Calculo la Z en el puerto de entrada a partir de la MAI
    >>> Zmai = calc_MAI_impedance_ij(Ymai, input_port[0], input_port[1], verbose=con_detalles)
    >>> print_latex(a_equal_b_latex_s('Z(s)', Zmai  ))
    Zmai  = (2*G*s + 2*s**2*(G*s + 1) + 1)/(2*G*s**2 + G + 2*s)

    """
    
    if not isinstance(file_schematic, str) or not os.path.exists(file_schematic):
        raise ValueError("file_schematic debe ser el nombre de un archivo que exista.")
    
    if not (isinstance(opamp_model, str) and opamp_model in opamp_models_str):
        raise ValueError(f'El argumento elemento debe ser un string contenido en {opamp_models_str}.')
    
    # Obtener la carpeta (directorio)
    folder_name = os.path.dirname(file_schematic)
    
    # Obtener el nombre del archivo con la extensión
    filename_with_extension = os.path.basename(file_schematic)
    
    # Separar el nombre del archivo de la extensión
    baseFileName, extension = os.path.splitext(filename_with_extension)
    
    fileName_netlist = os.path.join(folder_name, baseFileName + '.net')
    
    if os.path.exists(fileName_netlist) and (os.path.getmtime(file_schematic) < os.path.getmtime(fileName_netlist)):
        print(f'Utilizando netlist: {fileName_netlist}')
    else:
        
        if platform.system() == 'Windows':
            
            print_console_alert('Esta parte no ha sido probada')
            
            ltspice_bin = os.path.expanduser('c:\\Program Files\\LTC\\LTspiceXVII\\XVIIx64.exe')
            subprocess.run([ltspice_bin, '-netlist', file_schematic])
        else:
            home_directory = os.path.expanduser("~")
            ltspice_bin = os.path.expanduser('~/.wine/drive_c/Program Files/LTC/LTspiceXVII/XVIIx64.exe')
            
            # Configurar la variable de entorno WINEPREFIX
            os.environ['WINEPREFIX'] = os.path.join(home_directory, '.wine')    
            print(f'Actualizando netlist a partir de: {file_schematic} ...')
            subprocess.run(['wine', ltspice_bin, '-wine', '-netlist', file_schematic], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            print(f'Netlist generado: {fileName_netlist} ...')
    
    with open(fileName_netlist, 'r', encoding='latin-1') as file:
        # Read the content of the file
        example_net_list = file.read()
    
    def parse_param_line(param_line):
        # Eliminar el prefijo ".param" y los espacios en blanco adicionales
        param_line = param_line.replace('.param', '').strip()
        
        # Usar una expresión regular para encontrar los pares de parámetros y sus definiciones
        pattern = re.compile(r'(\w+)\s*=\s*([^=]+?)(?=\s+\w+\s*=|$)')
        
        # Encontrar todos los pares de parámetros y sus definiciones
        matches = pattern.findall(param_line)
        
        # Convertir los pares en un diccionario
        param_dict = {param.strip(): definition.strip() for param, definition in matches}
        
        return param_dict


    # parámetros que almacenaré del netlist
    dic_params = {}
    
    # initialize variables
    num_rlc = 0 # number of passive elements
    num_ind = 0 # number of inductors
    num_v = 0    # number of independent voltage sources
    num_i = 0    # number of independent current sources
    i_unk = 0  # number of current unknowns
    num_opamps = 0   # number of op amps
    num_vcvs = 0     # number of controlled sources of various types
    num_vccs = 0
    num_cccs = 0
    num_ccvs = 0
    num_cpld_ind = 0 # number of coupled inductors

    content = example_net_list.splitlines()

    content = [x.strip() for x in content]  #remove leading and trailing white space
    # remove empty lines
    while '' in content:
        content.pop(content.index(''))

    # remove comment lines, these start with a asterisk *
    content = [n for n in content if not n.startswith('*')]
    content = [n for n in content if not n.startswith('#')]
    # remove other comment lines, these start with a semicolon ;
    content = [n for n in content if not n.startswith(';')]
    # remove spice directives, these start with a period, .
    #content = [n for n in content if not n.startswith('.')]
    # converts 1st letter to upper case
    #content = [x.upper() for x in content] <- this converts all to upper case
    content = [x.capitalize() for x in content]
    # removes extra spaces between entries
    content = [' '.join(x.split()) for x in content]

    line_cnt = len(content) # number of lines in the netlist
    branch_cnt = 0  # number of branches in the netlist
    # check number of entries on each line, count each element type
    for i in range(line_cnt):
        x = content[i][0]
        tkn = content[i].split()
        tk_cnt = len(tkn) # split the line into a list of words

        if (x == 'R') or (x == 'L') or (x == 'C'):
            # if tk_cnt != 4:
            #     raise Exception("branch {:d} not formatted correctly, {:s} ".format(i,content[i]),
            #     "had {:d} items and should only be 4".format(tk_cnt))
            num_rlc += 1
            branch_cnt += 1
            if x == 'L':
                num_ind += 1
        
        elif x == 'V':
            if tk_cnt < 4:
                raise Exception("branch {:d} not formatted correctly, {:s} ".format(i,content[i]),
                "had {:d} items and should only be 4".format(tk_cnt))
            num_v += 1
            branch_cnt += 1
        elif x == 'I':
            if tk_cnt != 4:
                raise Exception("branch {:d} not formatted correctly, {:s} ".format(i,content[i]),
                "had {:d} items and should only be 4".format(tk_cnt))
            num_i += 1
            branch_cnt += 1
        elif x == 'O':
            if tk_cnt != 4:
                raise Exception("branch {:d} not formatted correctly, {:s} ".format(i,content[i]),
                "had {:d} items and should only be 4".format(tk_cnt))
            num_opamps += 1
        
        #subcircuits
        elif x == 'X':
            
            # opamp subcircuit
            if tkn[4] == 'opamp':
                pattern = r"X(\w+) (\w+) (\w+) (\w+) opamp aol=([\d\w]+) gbw=([\d\w]+)"
                
                # matchear la línea de los opamps ideales
                match = re.match(pattern, content[i])
                
                if not match:
                    raise Exception("branch {:d} not formatted correctly, {:s} ".format(i,content[i]),
                    "Se esperaba una línea estilo: XU1 nodo_v- nodo_v+ nodo_vout opamp Aol=100K GBW=10Meg")
                
                num_opamps += 1
                
        elif x == 'E':
            if (tk_cnt != 6):
                raise Exception("branch {:d} not formatted correctly, {:s} ".format(i,content[i]),
                "had {:d} items and should only be 6".format(tk_cnt))
            num_vcvs += 1
            branch_cnt += 1
        elif x == 'G':
            if (tk_cnt != 6):
                raise Exception("branch {:d} not formatted correctly, {:s} ".format(i,content[i]),
                "had {:d} items and should only be 6".format(tk_cnt))
            num_vccs += 1
            branch_cnt += 1
        elif x == 'F':
            if (tk_cnt != 5):
                raise Exception("branch {:d} not formatted correctly, {:s} ".format(i,content[i]),
                "had {:d} items and should only be 5".format(tk_cnt))
            num_cccs += 1
            branch_cnt += 1
        elif x == 'H':
            if (tk_cnt != 5):
                raise Exception("branch {:d} not formatted correctly, {:s} ".format(i,content[i]),
                "had {:d} items and should only be 5".format(tk_cnt))
            num_ccvs += 1
            branch_cnt += 1
        elif x == 'K':
            if (tk_cnt != 4):
                raise Exception("branch {:d} not formatted correctly, {:s} ".format(i,content[i]),
                "had {:d} items and should only be 4".format(tk_cnt))
            num_cpld_ind += 1

        # parámetros, includes, y otros comandos
        elif x == '.':
            pass
            
        else:
            raise Exception("unknown element type in branch {:d}: {:s}".format(i,content[i]))


    ''' The parser performs the following operations.
     1. puts branch elements into data frame  
     2. counts number of nodes  

     data frame labels:
     - element: type of element  
     - p node: positive node  
     - n node: negative node, for a current source, the arrow point terminal, LTspice 
     puts the inductor phasing dot on this terminal  
     - cp node: controlling positive node of branch  
     - cn node: controlling negative node of branch  
     - Vout: opamp output node  
     - value: value of element or voltage  
     - Vname: voltage source through which the controlling current flows. Need to 
     add a zero volt voltage source to the controlling branch.  
     - Lname1: name of coupled inductor 1  
     - Lname2: name of coupled inductor 2'''  


    # build the pandas data frame
    df_netlist = pd.DataFrame(columns=['element','p node','n node','cp node','cn node',
        'Vout','value','Vname','Lname1','Lname2'])

    # this data frame is for branches with unknown currents
    df_netlist_unknown_currents = pd.DataFrame(columns=['element','p node','n node'])

    # lista de nombre de los nodos
    node_names = ['0'] # nodo 0 es GND en LTspice
    cant_nodos = 1

    # node_names = [] # nodo 0 es GND en LTspice
    # cant_nodos = 0
    

    # ### Functions to load branch elements into data frame and check for gaps in node numbering


        
    def parse_value_or_expression(input_string):
        # Patrón para expresiones matemáticas entre llaves
        expression_pattern = r'^\{(.+)\}$'
    
        # Intentar parsear como expresión matemática
        match = re.match(expression_pattern, input_string)
        if match:
            expr_str = match.group(1)
            try:
                sp_expr = sp.sympify(expr_str)
                return(sp_expr) 
            except sp.sympifyError:
                raise ValueError(f"No se pudo parsear la expresión matemática: {expr_str}")
    
        # Si no es una expresión matemática, intentar parsear como número
        number_pattern = r'^-?\d+(\.\d+)?$'
        if re.match(number_pattern, input_string):
            return float(input_string) if '.' in input_string else int(input_string)
    
        raise ValueError(f"El string: {input_string} no es un número válido ni una expresión matemática entre llaves.")
        

    def parse_nodes(str_node):
        
        nonlocal cant_nodos, node_names
        
        pattern = r'^n(\d{3})$'
        match = re.match(pattern, str_node)
    
        if match:
            # valor numérico consecutivo que asigna LTspice a los nodos.
            
            bAux = np.asarray(node_names) == str_node
            if np.any(bAux):
                this_node_int = np.flatnonzero(bAux)[0]
            else:
                # agrego nodo a la lista    
                this_node_int = cant_nodos
                node_names += [str_node]
                cant_nodos += 1
        else:
            # etiqueta de nodo
            bAux = np.asarray(node_names) == str_node
            if np.any(bAux):
                this_node_int = np.flatnonzero(bAux)[0]
            else:
                # agrego nodo a la lista    
                this_node_int = cant_nodos
                node_names += [str_node]
                cant_nodos += 1
           
        #     raise ValueError("El string: {:s} no tiene el formato 'nXXX' donde XXX es un número entero con tres dígitos.".format(str_node))
    
        return this_node_int

    # Diccionario de conversión de sufijos a escala
    suffix_scale = {
        'T': 1e12, 't': 1e12,
        'G': 1e9, 'g': 1e9,
        'MEG': 1e6, 'meg': 1e6, 'Meg': 1e6,
        'K': 1e3, 'k': 1e3,
        'M': 1e-3, 'm': 1e-3,
        'U': 1e-6, 'u': 1e-6, 'µ': 1e-6,
        'N': 1e-9, 'n': 1e-9,
        'P': 1e-12, 'p': 1e-12,
        'F': 1e-15, 'f': 1e-15
    }

    def parse_opamp(opamp_tokens):
    
        value_str = opamp_tokens[5]
        pattern = re.compile(r'aol=(\d+(\.\d+)?)([TtGgMmUuNnPpFf]|[Kk]|MEG|meg|Meg)?')
        match = pattern.match(value_str)
        if not match:
            raise ValueError(f"El valor '{value_str}' no tiene un formato válido")
        
        value = float(match.group(1))
        suffix = match.group(3)
       
        if suffix:
            aol = value * suffix_scale[suffix]
        else:
            aol = np.nan
            
        value_str = opamp_tokens[6]
        pattern = re.compile(r'gbw=(\d+(\.\d+)?)([TtGgMmUuNnPpFf]|[Kk]|MEG|meg|Meg)?')
        match = pattern.match(value_str)
        if not match:
            raise ValueError(f"El valor '{value_str}' no tiene un formato válido")
        
        value = float(match.group(1))
        suffix = match.group(3)
       
        if suffix:
            gbw = value * suffix_scale[suffix]
        else:
            gbw = np.nan
        
        return [aol, gbw]

    # loads voltage or current sources into branch structure
    def indep_source(line_df_netlist, line_nu):
        tk = content[line_nu].split()
        df_netlist.loc[line_df_netlist,'element'] = tk[0]
        df_netlist.loc[line_df_netlist,'p node'] = parse_nodes(tk[1])
        df_netlist.loc[line_df_netlist,'n node'] = parse_nodes(tk[2])
        if tk[3] == 'ac':
            df_netlist.loc[line_df_netlist,'value'] = parse_value_or_expression(tk[4])
        else:
            df_netlist.loc[line_df_netlist,'value'] = parse_value_or_expression(tk[3])

    # loads passive elements into branch structure
    def rlc_element(line_df_netlist, line_nu):
        tk = content[line_nu].split()
        df_netlist.loc[line_df_netlist,'element'] = tk[0]
        df_netlist.loc[line_df_netlist,'p node'] = parse_nodes(tk[1])
        df_netlist.loc[line_df_netlist,'n node'] = parse_nodes(tk[2])
        df_netlist.loc[line_df_netlist,'value'] = parse_value_or_expression(tk[3])

    # loads multi-terminal elements into branch structure
    # O - Op Amps
    def opamp_sub_network(line_df_netlist, line_nu):
        tk = content[line_nu].split()
        
        # forzamos tipo de elemento O
        this_element = tk[0]
        if this_element[0] == 'O':
            
            bOpampInvertInput = False
        else:
            bOpampInvertInput = True
            this_element = 'O' + this_element

        df_netlist.loc[line_df_netlist,'element'] = this_element
        
        if bOpampInvertInput:
            df_netlist.loc[line_df_netlist,'p node'] = parse_nodes(tk[2])
            df_netlist.loc[line_df_netlist,'n node'] = parse_nodes(tk[1])
        else:
            df_netlist.loc[line_df_netlist,'p node'] = parse_nodes(tk[1])
            df_netlist.loc[line_df_netlist,'n node'] = parse_nodes(tk[2])
                
        df_netlist.loc[line_df_netlist,'Vout'] = parse_nodes(tk[3])

        if( len(tk) > 4 and tk[4] == 'opamp' ):
            # usamos el modelo que impone LTspice o el ideal
            aol, gbw = parse_opamp(tk)
        else:
            # usamos valores idealizadamente altos
            aol, gbw = [10e10 , 10e20]
        
        if opamp_model == 'OA_integrador':
            # integrator opamp model
            df_netlist.loc[line_df_netlist,'value'] = gbw/s
            
        elif opamp_model == 'OA_1polo':
            # one-pole opamp model
            df_netlist.loc[line_df_netlist,'value'] = gbw/(s+gbw/aol)
            
        else:
            # modelo ideal de opamp
            aop = sp.symbols('aop', Real = True)
            df_netlist.loc[line_df_netlist,'value'] = aop
                

    # G - VCCS
    def vccs_sub_network(line_df_netlist, line_nu):
        tk = content[line_nu].split()
        df_netlist.loc[line_df_netlist,'element'] = tk[0]
        df_netlist.loc[line_df_netlist,'p node'] = parse_nodes(tk[1])
        df_netlist.loc[line_df_netlist,'n node'] = parse_nodes(tk[2])
        df_netlist.loc[line_df_netlist,'cp node'] = parse_nodes(tk[3])
        df_netlist.loc[line_df_netlist,'cn node'] = parse_nodes(tk[4])
        df_netlist.loc[line_df_netlist,'value'] = parse_value_or_expression(tk[5])

    # E - VCVS
    # in sympy E is the number 2.718, replacing E with Ea otherwise, sp.sympify() errors out
    def vcvs_sub_network(line_df_netlist, line_nu):
        tk = content[line_nu].split()
        df_netlist.loc[line_df_netlist,'element'] = tk[0].replace('E', 'Ea')
        df_netlist.loc[line_df_netlist,'p node'] = parse_nodes(tk[1])
        df_netlist.loc[line_df_netlist,'n node'] = parse_nodes(tk[2])
        df_netlist.loc[line_df_netlist,'cp node'] = parse_nodes(tk[3])
        df_netlist.loc[line_df_netlist,'cn node'] = parse_nodes(tk[4])
        df_netlist.loc[line_df_netlist,'value'] = parse_value_or_expression(tk[5])

    # F - CCCS
    def cccs_sub_network(line_df_netlist, line_nu):
        tk = content[line_nu].split()
        df_netlist.loc[line_df_netlist,'element'] = tk[0]
        df_netlist.loc[line_df_netlist,'p node'] = parse_nodes(tk[1])
        df_netlist.loc[line_df_netlist,'n node'] = parse_nodes(tk[2])
        df_netlist.loc[line_df_netlist,'Vname'] = tk[3].capitalize()
        df_netlist.loc[line_df_netlist,'value'] = parse_value_or_expression(tk[4])

    # H - CCVS
    def ccvs_sub_network(line_df_netlist, line_nu):
        tk = content[line_nu].split()
        df_netlist.loc[line_df_netlist,'element'] = tk[0]
        df_netlist.loc[line_df_netlist,'p node'] = parse_nodes(tk[1])
        df_netlist.loc[line_df_netlist,'n node'] = parse_nodes(tk[2])
        df_netlist.loc[line_df_netlist,'Vname'] = tk[3].capitalize()
        df_netlist.loc[line_df_netlist,'value'] = parse_value_or_expression(tk[4])

    # K - Coupled inductors
    def cpld_ind_sub_network(line_df_netlist, line_nu):
        tk = content[line_nu].split()
        df_netlist.loc[line_df_netlist,'element'] = tk[0]
        df_netlist.loc[line_df_netlist,'Lname1'] = tk[1].capitalize()
        df_netlist.loc[line_df_netlist,'Lname2'] = tk[2].capitalize()
        df_netlist.loc[line_df_netlist,'value'] = parse_value_or_expression(tk[3])

    # function to scan df_netlist and get largest node number
    def count_nodes():
        # need to check that nodes are consecutive
        # fill array with node numbers
        p = np.zeros(cant_netlist_valid)
        for i in range(cant_netlist_valid):
            # need to skip coupled inductor 'K' statements
            if df_netlist.loc[i,'element'][0] != 'K': #get 1st letter of element name
                p[df_netlist['p node'][i]] = df_netlist['p node'][i]
                p[df_netlist['n node'][i]] = df_netlist['n node'][i]

        # find the largest node number
        if df_netlist['n node'].max() > df_netlist['p node'].max():
            largest = df_netlist['n node'].max()
        else:
            largest =  df_netlist['p node'].max()

        largest = int(largest)
        # check for unfilled elements, skip node 0
        for i in range(1,largest):
            if p[i] == 0:
                raise Exception('nodes not in continuous order, node {:.0f} is missing'.format(p[i-1]+1))
        return largest

    # load branch info into data frame
    jj = 0
    
    for i in range(line_cnt):
        x = content[i][0]
        
        tkn = content[i].split()
        tk_cnt = len(tkn) # split the line into a list of words

        if (x == 'R') or (x == 'L') or (x == 'C'):
            rlc_element(jj, i)
            jj += 1
            
        elif x == '.':
            # parámetros
            if tkn[0] == '.param':
                dic_this_pars = parse_param_line(content[i])
                dic_params.update(dic_this_pars)
                
        elif (x == 'V') or (x == 'I'):
            indep_source(jj, i)
            jj += 1
        elif x == 'O' or (x == 'X' and tkn[4] == 'opamp'):
            opamp_sub_network(jj, i)
            jj += 1
        elif x == 'E':
            vcvs_sub_network(jj, i)
            jj += 1
        elif x == 'G':
            vccs_sub_network(jj, i)
            jj += 1
        elif x == 'F':
            cccs_sub_network(jj, i)
            jj += 1
        elif x == 'H':
            ccvs_sub_network(jj, i)
            jj += 1
        elif x == 'K':
            cpld_ind_sub_network(jj, i)
            jj += 1
        else:
            raise Exception("unknown element type in branch {:d}, {:s}".format(i,content[i]))

    '''29 Nov 2023:  When the D matrix is built, independent voltage sources are processed
    in the data frame order when building the D matrix. If the voltage source followed element
    L, H, F, K types in the netlist, a row was inserted that put the voltage source in a different
    row in relation to its position in the Ev matrix. This would cause the node attached to 
    the terminal of the voltage source to be zero volts.  
    Solution - The following block of code was added to move voltage source types to the 
    beginning of the net list dataframe before any calculations are performed.''' 

    cant_netlist_valid = len(df_netlist)

    # Check for position of voltages sources in the dataframe.
    source_index = [] # keep track of voltage source row number
    other_index = [] # make a list of all other types
    for i in range(cant_netlist_valid):
        # process all the elements creating unknown currents
        x = df_netlist.loc[i,'element'][0]   #get 1st letter of element name
        if (x == 'V'):
            source_index.append(i)
        else:
            other_index.append(i)

    df_netlist = df_netlist.reindex(source_index+other_index,copy=True) # re-order the data frame
    df_netlist.reset_index(drop=True, inplace=True) # renumber the index

    # count number of nodes
    num_nodes = count_nodes()

    # Build df_netlist_unknown_currents: consists of branches with current unknowns, used for C & D matrices
    # walk through data frame and find these parameters
    count = 0
    for i in range(cant_netlist_valid):
        # process all the elements creating unknown currents
        x = df_netlist.loc[i,'element'][0]   #get 1st letter of element name
        if (x == 'L') or (x == 'V') or (x == 'O') or (x == 'E') or (x == 'H') or (x == 'F'):
            df_netlist_unknown_currents.loc[count,'element'] = df_netlist.loc[i,'element']
            df_netlist_unknown_currents.loc[count,'p node'] = df_netlist.loc[i,'p node']
            df_netlist_unknown_currents.loc[count,'n node'] = df_netlist.loc[i,'n node']
            count += 1

    # print the netlist report
    report = 'Net list report\n'
    report = report+('number of lines in netlist: {:d}\n'.format(line_cnt))
    report = report+'number of branches: {:d}\n'.format(branch_cnt)
    report = report+'number of nodes: {:d}\n'.format(num_nodes)
    # count the number of element types that affect the size of the B, C, D, E and J arrays
    # these are current unknows
    i_unk = num_v+num_opamps+num_vcvs+num_ccvs+num_cccs+num_ind
    report = report+'number of unknown currents: {:d}\n'.format(i_unk)
    report = report+'number of RLC (passive components): {:d}\n'.format(num_rlc)
    report = report+'number of inductors: {:d}\n'.format(num_ind)
    report = report+'number of independent voltage sources: {:d}\n'.format(num_v)
    report = report+'number of independent current sources: {:d}\n'.format(num_i)
    report = report+'number of op amps: {:d}\n'.format(num_opamps)
    report = report+'number of E - VCVS: {:d}\n'.format(num_vcvs)
    report = report+'number of G - VCCS: {:d}\n'.format(num_vccs)
    report = report+'number of F - CCCS: {:d}\n'.format(num_cccs)
    report = report+'number of H - CCVS: {:d}\n'.format(num_ccvs)
    report = report+'number of K - Coupled inductors: {:d}\n'.format(num_cpld_ind)

    # initialize some sp.Symbolic matrix with sp.zeros
    # A is formed by [[G, C] [B, D]]
    # Z = [I,E]
    # X = [V, J]
    V = sp.zeros(num_nodes,1)
    I = sp.zeros(num_nodes,1)
    G = sp.zeros(num_nodes,num_nodes)  # also called Yr, the reduced nodal matrix

    # count the number of element types that affect the size of the B, C, D, E and J arrays
    # these are element types that have unknown currents
    i_unk = num_v+num_opamps+num_vcvs+num_ccvs+num_ind+num_cccs
    # if i_unk == 0, just generate empty arrays
    B = sp.zeros(num_nodes,i_unk)
    C = sp.zeros(i_unk,num_nodes)
    D = sp.zeros(i_unk,i_unk)
    Ev = sp.zeros(i_unk,1)
    J = sp.zeros(i_unk,1)

    ''' The G matrix is n by n, where n is the number of nodes. 
    The matrix is formed by the interconnections between the resistors, 
    capacitors and VCCS type elements.  In the original paper G is called Yr, 
    where Yr is a reduced form of the nodal matrix excluding the contributions 
    due to voltage sources, current controlling elements, etc.  In python row 
    and columns are: G[row, column]'''
    for i in range(cant_netlist_valid):  # process each row in the data frame
        n1 = df_netlist.loc[i,'p node']
        n2 = df_netlist.loc[i,'n node']
        cn1 = df_netlist.loc[i,'cp node']
        cn2 = df_netlist.loc[i,'cn node']
        # process all the passive elements, save conductance to temp value
        x = df_netlist.loc[i,'element'][0]   #get 1st letter of element name
        if x == 'R':
            g = 1/sp.sympify(df_netlist.loc[i,'element'])
        if x == 'C':
            g = s*sp.sympify(df_netlist.loc[i,'element'])
        if x == 'G':   #vccs type element
            g = sp.sympify(df_netlist.loc[i,'element'].lower())  # use a sp.Symbol for gain value

        if (x == 'R') or (x == 'C'):
            # If neither side of the element is connected to ground
            # then subtract it from the appropriate location in the matrix.
            if (n1 != 0) and (n2 != 0):
                G[n1-1,n2-1] += -g
                G[n2-1,n1-1] += -g

            # If node 1 is connected to ground, add element to diagonal of matrix
            if n1 != 0:
                G[n1-1,n1-1] += g

            # same for for node 2
            if n2 != 0:
                G[n2-1,n2-1] += g

        if x == 'G':    #vccs type element
            # check to see if any terminal is grounded
            # then stamp the matrix
            if n1 != 0 and cn1 != 0:
                G[n1-1,cn1-1] += g

            if n2 != 0 and cn2 != 0:
                G[n2-1,cn2-1] += g

            if n1 != 0 and cn2 != 0:
                G[n1-1,cn2-1] -= g

            if n2 != 0 and cn1 != 0:
                G[n2-1,cn1-1] -= g

    '''The B matrix is an n by m matrix with only 0, 1 and -1 elements, where 
    n = number of nodes and m is the number of current unknowns, i_unk. There is 
    one column for each unknown current. The code loop through all the branches 
    and process elements that have stamps for the B matrix:  
     - Voltage sources (V)  
     - Opamps (O)  
     - CCVS (H)  
     - CCCS (F)  
     - VCVS (E)  
     - Inductors (L)    

    The order of the columns is as they appear in the netlist.  CCCS (F) does not get
    its own column because the controlling current is through a zero volt voltage source,
    called Vname and is already in the net list.'''
    sn = 0   # count source number as code walks through the data frame
    for i in range(cant_netlist_valid):
        n1 = df_netlist.loc[i,'p node']
        n2 = df_netlist.loc[i,'n node']
        n_vout = df_netlist.loc[i,'Vout'] # node connected to op amp output

        # process elements with input to B matrix
        x = df_netlist.loc[i,'element'][0]   #get 1st letter of element name
        if x == 'V':
            if i_unk > 1:  #is B greater than 1 by n?, V
                if n1 != 0:
                    B[n1-1,sn] = 1
                if n2 != 0:
                    B[n2-1,sn] = -1
            else:
                if n1 != 0:
                    B[n1-1] = 1
                if n2 != 0:
                    B[n2-1] = -1
            sn += 1   #increment source count
        if x == 'O':  # op amp type, output connection of the opamp goes in the B matrix
            B[n_vout-1,sn] = 1
            sn += 1   # increment source count
        if (x == 'H') or (x == 'F'):  # H: ccvs, F: cccs,
            if i_unk > 1:  #is B greater than 1 by n?, H, F
                # check to see if any terminal is grounded
                # then stamp the matrix
                if n1 != 0:
                    B[n1-1,sn] = 1
                if n2 != 0:
                    B[n2-1,sn] = -1
            else:
                if n1 != 0:
                    B[n1-1] = 1
                if n2 != 0:
                    B[n2-1] = -1
            sn += 1   #increment source count
        if x == 'E':   # vcvs type, only ik column is altered at n1 and n2
            if i_unk > 1:  #is B greater than 1 by n?, E
                if n1 != 0:
                    B[n1-1,sn] = 1
                if n2 != 0:
                    B[n2-1,sn] = -1
            else:
                if n1 != 0:
                    B[n1-1] = 1
                if n2 != 0:
                    B[n2-1] = -1
            sn += 1   #increment source count
        if x == 'L':
            if i_unk > 1:  #is B greater than 1 by n?, L
                if n1 != 0:
                    B[n1-1,sn] = 1
                if n2 != 0:
                    B[n2-1,sn] = -1
            else:
                if n1 != 0:
                    B[n1-1] = 1
                if n2 != 0:
                    B[n2-1] = -1
            sn += 1   #increment source count

    # check source count
    if sn != i_unk:
        raise Exception('source number, sn={:d} not equal to i_unk={:d} in matrix B'.format(sn,i_unk))

    ''' The C matrix is an m by n matrix with only 0, 1 and -1 elements (except for controlled sources).  
    The code is similar to the B matrix code, except the indices are swapped.   The code loops through 
    all the branches and process elements that have stamps for the C matrix:  
     - Voltage sources (V)  
     - Opamps (O)  
     - CCVS (H)  
     - CCCS (F)  
     - VCVS (E)  
     - Inductors (L)  

     Op Amp elements
     The op amp element is assumed to be an ideal op amp and use of this component is valid only when 
     used in circuits with a DC path (a short or a resistor) from the output terminal to the negative 
     input terminal of the op amp. No error checking is provided and if the condition is violated, 
     the results likely will be erroneous.   

     References use in the debugging of the opamp stamp:   
     1. Design of Analog Circuits Through sp.Symbolic Analysis, edited by Mourad Fakhfakh, Esteban Tlelo-Cuautle, Francisco V. Fernández   
     2. Computer Aided Design and Design Automation, edited by Wai-Kai Chen  

     find the the column position in the C and D matrix for controlled sources
     needs to return the node numbers and branch number of controlling branch'''
    def find_vname(name):
        # need to walk through data frame and find these parameters
        for i in range(len(df_netlist_unknown_currents)):
            # process all the elements creating unknown currents
            if name == df_netlist_unknown_currents.loc[i,'element']:
                n1 = df_netlist_unknown_currents.loc[i,'p node']
                n2 = df_netlist_unknown_currents.loc[i,'n node']
                return n1, n2, i  # n1, n2 & col_num are from the branch of the controlling element

        raise Exception('failed to find matching branch element in find_vname')

    # generate the C Matrix
    sn = 0   # count source number as code walks through the data frame
    for i in range(cant_netlist_valid):
        n1 = df_netlist.loc[i,'p node']
        n2 = df_netlist.loc[i,'n node']
        cn1 = df_netlist.loc[i,'cp node'] # nodes for controlled sources
        cn2 = df_netlist.loc[i,'cn node']
        n_vout = df_netlist.loc[i,'Vout'] # node connected to op amp output

        # process elements with input to B matrix
        x = df_netlist.loc[i,'element'][0]   #get 1st letter of element name
        if x == 'V':
            if i_unk > 1:  #is B greater than 1 by n?, V
                if n1 != 0:
                    C[sn,n1-1] = 1
                if n2 != 0:
                    C[sn,n2-1] = -1
            else:
                if n1 != 0:
                    C[n1-1] = 1
                if n2 != 0:
                    C[n2-1] = -1
            sn += 1   #increment source count

        if x == 'O':  # op amp type, input connections of the opamp go into the C matrix

            #is B greater than 1 by n?, O            
            if i_unk > 1:  
                # check to see if any terminal is grounded
                # then stamp the matrix
                if n1 != 0:
                    C[sn,n1-1] = 1
                if n2 != 0:
                    C[sn,n2-1] = -1

                # if opamp_model != 'OA_ideal':
                # la salida del opamp, tendrá un valor de A(s)(V+ - V-)
                # el valor de A(s) indica el modelo que usamos del opamp.
                # la matriz C modela A(s) como B = -1/A(s) según: 
                # (Vlach1994) Vlach, Jiri - Linear Circuit Theory_ Matrices in 
                # Computer Applications-Apple Academic Press.
                C[sn,n_vout-1] = -1/df_netlist.loc[i,'value']
                    
            else:
                if n1 != 0:
                    C[n1-1] = 1
                if n2 != 0:
                    C[n2-1] = -1

                # if opamp_model != 'OA_ideal':
                C[n_vout-1] = -1/df_netlist.loc[i,'value']
            
            sn += 1   # increment source count

        if x == 'F':  # need to count F (cccs) types
            sn += 1   #increment source count
        if x == 'H':  # H: ccvs
            if i_unk > 1:  #is B greater than 1 by n?, H
                # check to see if any terminal is grounded
                # then stamp the matrix
                if n1 != 0:
                    C[sn,n1-1] = 1
                if n2 != 0:
                    C[sn,n2-1] = -1
            else:
                if n1 != 0:
                    C[n1-1] = 1
                if n2 != 0:
                    C[n2-1] = -1
            sn += 1   #increment source count
        if x == 'E':   # vcvs type, ik column is altered at n1 and n2, cn1 & cn2 get value
            if i_unk > 1:  #is B greater than 1 by n?, E
                if n1 != 0:
                    C[sn,n1-1] = 1
                if n2 != 0:
                    C[sn,n2-1] = -1
                # add entry for cp and cn of the controlling voltage
                if cn1 != 0:
                    C[sn,cn1-1] = -sp.sympify(df_netlist.loc[i,'element'].lower())
                if cn2 != 0:
                    C[sn,cn2-1] = sp.sympify(df_netlist.loc[i,'element'].lower())
            else:
                if n1 != 0:
                    C[n1-1] = 1
                if n2 != 0:
                    C[n2-1] = -1
                vn1, vn2, df_netlist_unknown_currents_index = find_vname(df_netlist.loc[i,'Vname'])
                if vn1 != 0:
                    C[vn1-1] = -sp.sympify(df_netlist.loc[i,'element'].lower())
                if vn2 != 0:
                    C[vn2-1] = sp.sympify(df_netlist.loc[i,'element'].lower())
            sn += 1   #increment source count

        if x == 'L':
            if i_unk > 1:  #is B greater than 1 by n?, L
                if n1 != 0:
                    C[sn,n1-1] = 1
                if n2 != 0:
                    C[sn,n2-1] = -1
            else:
                if n1 != 0:
                    C[n1-1] = 1
                if n2 != 0:
                    C[n2-1] = -1
            sn += 1   #increment source count

    # check source count
    if sn != i_unk:
        raise Exception('source number, sn={:d} not equal to i_unk={:d} in matrix C'.format(sn,i_unk))

    ''' The D matrix is an m by m matrix, where m is the number of unknown currents.  
    m = i_unk = num_v+num_opamps+num_vcvs+num_ccvs+num_ind+num_cccs

    Stamps that affect the D matrix are: inductor, ccvs and cccs  
    inductors: minus sign added to keep current flow convention consistent  

    Coupled inductors notes:  
    Can the K statement be anywhere in the net list, even before Lx and Ly?   
    12/6/2017 doing some debugging on with coupled inductors  
    LTspice seems to put the phasing dot on the neg node when it generates the netlist   
    This code uses M for mutual inductance, LTspice uses k for the coupling coefficient.
    Nota: Lo cambiamos ya que M = k.sqrt(L1.L2)
    '''  
    

    # generate the D Matrix
    sn = 0   # count source number as code walks through the data frame
    for i in range(cant_netlist_valid):
        n1 = df_netlist.loc[i,'p node']
        n2 = df_netlist.loc[i,'n node']
        #cn1 = df_netlist.loc[i,'cp node'] # nodes for controlled sources
        #cn2 = df_netlist.loc[i,'cn node']
        #n_vout = df_netlist.loc[i,'Vout'] # node connected to op amp output

        # process elements with input to D matrix
        x = df_netlist.loc[i,'element'][0]   #get 1st letter of element name
        if (x == 'V') or (x == 'O') or (x == 'E'):  # need to count V, E & O types
            sn += 1   #increment source count

        if x == 'L':
            if i_unk > 1:  #is D greater than 1 by 1?
                D[sn,sn] += -s*sp.sympify(df_netlist.loc[i,'element'])
            else:
                D[sn] += -s*sp.sympify(df_netlist.loc[i,'element'])
            sn += 1   #increment source count

        if x == 'H':  # H: ccvs
            # if there is a H type, D is m by m
            # need to find the vn for Vname
            # then stamp the matrix
            vn1, vn2, df_netlist_unknown_currents_index = find_vname(df_netlist.loc[i,'Vname'])
            D[sn,df_netlist_unknown_currents_index] += -sp.sympify(df_netlist.loc[i,'element'].lower())
            sn += 1   #increment source count

        if x == 'F':  # F: cccs
            # if there is a F type, D is m by m
            # need to find the vn for Vname
            # then stamp the matrix
            vn1, vn2, df_netlist_unknown_currents_index = find_vname(df_netlist.loc[i,'Vname'])
            D[sn,df_netlist_unknown_currents_index] += -sp.sympify(df_netlist.loc[i,'element'].lower())
            D[sn,sn] = 1
            sn += 1   #increment source count

        if x == 'K':  # K: coupled inductors, KXX LYY LZZ value
            # if there is a K type, D is m by m
            vn1, vn2, ind1_index = find_vname(df_netlist.loc[i,'Lname1'])  # get i_unk position for Lx
            vn1, vn2, ind2_index = find_vname(df_netlist.loc[i,'Lname2'])  # get i_unk position for Ly
            # enter sM on diagonals = value*sqrt(LXX*LZZ)

            this_M = 'M{:s}'.format(df_netlist.loc[i,'element'].lower()[1:])
            D[ind1_index,ind2_index] += -s*sp.sympify(this_M)  # s*Mxx
            D[ind2_index,ind1_index] += -s*sp.sympify(this_M)  # -s*Mxx
            
            # el valor del trafo para a ser M = k*sqrt(L1.L2)
            if bAplicarValoresComponentes:
                inductor_rows = df_netlist[df_netlist['element'].isin([df_netlist.loc[i,'Lname1'], df_netlist.loc[i,'Lname2']])]            
                dic_params.update({this_M : df_netlist.loc[i,'value'] * sp.sqrt(inductor_rows['value'].values[0]*inductor_rows['value'].values[1])})
            else:
                dic_params.update({this_M : df_netlist.loc[i,'value'] * sp.sympify('sqrt({:s}*{:s})'.format(df_netlist.loc[i,'Lname1'], df_netlist.loc[i,'Lname2'] ))} )

            # D[ind1_index,ind2_index] += -s*sp.sympify('{:s}*sqrt({:s}*{:s})'.format(df_netlist.loc[i,'element'], df_netlist.loc[i,'Lname1'], df_netlist.loc[i,'Lname2'] ))  # -s*k*sqrt(L1.L2)
            # D[ind2_index,ind1_index] += -s*sp.sympify('{:s}*sqrt({:s}*{:s})'.format(df_netlist.loc[i,'element'], df_netlist.loc[i,'Lname1'], df_netlist.loc[i,'Lname2'] ))  # -s*k*sqrt(L1.L2)

    ''' The V matrix is an n by 1 matrix formed of the node voltages, where n is the number of nodes. Each element in V corresponds to the voltage at the node.  
    Maybe make small v's v_1 so as not to confuse v1 with V1.'''
    # generate the V matrix
    for i in range(num_nodes):
        # V[i] = sp.sympify('v{:d}'.format(i+1))
        V[i] = sp.Symbol( 'v_{:s}'.format(node_names[i+1]) )

    ''' The J matrix is an m by 1 matrix, where m is the number of unknown currents.
    i_unk = num_v+num_opamps+num_vcvs+num_ccvs+num_ind+num_cccs
    The J matrix is an m by 1 matrix, with one entry for each i_unk from a source'''
    for i in range(len(df_netlist_unknown_currents)):
        # process all the unknown currents
        J[i] = sp.sympify('I_{:s}'.format(df_netlist_unknown_currents.loc[i,'element']))

    ''' The I matrix is an n by 1 matrix, where n is the number of nodes. The value
    of each element of I is determined by the sum of current sources into the 
    corresponding node. If there are no current sources connected to the node, the value is zero.'''

    # generate the I matrix, current sources have n2 = arrow end of the element
    sni = 0   # count current source 
    for i in range(cant_netlist_valid):
        n1 = df_netlist.loc[i,'p node']
        n2 = df_netlist.loc[i,'n node']
        # process all the passive elements, save conductance to temp value
        x = df_netlist.loc[i,'element'][0]   #get 1st letter of element name
        if x == 'I':
            sni += 1   # count current source 
            g = sp.sympify(df_netlist.loc[i,'element'])
            # sum the current into each node
            if n1 != 0:
                I[n1-1] -= g
            if n2 != 0:
                I[n2-1] += g

    # The Ev matrix is m by 1 and holds the values of the independent voltage sources.
    sn = 0   # count source number
    for i in range(cant_netlist_valid):
        # process all the passive elements
        x = df_netlist.loc[i,'element'][0]   #get 1st letter of element name
        if x == 'V':
            Ev[sn] = sp.sympify(df_netlist.loc[i,'element'])
            sn += 1

    if (sn + sni) == 0:
        print_console_alert('No se han encontrado generadores en el esquem')
        
        

    ''' The Z matrix holds the independent voltage and current sources and is the combination of 2
    smaller matrices I and Ev. The Z matrix is (m+n) by 1, n is the number of nodes, and m is the
    number of independent voltage sources. The I matrix is n by 1 and contains the sum of the currents
    through the passive elements into the corresponding node (either zero, or the sum of independent
    current sources). The Ev matrix is m by 1 and holds the values of the independent voltage sources.'''
    Z = I[:] + Ev[:]  # the + operator in python concatenates the lists

    # Put matricies into SymPy 
    Z = sp.Matrix(Z)

    ''' The X matrix is an (n+m) by 1 vector that holds the unknown quantities (node voltages and the currents through
    the independent voltage sources). The top n elements are the n node voltages. The bottom m elements represent the
    currents through the m independent voltage sources in the circuit. The V matrix is n by 1 and holds the unknown voltages.
    The J matrix is m by 1 and holds the unknown currents through the voltage sources '''
    X = V[:] + J[:]  # the + operator in python concatenates the lists

    # Put matricies into SymPy 
    X = sp.Matrix(X)

    # The A matrix is (m+n) by (m+n) and will be developed as the combination of 4 smaller matrices, G, B, C, and D.
    n = num_nodes
    m = i_unk
    A = sp.zeros(m+n,m+n)
    for i in range(n):
        for j in range(n):
            A[i,j] = G[i,j]

    if i_unk > 1:
        for i in range(n):
            for j in range(m):
                A[i,n+j] = B[i,j]
                A[n+j,i] = C[j,i]

        for i in range(m):
            for j in range(m):
                A[n+i,n+j] = D[i,j]

    if i_unk == 1:
        for i in range(n):
            A[i,n] = B[i]
            A[n,i] = C[i]
        A[n,n] = D[0] # added 1/7/2024 while debugging source free circuit with one inductor

    # recolecto los símbolos 
    dic_comp_name_vals = dict(zip(df_netlist['element'][1:], df_netlist['value'][1:]))
    
    def convert_index_to_name(x): 
        if pd.isna(x):
            return x
        return node_names[int(x)]
    
    def translate_node_names(df_in):
    
        df_out = df_in.copy()
        df_out['p node'] = df_in['p node'].apply(convert_index_to_name)
        df_out['n node'] = df_in['n node'].apply(convert_index_to_name)
        df_out['Vout'] = df_in['Vout'].apply(convert_index_to_name)
        
        # print(df_out)
        
        return df_out
    
    df_netlist = translate_node_names(df_netlist)


    node_sym = [ ii for ii in X.free_symbols ]
    node_sym_names = [ str(ii) for ii in X.free_symbols ]
    _, v_in_idx, _ = np.intersect1d(node_sym_names, posibles_entradas, return_indices=True)
    _, v_out_idx, _ = np.intersect1d(node_sym_names, posibles_salidas, return_indices=True)
    
    if len(v_in_idx) > 0:
        v_in = node_sym[v_in_idx[0]]
    else:
        v_in = 0
        
    if len(v_out_idx) > 0:
        v_out = node_sym[v_out_idx[0]]
    else:
        v_out = 0

    # aplicar los valores de los componentes.
    if bAplicarValoresComponentes:
        A = A.subs(dic_comp_name_vals)
    
    # aplicar parametrizaciones del esquemático.
    if bAplicarParametros:
        A = A.subs(dic_params)
    
    mna_sym = [ ii for ii in A.free_symbols ]
    mna_sym_names = [ str(ii) for ii in A.free_symbols ]
    
    _, opamp_idx, _ = np.intersect1d( mna_sym_names, parametros_opamp, return_indices=True)
    
    # Los parámetros del opamp que figuren en las ecuaciones se resguardan. 
    if len(opamp_idx) > 0:
        aop = mna_sym[opamp_idx[0]]
    else:
        aop = 0
    
    # eps es un valor cercano a 0 útil para la simulación circuital. 
    # Externamente puede anularse para facilitar el análisis simbólico
    _, this_idx, _ = np.intersect1d( mna_sym_names, ('eps', ), return_indices=True)
    
    if len(this_idx) > 0:
        eps = mna_sym[this_idx[0]]
    else:
        eps = 0

    equ_smna = sp.Eq(A*X,Z)
    
    extra_results = {
        "A": A,
        "X": X,
        "Z": Z,
        "v_in": v_in,
        "v_out": v_out,
        "eps": eps,
        "aop": aop,
        "comp_name_values": dic_comp_name_vals,
        "df_netlist": df_netlist,
        "df_unk_currents": df_netlist_unknown_currents,
        "dic_params": dic_params
    }
    
    return equ_smna, extra_results

