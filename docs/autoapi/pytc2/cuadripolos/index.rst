pytc2.cuadripolos
=================

.. py:module:: pytc2.cuadripolos

.. autoapi-nested-parse::

   Created on Thu Mar  2 14:14:13 2023

   @author: mariano



Attributes
----------

.. autoapisummary::

   pytc2.cuadripolos.opamp_models_str
   pytc2.cuadripolos.parametros_opamp
   pytc2.cuadripolos.posibles_entradas
   pytc2.cuadripolos.posibles_salidas


Functions
---------

.. autoapisummary::

   pytc2.cuadripolos.S2Ts_s
   pytc2.cuadripolos.Ts2S_s
   pytc2.cuadripolos.Ts2Tabcd_s
   pytc2.cuadripolos.Tabcd2S_s
   pytc2.cuadripolos.S2Tabcd_s
   pytc2.cuadripolos.Y2Tabcd_s
   pytc2.cuadripolos.Z2Tabcd_s
   pytc2.cuadripolos.Tabcd2Z_s
   pytc2.cuadripolos.Tabcd2Y_s
   pytc2.cuadripolos.I2Tabcd_s
   pytc2.cuadripolos.Model_conversion
   pytc2.cuadripolos.y2mai
   pytc2.cuadripolos.may2y
   pytc2.cuadripolos.Y2Tabcd
   pytc2.cuadripolos.Z2Tabcd
   pytc2.cuadripolos.Tabcd2Y
   pytc2.cuadripolos.I2Tabcd
   pytc2.cuadripolos.SparZ_s
   pytc2.cuadripolos.SparY_s
   pytc2.cuadripolos.TabcdLYZ_s
   pytc2.cuadripolos.TabcdLZY_s
   pytc2.cuadripolos.TabcdZ_s
   pytc2.cuadripolos.TabcdY_s
   pytc2.cuadripolos.TabcdLYZ
   pytc2.cuadripolos.TabcdLZY
   pytc2.cuadripolos.TabcdZ
   pytc2.cuadripolos.TabcdY
   pytc2.cuadripolos.calc_MAI_ztransf_ij_mn
   pytc2.cuadripolos.calc_MAI_vtransf_ij_mn
   pytc2.cuadripolos.calc_MAI_impedance_ij
   pytc2.cuadripolos.smna


Module Contents
---------------

.. py:function:: S2Ts_s(Spar)

   Convierte una matriz de parámetros de dispersión (S) simbólica
   en el modelo de parámetros de transferencia de dispersión (Ts).

   Esta función toma una matriz simbólica que representa los parámetros de dispersión (S) de un sistema y calcula la matriz de parámetros de transferencia de dispersión (Ts) correspondiente.


   :param Spar: Matriz de parámetros de dispersión S.
   :type Spar: Symbolic Matrix

   :returns: **Ts** -- Matriz de parámetros de transferencia de dispersión Ts.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Spar no es una instancia de Symbolic Matrix.
       Si Spar no tiene el formato correcto [ [Spar_11, Spar_12], [Spar_21, Spar_22] ].
       Si Spar_12 es nulo.

   .. seealso:: :func:`Ts2S_s`, :func:`S2Tabcd_s`, :func:`Model_conversion`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import S2Ts_s
   >>> Spar = sp.Matrix([[sp.symbols('S11'), sp.symbols('S12')],
   ...                   [sp.symbols('S21'), sp.symbols('S22')]])
   >>> Ts = S2Ts_s(Spar)
   >>> print(Ts)
   Matrix([[1/S21, -S22/S21], [S11/S21, -S11*S22/S21 + S12]])

   .. rubric:: Notes

   - La matriz Spar debe tener la forma [ [Spar_11, Spar_12], [Spar_21, Spar_22] ].
   - Spar_12 no puede ser nulo.
   - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.


.. py:function:: Ts2S_s(Ts)

   Convierte una matriz de transferencia de scattering (Ts) simbólica
   al modelo de parámetros scattering (S).

   :param Ts: Matriz de parámetros S.
   :type Ts: Symbolic Matrix

   :returns: **Spar** -- Matriz de parámetros de scattering.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Ts no es una instancia de Symbolic Matrix.
       Si Ts no tiene el formato correcto [ [Ts_11, Ts_12], [Ts_21, Ts_22] ].
       Si Ts_11 es nulo.

   .. seealso:: :func:`S2Ts_s`, :func:`S2Tabcd_s`, :func:`Model_conversion`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import Ts2S_s
   >>> Ts = sp.Matrix([[sp.symbols('Ts11'), sp.symbols('Ts12')],
   ...                 [sp.symbols('Ts21'), sp.symbols('Ts22')]])
   >>> Spar = Ts2S_s(Ts)
   >>> print(Spar)
   Matrix([[Ts21/Ts11, Ts22 - Ts12*Ts21/Ts11], [1/Ts11, -Ts12/Ts11]])

   .. rubric:: Notes

   - La matriz Ts debe tener la forma [ [Ts_11, Ts_12], [Ts_21, Ts_22] ].
   - Ts_11 no puede ser nulo.
   - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.


.. py:function:: Ts2Tabcd_s(Ts, Z01=sp.Rational('1'), Z02=sp.Rational('1'))

   Converts a symbolic scattering parameter matrix (Ts) to the symbolic ABCD or Tabcd model.

   This function converts a symbolic scattering parameter matrix (Ts) to the symbolic ABCD or Tabcd model.

   :param Ts: The Ts parameter matrix.
   :type Ts: Symbolic Matrix
   :param Z0: The reference impedance, defaults to 1.
   :type Z0: sp.Expr, optional

   :returns: **Tabcd** -- The ABCD parameter matrix.
   :rtype: Symbolic Matrix

   :raises ValueError: If Ts is not an instance of sp.Matrix.
       If Z0 is not an instance of sp.Expr.

   .. seealso:: :func:`Ts2S_s`, :func:`S2Tabcd_s`, :func:`Model_conversion`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import Ts2Tabcd_s
   >>> Z0 = sp.symbols('Z0')
   >>> Ts = sp.Matrix([[sp.symbols('Ts_11'), sp.symbols('Ts_12')],
   ...                 [sp.symbols('Ts_21'), sp.symbols('Ts_22')]])
   >>> Tabcd = Ts2Tabcd_s(Ts, Z0)
   >>> print(Tabcd)
   Matrix([[Ts_11/2 - Ts_12/2 - Ts_21/2 + Ts_22/2, Z0*(Ts_11 - Ts_12 + Ts_21 - Ts_22)/2], [(Ts_11 + Ts_12 - Ts_21 - Ts_22)/(2*Z0), Ts_11/2 - Ts_12/2 - Ts_21/2 + Ts_22/2]])


.. py:function:: Tabcd2S_s(Tabcd, Z01=sp.Rational('1'), Z02=sp.Rational('1'))

   Convierte una matriz de parámetros ABCD (Tabcd) simbólica
   al modelo de parámetros scattering (S).

   :param Tabcd: Matriz de parámetros ABCD.
   :type Tabcd: Symbolic Matrix
   :param Z0: Impedancia característica del medio. Por defecto es 1.
   :type Z0: sympy expression, optional

   :returns: **Spar** -- Matriz de parámetros de scattering.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Tabcd no es una instancia de Symbolic Matrix.
       Si Tabcd no tiene el formato correcto [ [A, B], [C, D] ].
       Si la matriz Tabcd no es invertible.

   .. seealso:: :func:`Ts2S_s`, :func:`S2Tabcd_s`, :func:`Model_conversion`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import Tabcd2S_s
   >>> Tabcd = sp.Matrix([[sp.symbols('A'), sp.symbols('B')],
   ...                    [sp.symbols('C'), sp.symbols('D')]])
   >>> Spar = Tabcd2S_s(Tabcd)
   >>> print(Spar)
   Matrix([[(A + B - C - D)/(A + B + C + D), 2*(A*D - B*C)/(A + B + C + D)], [2/(A + B + C + D), (-A + B - C + D)/(A + B + C + D)]])

   .. rubric:: Notes

   - La matriz Tabcd debe tener el formato [ [A, B], [C, D] ].
   - La matriz Tabcd debe ser invertible para realizar la conversión correctamente.
   - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.


.. py:function:: S2Tabcd_s(Spar, Z01=sp.Rational('1'), Z02=sp.Rational('1'))

   Convierte una matriz de parámetros scattering (S) simbólica
   al modelo de parámetros ABCD (Tabcd).

   :param Spar: Matriz de parámetros S.
   :type Spar: Symbolic Matrix
   :param Z0: Impedancia característica del medio. Por defecto es 1.
   :type Z0: sympy expression, optional

   :returns: **Tabcd** -- Matriz de parámetros ABCD.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Spar no es una instancia de Symbolic Matrix.
       Si Spar no tiene el formato correcto [ [S11, S12], [S21, S22] ].
       Si Spar[1, 0] es nulo.

   .. seealso:: :func:`Tabcd2S_s`, :func:`Y2Tabcd_s`, :func:`Model_conversion`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import S2Tabcd_s
   >>> Spar = sp.Matrix([[sp.symbols('S11'), sp.symbols('S12')],
   ...                   [sp.symbols('S21'), sp.symbols('S22')]])
   >>> Tabcd = S2Tabcd_s(Spar)
   >>> print(Tabcd)
   Matrix([[(-S11*S22 - S11 + S12*S21 + S22 + 1)/(2*S21), (S11*S22 + S11 - S12*S21 + S22 + 1)/(2*S21)], [(S11*S22 - S11 - S12*S21 - S22 + 1)/(2*S21), (-S11*S22 - S11 + S12*S21 + S22 + 1)/(2*S21)]])

   .. rubric:: Notes

   - La matriz Spar debe tener el formato [ [S11, S12], [S21, S22] ].
   - Spar[1, 0] no puede ser nulo.
   - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.


.. py:function:: Y2Tabcd_s(YY)

   Convierte una matriz de admitancia de dos puertos (YY) simbólica
   al modelo de parámetros ABCD (Tabcd).

   :param YY: Matriz de admitancia de dos puertos.
   :type YY: Symbolic Matrix

   :returns: **TT** -- Matriz de parámetros ABCD.
   :rtype: Symbolic Matrix

   :raises ValueError: Si YY no es una instancia de Symbolic Matrix.
       Si YY no tiene el formato correcto [ [Y11, Y12], [Y21, Y22] ].
       Si Y21 es nulo.

   .. seealso:: :func:`Ts2S_s`, :func:`Tabcd2Y_s`, :func:`Model_conversion`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import Y2Tabcd_s
   >>> YY = sp.Matrix([[sp.symbols('Y11'), sp.symbols('Y12')],
   ...                 [sp.symbols('Y21'), sp.symbols('Y22')]])
   >>> TT = Y2Tabcd_s(YY)
   >>> print(TT)
   Matrix([[-Y22/Y21, -1/Y21], [-(Y11*Y22 - Y12*Y21)/Y21, -Y22/Y21]])

   .. rubric:: Notes

   - La matriz YY debe tener el formato [ [Y11, Y12], [Y21, Y22] ].
   - YY[1, 0] no puede ser nulo.
   - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.


.. py:function:: Z2Tabcd_s(ZZ)

   Convierte la matriz de impedancia (ZZ) simbólica
   al modelo de parámetros ABCD (Tabcd).

   :param ZZ: Matriz de impedancia.
   :type ZZ: Symbolic Matrix

   :returns: **TT** -- Matriz de parámetros ABCD.
   :rtype: Symbolic Matrix

   :raises ValueError: Si ZZ no es una instancia de Symbolic Matrix.
       Si ZZ no tiene el formato correcto [ [Z11, Z12], [Z21, Z22] ].
       Si Z21 es nulo.

   .. seealso:: :func:`Tabcd2Z_s`, :func:`Tabcd2Y_s`, :func:`Model_conversion`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import Z2Tabcd_s
   >>> ZZ = sp.Matrix([[sp.symbols('Z11'), sp.symbols('Z12')],
   ...                 [sp.symbols('Z21'), sp.symbols('Z22')]])
   >>> TT = Z2Tabcd_s(ZZ)
   >>> print(TT)
   Matrix([[Z11/Z21, (Z11*Z22 - Z12*Z21)/Z21], [1/Z21, Z22/Z21]])

   .. rubric:: Notes

   - La matriz ZZ debe tener el formato [ [Z11, Z12], [Z21, Z22] ].
   - Z21 no puede ser nulo.
   - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.


.. py:function:: Tabcd2Z_s(TT)

   Convierte una matriz de parámetros ABCD (TT) simbólica
   al modelo de impedancia de dos puertos (ZZ).

   :param TT: Matriz de parámetros ABCD.
   :type TT: Symbolic Matrix

   :returns: **ZZ** -- Matriz de impedancia de dos puertos.
   :rtype: Symbolic Matrix

   :raises ValueError: Si TT no es una instancia de Symbolic Matrix.
       Si TT no tiene el formato correcto [ [A, B], [C, D] ].
       Si C es nulo.

   .. seealso:: :func:`Z2Tabcd_s`, :func:`Tabcd2Y_s`, :func:`Model_conversion`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import Tabcd2Z_s
   >>> TT = sp.Matrix([[sp.symbols('A'), sp.symbols('B')],
   ...                 [sp.symbols('C'), sp.symbols('D')]])
   >>> ZZ = Tabcd2Z_s(TT)
   >>> print(ZZ)
   Matrix([[A/C, (A*D - B*C)/C], [1/C, D/C]])

   .. rubric:: Notes

   - La matriz TT debe tener el formato [ [A, B], [C, D] ].
   - C no puede ser nulo.
   - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.


.. py:function:: Tabcd2Y_s(TT)

   Convierte una matriz de parámetros ABCD (TT) simbólica
   al modelo de admitancia de dos puertos (YY).

   :param TT: Matriz de parámetros ABCD.
   :type TT: Symbolic Matrix

   :returns: **YY** -- Matriz de admitancia de dos puertos.
   :rtype: Symbolic Matrix

   :raises ValueError: Si TT no es una instancia de Symbolic Matrix.
       Si TT no tiene el formato correcto [ [A, B], [C, D] ].
       Si B es nulo.

   .. seealso:: :func:`Y2Tabcd_s`, :func:`Tabcd2Z_s`, :func:`Model_conversion`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import Tabcd2Y_s
   >>> TT = sp.Matrix([[sp.symbols('A'), sp.symbols('B')],
   ...                 [sp.symbols('C'), sp.symbols('D')]])
   >>> YY = Tabcd2Y_s(TT)
   >>> print(YY)
   Matrix([[D/B, -(A*D - B*C)/B], [-1/B, A/B]])

   .. rubric:: Notes

   - La matriz TT debe tener el formato [ [A, B], [C, D] ].
   - B no puede ser nulo.
   - Esta función está diseñada para trabajar con matrices simbólicas utilizando el módulo SymPy.


.. py:function:: I2Tabcd_s(gamma, z01, z02=None)

   Convierte una ganancia compleja expresada en neppers (gamma)
   y la impedancia de referencia (z01,2) en una matriz de parámetros ABCD (TT).

   :param gamma: Ganancia compleja expresada en neppers (Re{gamma}) y radianes (Im{gamma}).
   :type gamma: Symbol
   :param z01: Impedancia de referencia del puerto 1.
   :type z01: Symbol
   :param z02: Impedancia de referencia del puerto 2. Si no se proporciona, se asume z02 = z01.
   :type z02: Symbol, opcional

   :returns: **TT** -- Matriz ABCD en función de los parámetros imagen.
   :rtype: Symbolic Matrix

   :raises ValueError: Si z01 no es un símbolo o no es un número real positivo.
       Si z02 no es un símbolo o no es un número real positivo.
       Si gamma no es un número complejo.

   .. seealso:: :func:`Y2Tabcd_s`, :func:`Tabcd2Z_s`, :func:`Model_conversion`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import I2Tabcd_s
   >>> gamma = sp.symbols('gamma')
   >>> z01 = sp.symbols('z01')
   >>> z02 = sp.symbols('z02')
   >>> TT = I2Tabcd_s(gamma, z01, z02)
   >>> print(TT)
   Matrix([[sqrt(z01/z02)*cosh(gamma), sqrt(z01*z02)*sinh(gamma)], [sinh(gamma)/sqrt(z01*z02), sqrt(z02/z01)*cosh(gamma)]])

   .. rubric:: Notes

   - Esta función está diseñada para trabajar con expresiones simbólicas utilizando el módulo SymPy.


.. py:function:: Model_conversion(src_model, dst_model)

   Convierte modelos de cuadripolos lineales de un formato a otro.

   :param src_model: Diccionario que describe el modelo de origen.
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
   :type src_model: dict
   :param dst_model: Diccionario que describe el modelo de salida.
                     Debe tener las mismas claves que src_model.
   :type dst_model: dict

   :returns: Diccionario que contiene la matriz convertida y el nombre del modelo resultante.
   :rtype: dict

   :raises ValueError: Si los modelos de origen y destino son iguales.
       Si falta alguna clave en src_model o dst_model.
       Si la variable independiente no es un símbolo o no es un número real positivo.

   .. seealso:: :func:`Y2Tabcd_s`, :func:`Tabcd2Z_s`, :func:`S2Ts_s`

   .. rubric:: Example

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

   .. rubric:: Notes

   - Esta función está diseñada para trabajar con expresiones simbólicas utilizando el módulo SymPy.


.. py:function:: y2mai(YY)

   Convierte una matriz de admitancia definida (YY) a una matriz admitancia indefinida (Ymai).

   :param YY: Matriz admitancia definida.
   :type YY: sympy.Matrix

   :returns: **Ymai** -- Matriz admitancia indefinida.
   :rtype: sympy.Matrix

   :raises ValueError: Si YY no es una instancia de sympy.Matrix.

   .. seealso:: :func:`may2y`, :func:`Y2Tabcd`, :func:`I2Tabcd`

   .. rubric:: Example

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import y2mai
   >>> YY = sp.Matrix([[sp.symbols('Y11'), sp.symbols('Y12')],
   ...                 [sp.symbols('Y21'), sp.symbols('Y22')]])
   >>> Ymai = y2mai(YY)
   >>> print(Ymai)
   Matrix([[Y11, Y12, -Y11 - Y12], [Y21, Y22, -Y21 - Y22], [-Y11 - Y21, -Y12 - Y22, Y11 + Y12 + Y21 + Y22]])

   .. rubric:: Notes

   - Esta función suma las corrientes de entrada y salida para obtener la matriz admitancia indefinida.
   - Se espera que YY sea una instancia de sympy.Matrix.


.. py:function:: may2y(Ymai, nodes2del)

   Convierte una matriz admitancia indefinida (Ymai) a una matriz admitancia (YY) luego de eliminar filas y columnas indicadas en nodes2del.

   :param Ymai: Matriz admitancia indefinida.
   :type Ymai: sympy.Matrix
   :param nodes2del: Índices de las filas y columnas que se eliminarán.
   :type nodes2del: list or integer

   :returns: **YY** -- Matriz admitancia.
   :rtype: sympy.Matrix

   :raises ValueError: Si Ymai no es una instancia de sympy.Matrix.
       Si nodes2del no es una lista o un entero.
       Si los elementos de nodes2del no son enteros o están fuera del rango de índices de Ymai.

   .. seealso:: :func:`y2mai`, :func:`Y2Tabcd`, :func:`I2Tabcd`

   .. rubric:: Example

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import may2y
   >>> Ymai = sp.Matrix([[sp.symbols('Y11'), sp.symbols('Y12'), sp.symbols('Y13')],
   ...                 [sp.symbols('Y21'), sp.symbols('Y22'), sp.symbols('Y23')],
   ...                 [sp.symbols('Y31'), sp.symbols('Y32'), sp.symbols('Y33')]])
   >>> nodes2del = [0, 2]
   >>> YY = may2y(Ymai, nodes2del)
   >>> print(YY)
   Matrix([[Y22]])

   .. rubric:: Notes

   - Esta función elimina las filas y columnas indicadas en nodes2del de Ymai para obtener la matriz admitancia YY.
   - Se espera que Ymai sea una instancia de sympy.Matrix.
   - nodes2del puede ser una lista de índices o un solo entero.
   - Los índices en nodes2del deben ser enteros y estar dentro del rango de índices de Ymai.


.. py:function:: Y2Tabcd(YY)

   Convierte una matriz de admitancia de dos puertos (YY) a la matriz de parámetros ABCD (TT).

   :param YY: Matriz de admitancia de dos puertos.
   :type YY: numpy.ndarray

   :returns: **TT** -- Matriz de parámetros ABCD.
   :rtype: numpy.ndarray

   :raises ValueError: Si YY no es una matriz de 2x2.
       Si Y21 es cero.

   .. seealso:: :func:`Z2Tabcd`, :func:`Tabcd2Y`, :func:`y2mai`

   .. rubric:: Example

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

   .. rubric:: Notes

   - Esta función asume que YY tiene el formato [ [Y11, Y12], [Y21, Y22] ].
   - YY[1, 0] no puede ser cero para evitar una división por cero.


.. py:function:: Z2Tabcd(ZZ)

   Convierte una matriz de impedancia de dos puertos (ZZ) a la matriz de parámetros ABCD (TT).

   :param ZZ: Matriz de impedancia de dos puertos.
   :type ZZ: numpy.ndarray

   :returns: **TT** -- Matriz de parámetros ABCD.
   :rtype: numpy.ndarray

   :raises ValueError: Si ZZ no es una matriz de 2x2.
       Si Z21 es cero.

   .. seealso:: :func:`Y2Tabcd`, :func:`Tabcd2Z`, :func:`may2y`

   .. rubric:: Example

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

   .. rubric:: Notes

   - Esta función asume que ZZ tiene el formato [ [Z11, Z12], [Z21, Z22] ].
   - ZZ[1, 0] no puede ser cero para evitar una división por cero.


.. py:function:: Tabcd2Y(TT)

   Convierte una matriz de parámetros ABCD (TT) a la matriz de admitancia de dos puertos (YY).

   :param TT: Matriz de parámetros ABCD.
   :type TT: numpy.ndarray

   :returns: **YY** -- Matriz de admitancia de dos puertos.
   :rtype: numpy.ndarray

   :raises ValueError: Si TT no es una matriz de 2x2.
       Si B es cero.

   .. seealso:: :func:`Y2Tabcd`, :func:`Tabcd2Z`, :func:`may2y`

   .. rubric:: Example

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

   .. rubric:: Notes

   - Esta función asume que TT tiene el formato [ [A, B], [C, D] ].
   - B no puede ser cero para evitar una división por cero.


.. py:function:: I2Tabcd(gamma, z01, z02=None)

   Convierte una ganancia compleja expresada en neppers (gamma)
   y la impedancia de referencia (z01,2) en una matriz de parámetros ABCD (TT).

   :param gamma: Ganancia compleja expresada en neppers (Re{gamma}) y radianes (Im{gamma}).
   :type gamma: float or complex
   :param z01: Impedancia de referencia del puerto 1.
   :type z01: float
   :param z02: Impedancia de referencia del puerto 2. Si no se proporciona, se asume z02 = z01.
   :type z02: float, opcional

   :returns: **TT** -- Matriz ABCD en función de los parámetros imagen.
   :rtype: numpy.ndarray

   :raises ValueError: Si z01 no es un número real positivo.
       Si z02 no es un número real positivo.

   .. seealso:: :func:`y2mai`, :func:`Tabcd2Y`, :func:`Y2Tabcd`

   .. rubric:: Examples

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


   .. rubric:: Notes

   - Esta función calcula la matriz de parámetros ABCD en función de una ganancia compleja gamma y las impedancias de referencia z01 y z02.
   - Si z02 no se proporciona, se asume que z02 = z01.
   - Se espera que z01 y z02 sean números reales positivos.


.. py:function:: SparZ_s(Zexc, Z01=sp.Rational(1), Z02=None)

   Convierte una matriz de transferencia de scattering (Ts) simbólica
   al modelo de parámetros scattering (S).

   :param Zexc: Función de excitación de la impedancia a representar.
   :type Zexc: sympy.Symbol
   :param Z01: Impedancia de referencia en el plano 1. Por defecto es 1.
   :type Z01: sympy.Symbol, optional
   :param Z02: Impedancia de referencia en el plano 2. Por defecto es 1.
   :type Z02: sympy.Symbol, optional

   :returns: **Spar** -- Matriz de parámetros de scattering de Z.
   :rtype: sympy.Matrix

   :raises ValueError: Si Zexc no es una instancia de Symbolic.
       Si Z01 no es una instancia de Symbolic.
       Si Z02 no es una instancia de Symbolic.

   .. seealso:: :func:`SparY_s`, :func:`TabcdLYZ_s`, :func:`TabcdZ_s`

   .. rubric:: Examples

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

   .. rubric:: Notes

   - Esta función está diseñada para trabajar con impedancias simbólicas utilizando el módulo SymPy.


.. py:function:: SparY_s(Yexc, Y01=sp.Rational('1'), Y02=None)

   Convierte una matriz de transferencia de scattering (Ts) simbólica
   al modelo de parámetros scattering (S).

   :param Yexc: Función de excitación de la admitancia a representar.
   :type Yexc: Symbolic impedance
   :param Y01: Admitancia de referencia en el plano 1. Por defecto es 1.
   :type Y01: Symbolic impedance, optional
   :param Y02: Admitancia de referencia en el plano 2. Por defecto es 1.
   :type Y02: Symbolic impedance, optional

   :returns: **Spar** -- Matriz de parámetros de scattering de Y.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Yexc no es una instancia de Symbolic.
       Si Y01 no es una instancia de Symbolic.
       Si Y02 no es una instancia de Symbolic.

   .. seealso:: :func:`SparZ_s`, :func:`TabcdLYZ_s`, :func:`TabcdLZY`

   .. rubric:: Examples

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

   .. rubric:: Notes

   - Esta función está diseñada para trabajar con admitancias simbólicas utilizando el módulo SymPy.


.. py:function:: TabcdLYZ_s(Yexc, Zexc)

   Implementa una matriz de transferencia ABCD (Tabcd) a partir de
   un cuadripolo constituido por una Y en derivación seguida  por
   una Z en serie.

   :param Yexc: Función de excitación de la admitancia a representar.
   :type Yexc: Symbolic admitance
   :param Zexc: Función de excitación de la impedancia a representar.
   :type Zexc: Symbolic impedance

   :returns: **Tabcd** -- Matriz de parámetros ABCD.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Yexc no es una instancia de Symbolic.
       Si Zexc no es una instancia de Symbolic.

   .. seealso:: :func:`SparZ_s`, :func:`TabcdZ`, :func:`TabcdLZY`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import TabcdLYZ_s
   >>> Y = sp.symbols('Y')
   >>> Z = sp.symbols('Z')
   >>> TT = TabcdLYZ_s(Y, Z)
   >>> print(TT)
   Matrix([[1, Z], [Y, Y*Z + 1]])


.. py:function:: TabcdLZY_s(Zexc, Yexc)

   Implementa una matriz de transferencia ABCD (Tabcd) a partir de
   un cuadripolo constituido por una Z en serie seguida de una Y en
   derivación.

   :param Zexc: Función de excitación de la impedancia a representar.
   :type Zexc: Symbolic impedance
   :param Yexc: Función de excitación de la admitancia a representar.
   :type Yexc: Symbolic admitance

   :returns: **Tabcd** -- Matriz de parámetros ABCD.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Zexc no es una instancia de Symbolic.
       Si Yexc no es una instancia de Symbolic.

   .. seealso:: :func:`SparZ_s`, :func:`TabcdLYZ_s`, :func:`TabcdY_s`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import TabcdLZY_s
   >>> Y = sp.symbols('Y')
   >>> Z = sp.symbols('Z')
   >>> TT = TabcdLZY_s(Z, Y)
   >>> print(TT)
   Matrix([[Y*Z + 1, Z], [Y, 1]])


.. py:function:: TabcdZ_s(Zexc)

   Implementa una matriz de transferencia ABCD (Tabcd) a partir de
   un cuadripolo constituido únicamente por una Z en serie.

   :param Zexc: Función de excitación de la impedancia a representar.
   :type Zexc: Symbolic impedance

   :returns: **Tabcd** -- Matriz de parámetros ABCD.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Zexc no es una instancia de Symbolic.

   .. seealso:: :func:`SparZ_s`, :func:`TabcdLYZ_s`, :func:`TabcdY_s`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import TabcdZ_s
   >>> Z = sp.symbols('Z')
   >>> TT = TabcdZ_s(Z)
   >>> print(TT)
   Matrix([[1, Z], [0, 1]])


.. py:function:: TabcdY_s(Yexc)

   Implementa una matriz de transferencia ABCD (Tabcd) a partir de
   un cuadripolo constituido únicamente por una Y en derivación.

   :param Yexc: Función de excitación de la admitancia a representar.
   :type Yexc: Symbolic admitance

   :returns: **Tabcd** -- Matriz de parámetros ABCD.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Yexc no es una instancia de Symbolic.

   .. seealso:: :func:`SparZ_s`, :func:`TabcdLYZ_s`, :func:`TabcdY_s`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.cuadripolos import TabcdY_s
   >>> Y = sp.symbols('Y')
   >>> TT = TabcdY_s(Y)
   >>> print(TT)
   Matrix([[1, 0], [Y, 1]])


.. py:function:: TabcdLYZ(Yexc, Zexc)

   Implementa una matriz de transferencia ABCD (Tabcd) a partir de
   un cuadripolo constituido por una Y en derivación seguida  por
   una Z en serie.

   :param Yexc: Función de excitación de la admitancia a representar.
   :type Yexc: Symbolic admitance
   :param Zexc: Función de excitación de la impedancia a representar.
   :type Zexc: Symbolic impedance

   :returns: **Tabcd** -- Matriz de parámetros ABCD.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Yexc no es una instancia de Symbolic.
       Si Zexc no es una instancia de Symbolic.

   .. rubric:: Examples

   >>> from pytc2.cuadripolos import TabcdLYZ
   >>> TT = TabcdLYZ(Yexc=2., Zexc=3.)
   >>> print(TT)
   [[1 3]
    [2 7]]

   >>> # Recordar la definición de la matriz como:
   ([[1, Z], [Y, Y*Z + 1]])


.. py:function:: TabcdLZY(Zexc, Yexc)

   Implementa una matriz de transferencia ABCD (Tabcd) a partir de
   un cuadripolo constituido por una Z en serie seguida una Y en
   derivación.

   :param Zexc: Función de excitación de la impedancia a representar.
   :type Zexc: Symbolic impedance
   :param Yexc: Función de excitación de la admitancia a representar.
   :type Yexc: Symbolic admitance

   :returns: **Tabcd** -- Matriz de parámetros ABCD.
   :rtype: Symbolic Matrix

   :raises ValueError: Si Zexc no es una instancia de Symbolic.
       Si Yexc no es una instancia de Symbolic.

   .. rubric:: Examples

   >>> from pytc2.cuadripolos import TabcdLZY
   >>> TT = TabcdLZY(Yexc=2., Zexc=3.)
   >>> print(TT)
   [[7. 3.]
    [2. 1.]]

   >>> # Recordar la definición de la matriz como:
   [[Y*Z + 1, Z], [Y, 1]]


.. py:function:: TabcdZ(Zexc)

   Implementa una matriz de transferencia ABCD (Tabcd) a partir de
   un cuadripolo constituido únicamente por una Z en serie.

   :param Zexc: Función de excitación de la impedancia a representar.
   :type Zexc: Symbolic impedance

   :returns: **Tabcd** -- Matriz de parámetros ABCD.
   :rtype: np.ndarray

   :raises ValueError: Si Zexc no es una instancia de Symbolic.

   .. rubric:: Examples

   >>> from pytc2.cuadripolos import TabcdZ
   >>> TT = TabcdZ(Zexc=3.)
   >>> print(TT)
   [[1. 3.]
    [0. 1.]]

   >>> # Recordar la definición de la matriz como:
   [[1, Z],
    [0, 1]]


.. py:function:: TabcdY(Yexc)

   Implementa una matriz de transferencia ABCD (Tabcd) a partir de
   un cuadripolo constituido únicamente por una Y en derivación.

   :param Yexc: Función de excitación de la admitancia a representar.
   :type Yexc: Symbolic admitance

   :returns: **Tabcd** -- Matriz de parámetros ABCD.
   :rtype: np.ndarray

   :raises ValueError: Si Yexc no es una instancia de Symbolic.

   .. rubric:: Examples

   >>> from pytc2.cuadripolos import TabcdY
   >>> TT = TabcdY(Yexc=2.)
   >>> print(TT)
   [[1. 0.]
    [2. 1.]]

   >>> # Recordar la definición de la matriz como:
   [[1, 0],
    [Y, 1]]


.. py:function:: calc_MAI_ztransf_ij_mn(Ymai, ii=2, jj=3, mm=0, nn=1, verbose=False)

   Calculates the impedance transfer V_ij / I_mn.

   This function calculates the impedance transfer V_ij / I_mn of a given
   multiport network represented by its admittance matrix.

   :param Ymai: The indefinite admittance matrix.
   :type Ymai: sp.Matrix
   :param ii: The index i of the output element, defaults to 2.
   :type ii: int, optional
   :param jj: The index j of the output element, defaults to 3.
   :type jj: int, optional
   :param mm: The index m of the input element, defaults to 0.
   :type mm: int, optional
   :param nn: The index n of the input element, defaults to 1.
   :type nn: int, optional
   :param verbose: If True, prints intermediate calculations, defaults to False.
   :type verbose: bool, optional

   :returns: **Tz** -- The impedance transfer.
   :rtype: sp.Expr

   :raises ValueError: If any of the indices is not an integer.
       If Ymai is not an instance of sp.Matrix.

   .. rubric:: Examples

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


.. py:function:: calc_MAI_vtransf_ij_mn(Ymai, ii=2, jj=3, mm=0, nn=1, verbose=False)

   Calculates the voltage transfer V_ij / V_mn.

   This function calculates the voltage transfer V_ij / V_mn of a given
   multiport network represented by its admittance matrix.

   :param Ymai: The indefinite admittance matrix.
   :type Ymai: sp.Matrix
   :param ii: The index i of the output element, defaults to 2.
   :type ii: int, optional
   :param jj: The index j of the output element, defaults to 3.
   :type jj: int, optional
   :param mm: The index m of the input element, defaults to 0.
   :type mm: int, optional
   :param nn: The index n of the input element, defaults to 1.
   :type nn: int, optional
   :param verbose: If True, prints intermediate calculations, defaults to False.
   :type verbose: bool, optional

   :returns: **Av** -- The voltage transfer.
   :rtype: sp.Expr

   :raises ValueError: If any of the indices is not an integer.
       If Ymai is not an instance of sp.Matrix.

   .. rubric:: Examples

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


.. py:function:: calc_MAI_impedance_ij(Ymai, ii=0, jj=1, verbose=False)

   Calculates the impedance transfer V_ij / V_mn.

   This function calculates the impedance transfer V_ij / V_mn of a given
   multiport network represented by its admittance matrix.

   :param Ymai: The indefinite admittance matrix.
   :type Ymai: sp.Matrix
   :param ii: The index i of the output element, defaults to 0.
   :type ii: int, optional
   :param jj: The index j of the output element, defaults to 1.
   :type jj: int, optional
   :param verbose: If True, prints intermediate calculations, defaults to False.
   :type verbose: bool, optional

   :returns: **ZZ** -- The impedance transfer.
   :rtype: sp.Expr

   :raises ValueError: If ii or jj is not an integer.
       If Ymai is not an instance of sp.Matrix.

   .. rubric:: Examples

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


.. py:data:: opamp_models_str
   :value: ['OA_ideal', 'OA_1polo', 'OA_integrador']


.. py:data:: parametros_opamp
   :value: ('aop', 'gbw', 'aol')


.. py:data:: posibles_entradas
   :value: ('v_v1', 'v_vi', 'v_vin')


.. py:data:: posibles_salidas
   :value: ('v_v2', 'v_vo', 'v_vout')


.. py:function:: smna(file_schematic, opamp_model='OA_ideal', bAplicarValoresComponentes=True, bAplicarParametros=True)

   Realiza el análisis nodal modificado (simbólico) al circuito definido en
   el archivo *file_schematic*. Los formatos aceptados son LTspice
   y Netlist (ver ejemplos).

   El código se basa en el [publicado](https://tiburonboy.github.io/Symbolic-Modified-Nodal-Analysis-using-Python/Introduction.html) por Tony a.k.a @tiburonboy.

   :param file_schematic: Los formatos de esquemáticos aceptados son de LTspice y Netlist.
   :type file_schematic: nombre de archivo del circuito
   :param opamp_model: Cómo se van a tratar a los OpAmps que haya en el esquemático. Las opciones son:

                       * 'OA_ideal': El opamp tiene ganancia y ancho de banda infinito.
                       * 'OA_1polo': El opamp tiene una ganancia con un solo polo situado en :math:`:math:`\frac{G_{BW}}{s+\frac{G_{BW}}{A_{OL}}}``.
                       * 'OA_integrador': El opamp se comporta como un integrador: :math:`\frac{G_{BW}}{s}`.
   :type opamp_model: string
   :param bAplicarValoresComponentes: Se aplicarán los valores de cada componente a la ecuación MNA.
   :type bAplicarValoresComponentes: bool, optional
   :param bAplicarParametros: Se aplicarán los valores de cada parámerto hallado en el esquemático a la ecuación MNA.
   :type bAplicarParametros: bool, optional

   :returns: * **report** (*text string*) -- The net list report.
             * **df_netlist** (*pandas dataframe*) -- circuit net list info loaded into a dataframe
             * **df_netlist_unknown_currents** (*pandas dataframe*) -- branches with unknown currents
             * **A** (*sympy matrix*) -- The A matrix is (m+n) by (m+n) and is the combination of 4 smaller matrices, G, B, C, and D.
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
             * **X** (*list*) -- The X matrix is an (n+m) by 1 vector that holds the unknown quantities (node voltages
               and the currents through the independent voltage sources). The top n elements are the n node
               voltages. The bottom m elements represent the currents through the m independent voltage
               sources in the circuit. The V matrix is n by 1 and holds the unknown voltages. The J matrix
               is m by 1 and holds the unknown currents through the voltage sources
             * **Z** (*list*) -- The Z matrix holds the independent voltage and current sources and is the combination
               of 2 smaller matrices I and Ev. The Z matrix is (m+n) by 1, n is the number of nodes,
               and m is the number of independent voltage sources. The I matrix is n by 1 and contains
               the sum of the currents through the passive elements into the corresponding node (either
               zero, or the sum of independent current sources). The Ev matrix is m by 1 and holds the
               values of the independent voltage sources.

   :raises ValueError: Si el archivo no existe o no es un nombre válido.

   .. rubric:: Examples

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


