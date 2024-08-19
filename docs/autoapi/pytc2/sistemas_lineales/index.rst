pytc2.sistemas_lineales
=======================

.. py:module:: pytc2.sistemas_lineales

.. autoapi-nested-parse::

   Created on Thu Mar  2 14:15:44 2023

   @author: mariano



Attributes
----------

.. autoapisummary::

   pytc2.sistemas_lineales.phase_change_thr


Functions
---------

.. autoapisummary::

   pytc2.sistemas_lineales.tfcascade
   pytc2.sistemas_lineales.tfadd
   pytc2.sistemas_lineales.sos2tf_analog
   pytc2.sistemas_lineales.tf2sos_analog
   pytc2.sistemas_lineales.zpk2sos_analog
   pytc2.sistemas_lineales.pretty_print_lti
   pytc2.sistemas_lineales.parametrize_sos
   pytc2.sistemas_lineales.pretty_print_bicuad_omegayq
   pytc2.sistemas_lineales.pretty_print_SOS
   pytc2.sistemas_lineales.analyze_sys
   pytc2.sistemas_lineales.pzmap
   pytc2.sistemas_lineales.group_delay
   pytc2.sistemas_lineales.GroupDelay
   pytc2.sistemas_lineales.bodePlot
   pytc2.sistemas_lineales.plot_plantilla
   pytc2.sistemas_lineales._nearest_real_complex_idx
   pytc2.sistemas_lineales._cplxreal
   pytc2.sistemas_lineales._one_sos2tf
   pytc2.sistemas_lineales._build_poly_str
   pytc2.sistemas_lineales._build_omegayq_str
   pytc2.sistemas_lineales._complementaryColor


Module Contents
---------------

.. py:data:: phase_change_thr

   Representa la máxima variación en una función de fase de muestra a muestra.
   Es útil para detectar cuando una función atraviesa un cero y se produce un
   salto de :math:`\pi` radianes. Por cuestiones de muestreo y de variaciones
   muy rápidas de fase, a veces conviene que sea un poco menor a :math:`\pi`.

.. py:function:: tfcascade(tfa, tfb)

   Realiza la cascada de dos funciones de transferencia.

   Esta función toma dos funciones de transferencia, 'tfa' y 'tfb', y calcula la función de transferencia resultante
   de su cascada. La cascada de dos funciones de transferencia es el producto de ambas en el dominio de Laplace.

   :param tfa: Primera función de transferencia.
   :type tfa: TransferFunction
   :param tfb: Segunda función de transferencia.
   :type tfb: TransferFunction

   :returns: Función de transferencia resultante de la cascada.
   :rtype: TransferFunction

   :raises ValueError: Si 'tfa' o 'tfb' no son instancias de TransferFunction.

   .. seealso:: :func:`tfadd`, :func:`pretty_print_lti`

   .. rubric:: Examples

   >>> from scipy.signal import TransferFunction
   >>> tfa = TransferFunction([1, 2], [1, 3, 2])
   >>> tfb = TransferFunction([1], [1, 4])
   >>> tfc = tfcascade(tfa, tfb)
   >>> print(tfc)
   TransferFunction([1, 2], [1, 7, 14, 8])


.. py:function:: tfadd(tfa, tfb)

   Suma dos funciones de transferencia.

   Esta función toma los coeficientes de dos funciones de transferencia, 'tfa' y 'tfb', y devuelve los coeficientes
   de la función de transferencia resultante que es la suma de las dos funciones de transferencia dadas. La función
   resultante se devuelve como un objeto de tipo TransferFunction.

   :param tfa: Coeficientes de la primera función de transferencia.
   :type tfa: TransferFunction
   :param tfb: Coeficientes de la segunda función de transferencia.
   :type tfb: TransferFunction

   :returns: Función de transferencia resultante.
   :rtype: TransferFunction

   :raises ValueError: Si 'tfa' o 'tfb' no son instancias de TransferFunction.

   .. seealso:: :func:`tfcascade`, :func:`pretty_print_lti`

   .. rubric:: Examples

   >>> from scipy.signal import TransferFunction
   >>> tfa = TransferFunction([1, 2], [3, 4])
   >>> tfb = TransferFunction([5, 6], [7, 8])
   >>> tfadd(tfa, tfb)


.. py:function:: sos2tf_analog(mySOS)

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


   :param mySOS: Matriz que define las secciones de segundo orden (SOS) del sistema.
   :type mySOS: array_like

   :returns: Función de transferencia analógica resultante.
   :rtype: TransferFunction

   :raises ValueError: Si 'mySOS' no es una matriz 2D o si las filas de la matriz no tienen exactamente 6 elementos.

   .. seealso:: :func:`tf2sos_analog`, :func:`pretty_print_SOS`, :func:`pretty_print_lti`

   .. rubric:: Examples

   >>> import numpy as np
   >>> from pytc2.sistemas_lineales import sos2tf_analog
   >>> mySOS = np.array([[1, 0.5, 1, 1, 0.2, 1],
   ...                    [1, 1, 1, 1, 1, 1]])
   >>> tf_analog = sos2tf_analog(mySOS)
   >>> print(tf_analog)
   TransferFunctionContinuous(
   array([1., 1., 1.]),
   array([1. , 1.5, 1.7, 1.5, 1.2, 1. ]),
   dt: None
   )


.. py:function:: tf2sos_analog(num, den=[])

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


   :param num: Coeficientes numéricos de la función de transferencia.
   :type num: array_like, TransferFunction
   :param den: Coeficientes denóminos de la función de transferencia.
   :type den: array_like, opcional

   :returns: Matriz que define las secciones de segundo orden (SOS) del sistema analógico.
   :rtype: array_like

   :raises ValueError: Si 'num' o 'den' no son instancias de arrays de numpy.

   .. seealso:: :func:`sos2tf_analog`, :func:`pretty_print_SOS`, :func:`pretty_print_lti`

   .. rubric:: Examples

   >>> from pytc2.sistemas_lineales import tf2sos_analog
   >>> num = [1, 2, 3]
   >>> den = [4, 5, 6]
   >>> sos_analog = tf2sos_analog(num, den)
   >>> print(sos_analog)
   [[1. 2. 3. 4. 5. 6.]]


.. py:function:: zpk2sos_analog(zz, pp, kk)

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

   El algoritmo utilizado para convertir de ZPK a formato SOS sigue las sugerencias del libro :ref:`Design of Analog Filters de R. Schaumann <schau13>` , Cap. 5:

   1. Asignar ceros a los polos más cercanos.
   2. Ordenar las secciones por Q creciente.
   3. Ordenar las ganancias para maximizar el rango dinámico. Ver :ref:`Schaumann R. <schau13>` cap. 5.


   :param zz: Zeros de la función de transferencia.
   :type zz: array_like
   :param pp: Polos de la función de transferencia.
   :type pp: array_like
   :param kk: Ganancia del sistema.
   :type kk: float

   :returns: Matriz que define las secciones de segundo orden (SOS) del sistema analógico.
   :rtype: array_like

   :raises AssertionError: Si hay más ceros que polos.
   :raises ValueError: Si la factorización de la función de transferencia es incorrecta.

   .. seealso:: :func:`sos2tf_analog`, :func:`pretty_print_SOS`, :func:`pretty_print_lti`, :mod:`scipy.signal`

   .. rubric:: Examples

   >>> from pytc2.sistemas_lineales import zpk2sos_analog
   >>> zz = [1, 2, 3]
   >>> pp = [4, 5, 6]
   >>> kk = 2.5
   >>> sos_analog = zpk2sos_analog(zz, pp, kk)
   >>> print(sos_analog)
   [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]

   Notes:
   -------
   .. _schau13:

   Schaumann, Rolf, Haiqiao Xiao, and Van Valkenburg Mac. Design of analog filters 2nd. Edition. Oxford University Press, 2013. ISBN   0195373944, 9780195373943.



.. py:function:: pretty_print_lti(num, den=None, displaystr=True)

   Genera una representación matemática de una función de transferencia lineal en función de sus coeficientes numéricos.

   Esta función toma los coeficientes del numerador y, opcionalmente, los del denominador de una función de transferencia lineal y genera una representación matemática en función de estos coeficientes. Los parámetros opcionales permiten especificar si se debe mostrar o devolver la cadena formateada.


   :param num: Coeficientes del numerador de la función de transferencia.
   :type num: array_like, lista o TransferFunction
   :param den: Coeficientes del denominador de la función de transferencia. Por defecto es None.
   :type den: array_like, opcional
   :param displaystr: Indica si mostrar el resultado como salida o devolverlo como una cadena de texto. Por defecto es True.
   :type displaystr: bool, opcional

   :returns: Si displaystr es True, muestra la cadena formateada, si no, devuelve la cadena.
   :rtype: None or str

   :raises ValueError: Si los coeficientes numéricos no son de tipo array_like, lista, o un objeto TransferFunction.
       Si los coeficientes del denominador son proporcionados pero no son de tipo array_like.
       Si el argumento displaystr no es de tipo bool.

   .. seealso:: :func:`pretty_print_bicuad_omegayq`, :func:`_build_poly_str`

   .. rubric:: Examples

   >>> from pytc2.sistemas_lineales import pretty_print_lti
   >>> num = [1, 2, 3]
   >>> den = [4, 5, 6]
   >>> pretty_print_lti(num, den)
   [Devuelve la cadena formateada en LaTex de la función de transferencia]


.. py:function:: parametrize_sos(num, den=sp.Rational(1))

   Parametriza una función de transferencia de segundo orden en función de sus coeficientes.


   :param num: Coeficientes del numerador.
   :type num: Poly
   :param den: Coeficientes del denominador.
   :type den: Poly

   :returns:

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
   :rtype: tuple

   :raises ValueError: Si los coeficientes numéricos no son de tipo Poly.
       Si los coeficientes del denominador no son proporcionados o no son de tipo Poly.

   .. seealso:: :func:`pretty_print_bicuad_omegayq`, :func:`pretty_print_SOS`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.sistemas_lineales import parametrize_sos
   >>> from pytc2.general import s, print_latex, a_equal_b_latex_s
   >>> a, b, c, d, e , f = sp.symbols('a, b, c, d, e , f', real=True, positive=True)
   >>> num = sp.Poly((a*s + b),s)
   >>> den = sp.Poly((c*s + d),s)
   >>> num_bili1, den_bili1, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)
   >>> ￼print(num_bili1)
   Poly(s + b/a, s, domain='ZZ(a,b)')
   >>> ￼print(den_bili1)
   Poly(s + d/c, s, domain='ZZ(c,d)')
   >>> num = sp.Poly((a*s**2 + b*s + c),s)
   >>> den = sp.Poly((d*s**2 + e*s + f),s)
   >>> num_sos1, den_sos1, w_on, Q_n, w_od, Q_d, K = parametrize_sos(num, den)
   >>> print(w_on)
   sqrt(c)/sqrt(a)
   >>> print(Q_n)
   sqrt(a)*sqrt(c)/b


.. py:function:: pretty_print_bicuad_omegayq(num, den=None, displaystr=True)

   Genera una representación matemática de un sistema de segundo orden en función de su frecuencia natural (omega) y su factor de calidad (Q).

   Esta función toma los coeficientes del numerador y, opcionalmente, los del denominador de un sistema de segundo orden y genera una representación matemática en función de la frecuencia natural (omega) y el factor de calidad (Q). Los parámetros opcionales permiten especificar si se debe mostrar o devolver la cadena formateada.


   :param num: Los coeficientes del numerador del sistema de segundo orden.
   :type num: array_like
   :param den: Los coeficientes del denominador del sistema de segundo orden. Por defecto es None.
   :type den: array_like, opcional
   :param displaystr: Indica si mostrar el resultado como salida o devolverlo como una cadena de texto. Por defecto es True.
   :type displaystr: bool, opcional

   :returns: Si displaystr es True, muestra la cadena formateada, si no, devuelve la cadena.
   :rtype: None or str

   :raises ValueError: Si los coeficientes numéricos no son proporcionados.
       Si los coeficientes numéricos no tienen una longitud de 3 elementos.

   .. seealso:: :func:`pretty_print_SOS`, :func:`_build_omegayq_str`

   .. rubric:: Examples

   >>> from pytc2.sistemas_lineales import pretty_print_bicuad_omegayq
   >>> pretty_print_bicuad_omegayq([1, 2, 1], [1, 1, 1])
   [ Expresión formateada en LaTex del sistema de segundo orden]


.. py:function:: pretty_print_SOS(mySOS, mode='default', displaystr=True)

   Imprime de forma "bonita" una expresión que define a un sistema de segundo orden (SOS)

   Esta función toma una matriz que define las secciones de segundo orden (SOS) y muestra la representación matemática de la cadena de sistemas de segundo orden. Los parámetros opcionales permiten especificar el modo de impresión y si se debe mostrar o devolver la cadena formateada.

   Los SOS siempre deben definirse como::

       mySOS= ( [ a1_1 a2_1 a3_1 b1_1 b2_1 b3_1 ]
                [ a1_2 a2_2 a3_2 b1_2 b2_2 b3_2 ]
                ...
                [ a1_N a2_N a3_N b1_N b2_N b3_N ]
               )

   donde cada sección o línea de `mySOS` significa matemáticamente

   .. math:: T_i = \frac{a_{1i} \, s^2 + a_{2i} \, s + a_{3i}}{b_{1i} \, s^2 + b_{2i} \, s + b_{3i}}


   :param mySOS: La matriz que define los coeficientes de las secciones de segundo orden.
   :type mySOS: matriz numpy
   :param mode: El modo de impresión. Puede ser 'default' o 'omegayq'. Por defecto es 'default'.
   :type mode: str, opcional
   :param displaystr: Indica si mostrar el resultado como salida o devolverlo como una cadena de texto. Por defecto es True.
   :type displaystr: bool, opcional

   :returns: Si displaystr es True, muestra la cadena formateada, si no, devuelve la cadena.
   :rtype: None or str

   :raises ValueError: Si el modo de impresión no es válido.
       Si mySOS no es una matriz numpy.
       Si mySOS no tiene exactamente 6 columnas.
       Si displaystr no es un booleano.

   .. seealso:: :func:`parametrize_sos`, :func:`pretty_print_lti`, :func:`pretty_print_bicuad_omegayq`

   .. rubric:: Examples

   >>> import numpy as np
   >>> from pytc2.sistemas_lineales import pretty_print_SOS
   >>> mySOS = np.array([[1, 2, 1, 1, 1, 1], [1, 3, 1, 1, 4, 1]])
   >>> pretty_print_SOS(mySOS)
   [ Expresión formateada en LaTex de las SOS ]


.. py:function:: analyze_sys(all_sys, sys_name=None, img_ext='none', same_figs=True, annotations=True, xaxis='omega', fs=2 * np.pi)

   Analiza el comportamiento de un sistema lineal en términos de:

         * Respuesta de magnitud y fase o gráfico de Bode
         * Mapa de polos y ceros
         * Retardo de grupo

       La función admite el sistema a analizar (*all_sys*) como:

           * uno o una lista de objetos TransferFunction
           * una matriz que define varias secciones de segundo orden (SOS).

       Si *all_sys* es una matriz SOS, la función muestra cada una de las SOS
       y el sistema resultante de la cascada de todas las SOS.

   Esta función toma un sistema lineal (ya sea una lista de objetos
   TransferFunction o una matriz que define una cascada de secciones de segundo
   orden) y realiza un análisis completo del comportamiento del sistema,
   incluyendo trazado de gráficos de Bode, mapa de polos y ceros, y gráfico de
   retardo de grupo. Los parámetros opcionales permiten personalizar el análisis
   según las necesidades del usuario.


   :param all_sys: El sistema lineal a analizar como objeto/s *TransferFunction*. Ya sea una
                   lista de objetos TransferFuncion de scipy.signal o una matriz que define
                   una cascada de SOS.
   :type all_sys: TransferFunction o lista, tupla de TransferFunction o matriz numérica (Nx6)
   :param sys_name: Las etiquetas o descripción del sistema. Por defecto es None.
   :type sys_name: string o lista, opcional
   :param img_ext: Cuando es diferente de 'none', la función guarda los resultados del
                   gráfico en un archivo con la extensión indicada. Por defecto es 'none'.
   :type img_ext: string ['none', 'png', 'svg'], opcional
   :param same_figs: Usa siempre los mismos números de figura para trazar resultados.
                     Cuando es False, cada llamada produce un nuevo grupo de figuras en un
                     contenedor de gráficos separado. Por defecto es True.
   :type same_figs: booleano, opcional
   :param annotations: Agrega anotaciones al gráfico del mapa de polos y ceros. Cuando es True,
                       cada singularidad estará acompañada del valor de omega (es decir, la
                       distancia radial al origen) y Q (es decir, una medida de proximidad al
                       eje jw). Por defecto es True.
   :type annotations: booleano, opcional
   :param xaxis: El significado del eje X: "omega" se mide en radianes/s y se prefiere
                 para sistemas analógicos. "freq" se mide en Hz (1/s) y es válido tanto
                 para sistemas digitales como analógicos. "norm" es una versión
                 normalizada con la norma definida por fs. Por defecto es 'omega'.
   :type xaxis: string, opcional ['omega', 'freq', 'norm']
   :param fs: La frecuencia de muestreo del sistema digital o la norma para xaxis
              igual a "norm". Solo es válido si digital es True. Por defecto es None
              (definido en 1/dlti.dt).
   :type fs: valor real, opcional

   :raises ValueError: Si la extensión de imagen no es válida.
       Si sys_name no es una lista o un string.
       Si all_sys no es una lista o una matriz.
       Si xaxis no es válido.

   :returns: **return_values** -- Lista con tres pares de manijas de figuras y ejes de cada gráfico
             mostrado.
   :rtype: lista

   .. seealso:: :func:`pretty_print_bicuad_omegayq`, :func:`bodePlot`, :func:`pzmap`

   .. rubric:: Examples

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
   [ Expresión formateada en LaTex ]
   >>> analyze_sys([H1], sys_name='mi Ejemplo')
   [ Tres gráficas: respuesta en frec (mód, fase y retardo) y pzmap ]
   >>> # Compara el sistema anterior con otros dos con valores diferentes de Q
   >>> Q = 5
   >>> w0 = 1
   >>> num = np.array([w0**2])
   >>> den = np.array([1., w0 / Q, w0**2])
   >>> H2 = sig.TransferFunction(num, den)
   >>> analyze_sys([H1, H2], sys_name=['H1', 'H2'])


.. py:function:: pzmap(myFilter, annotations=False, filter_description=None, fig_id='none', axes_hdl='none', digital=False, fs=2 * np.pi)

   Grafica el mapa de polos y ceros de un filtro dado.

   :param myFilter: Objeto del filtro.
   :type myFilter: LTI object
   :param annotations: Indica si se deben añadir anotaciones a los polos y ceros.
                       El valor predeterminado es False.
   :type annotations: bool, opcional
   :param filter_description: Descripción del filtro. El valor predeterminado es None.
   :type filter_description: str, opcional
   :param fig_id: Identificador de la figura. Si se establece en 'none', se creará una nueva figura.
                  El valor predeterminado es 'none'.
   :type fig_id: str or int, opcional
   :param axes_hdl: Identificador o handle  del eje. Si se establece en 'none', se utilizará el eje actual.
                    El valor predeterminado es 'none'.
   :type axes_hdl: str or axes handle, opcional
   :param digital: Indica si el filtro es digital. El valor predeterminado es False.
   :type digital: bool, opcional
   :param fs: Frecuencia de muestreo. El valor predeterminado es 2*pi.
   :type fs: float, opcional

   :returns: * **fig_id** (*int*) -- Identificador de la figura creada.
             * **axes_hdl** (*axes handle*) -- handle  del eje utilizado para el gráfico.

   :raises ValueError: Si `fig_id` no es un string o un entero.
       Si `axes_hdl` no es un string o una handle  de eje válida.
       Si `digital` no es un booleano.
       Si `fs` no es un valor numérico.

   .. seealso:: :func:`analyze_sys`, :func:`bodePlot`, :func:`pzmap`

   .. rubric:: Examples

   >>> # Analiza un sistema con w0 = 1 rad/s y Q = sqrt(2)/2
   >>> import numpy as np
   >>> from scipy import signal as sig
   >>> from pytc2.sistemas_lineales import pzmap
   >>> Q = np.sqrt(2)/2
   >>> w0 = 1
   >>> num = np.array([w0**2])
   >>> den = np.array([1., w0 / Q, w0**2])
   >>> H1 = sig.TransferFunction(num, den)
   >>> fig_id, ax_hdl = pzmap(H1, annotations=True, filter_description='Filtro Pasabajos')


.. py:function:: group_delay(freq, phase)

   Calcula el retardo de grupo para una función de fase.


   :param freq: La grilla de frecuencia a la que se calcula la fase.
   :type freq: array_like
   :param phase: La fase de la función para la cual se calcula el retardo de grupo.
   :type phase: array_like

   :returns: **gd** -- Estimación del retardo de grupo, que es la derivada de la fase
             respecto a la frecuencia cambiada de signo.
   :rtype: array_like

   :raises ValueError: Si `freq` y `phase` no tienen la misma longitud.
       Si `freq` o `phase` no son arreglos NumPy.

   .. rubric:: Examples

   >>> from pytc2.sistemas_lineales import group_delay
   >>> import numpy as np
   >>> freq = np.linspace(0, 10, 10)
   >>> phase = np.sin(freq)
   >>> group_delay(freq, phase)
   array([-0.80657298,  0.09087493,  0.88720922,  0.69637424, -0.26929404,
          -0.93532747, -0.56065199,  0.43784299,  0.94916411,  0.94916411])


.. py:function:: GroupDelay(myFilter, fig_id='none', filter_description=None, npoints=1000, digital=False, xaxis='omega', unwrap_phase=False, fs=2 * np.pi)

   Calcula y grafica el retardo de grupo de un filtro.


   :param myFilter: Coeficientes del filtro o objeto TransferFunction del filtro.
   :type myFilter: array_like o scipy.signal.TransferFunction
   :param fig_id: Identificador de la figura. Si es 'none', crea una nueva figura. Por defecto es 'none'.
   :type fig_id: str o int, opcional
   :param filter_description: Descripción del filtro. Por defecto es None.
   :type filter_description: str, opcional
   :param npoints: Número de puntos para muestrear el eje de frecuencia. Por defecto es 1000.
   :type npoints: int, opcional
   :param digital: Indicador de si el filtro es digital. Por defecto es False.
   :type digital: bool, opcional
   :param xaxis: Tipo de eje x ('omega', 'freq', 'norm'). Por defecto es 'omega'.
   :type xaxis: str, opcional
   :param unwrap_phase: Evita que la respuesta de fase tenga saltos, habitualmente producidos
                        al haber ceros sobre el eje j.omega o la circunsferencia unitaria.
                        Por defecto es False.
   :type unwrap_phase: bool, opcional
   :param fs: Frecuencia de muestreo. Por defecto es 2*pi.
   :type fs: float, opcional

   :returns: * **fig_id** (*int*) -- Identificador de la figura.
             * **axes_hdl** (*Axes*) -- Manejador de ejes de la figura.

   :raises ValueError: Si myFilter no es un array NumPy ni un objeto TransferFunction.
       Si fig_id no es de tipo str, int o None.
       Si npoints no es un entero.
       Si digital no es un booleano.
       Si xaxis no es uno de los valores permitidos: 'omega', 'freq', 'norm'.
       Si fs no es un número.

   .. seealso:: :func:`analyze_sys`, :func:`bodePlot`, :func:`pzmap`

   .. rubric:: Example

   >>> # Analiza un sistema con w0 = 1 rad/s y Q = sqrt(2)/2
   >>> import numpy as np
   >>> from scipy import signal as sig
   >>> from pytc2.sistemas_lineales import GroupDelay
   >>> Q = np.sqrt(2)/2
   >>> w0 = 1
   >>> num = np.array([w0**2])
   >>> den = np.array([1., w0 / Q, w0**2])
   >>> H1 = sig.TransferFunction(num, den)
   >>> fig_id, axes_hdl = GroupDelay(H1, fig_id=1, filter_description='Filtro pasa bajos', npoints=1000, digital=False, xaxis='omega', fs=2*np.pi)


.. py:function:: bodePlot(myFilter, fig_id='none', axes_hdl='none', filter_description=None, npoints=1000, digital=False, xaxis='omega', unwrap_phase=False, fs=2 * np.pi)

   Grafica el diagrama de Bode (magnitud y fase) de un filtro.


   :param myFilter: Coeficientes del filtro o objeto TransferFunction del filtro.
   :type myFilter: array_like o scipy.signal.TransferFunction
   :param fig_id: Identificador de la figura. Si es 'none', crea una nueva figura. Por defecto es 'none'.
   :type fig_id: str o int, opcional
   :param axes_hdl: Manejador de ejes de la figura. Si es 'none', crea nuevos ejes. Por defecto es 'none'.
   :type axes_hdl: str o array_like de Axes, opcional
   :param filter_description: Descripción del filtro. Por defecto es None.
   :type filter_description: str, opcional
   :param npoints: Número de puntos para muestrear el eje de frecuencia. Por defecto es 1000.
   :type npoints: int, opcional
   :param digital: Indicador de si el filtro es digital. Por defecto es False.
   :type digital: bool, opcional
   :param xaxis: Tipo de eje x ('omega', 'freq', 'norm'). Por defecto es 'omega'.
   :type xaxis: str, opcional
   :param unwrap_phase: Evita que la respuesta de fase tenga saltos, habitualmente producidos
                        al haber ceros sobre el eje j.omega o la circunsferencia unitaria.
                        Por defecto es False.
   :type unwrap_phase: bool, opcional
   :param fs: Frecuencia de muestreo. Por defecto es 2*pi.
   :type fs: float, opcional

   :returns: * **fig_id** (*int*) -- Identificador de la figura.
             * **axes_hdl** (*array_like de Axes*) -- Manejadores de ejes de la figura.

   :raises ValueError: Si myFilter no es un array NumPy ni un objeto TransferFunction.
       Si los argumentos fig_id o axes_hdl no son válidos.
       Si xaxis no es uno de los valores permitidos: 'omega', 'freq', 'norm'.

   .. seealso:: :func:`analyze_sys`, :func:`GroupDelay`, :func:`pzmap`

   .. rubric:: Example

   >>> # Analiza un sistema con w0 = 1 rad/s y Q = sqrt(2)/2
   >>> import numpy as np
   >>> from scipy import signal as sig
   >>> from pytc2.sistemas_lineales import bodePlot
   >>> Q = np.sqrt(2)/2
   >>> w0 = 1
   >>> num = np.array([w0**2])
   >>> den = np.array([1., w0 / Q, w0**2])
   >>> H1 = sig.TransferFunction(num, den)
   >>> fig_id, axes_hdl = bodePlot(H1, fig_id=1, axes_hdl='none', filter_description='Filtro pasa bajos', npoints=1000, digital=False, xaxis='omega', fs=2*np.pi)


.. py:function:: plot_plantilla(filter_type='', fpass=0.25, ripple=0.5, fstop=0.6, attenuation=40, fs=2)

   Plotea una plantilla de diseño de filtro digital.

   :param filter_type: Tipo de filtro ('lowpass', 'highpass', 'bandpass', 'bandstop'). Por defecto es 'lowpass'.
   :type filter_type: str, opcional
   :param fpass: Frecuencia de paso o tupla de frecuencias de paso para los filtros 'bandpass' o 'bandstop'.
   :type fpass: float o tupla, opcional
   :param ripple: Máxima ondulación en la banda de paso (en dB). Por defecto es 0.5 dB.
   :type ripple: float, opcional
   :param fstop: Frecuencia de detención o tupla de frecuencias de detención para los filtros 'bandpass' o 'bandstop'.
   :type fstop: float o tupla, opcional
   :param attenuation: Atenuación mínima en la banda de detención (en dB). Por defecto es 40 dB.
   :type attenuation: float, opcional
   :param fs: Frecuencia de muestreo. Por defecto es 2.
   :type fs: float, opcional

   :rtype: None

   :raises ValueError: Si los argumentos no son del tipo o valor correcto.

   .. seealso:: :func:`analyze_sys`

   .. rubric:: Example

   >>> # Analiza un sistema con w0 = 1 rad/s y Q = sqrt(2)/2
   >>> import numpy as np
   >>> from scipy import signal as sig
   >>> import matplotlib.pyplot as plt
   >>> from pytc2.sistemas_lineales import bodePlot, plot_plantilla
   >>> Q = np.sqrt(2)/2
   >>> w0 = 1
   >>> num = np.array([w0**2])
   >>> den = np.array([1., w0 / Q, w0**2])
   >>> H1 = sig.TransferFunction(num, den)
   >>> fig_id, axes_hdl = bodePlot(H1, fig_id=1, axes_hdl='none', filter_description='Filtro pasa bajos', npoints=1000, digital=False, xaxis='omega', fs=2*np.pi)
   >>> plt.sca(axes_hdl[0])
   >>> plot_plantilla(filter_type='lowpass', fpass=1.0, ripple=3, fstop=3.0, attenuation=20, fs=2)


.. py:function:: _nearest_real_complex_idx(fro, to, which)

   Obtiene el índice del siguiente elemento real o complejo más cercano basado en la distancia.

   :param fro: Arreglo de partida que contiene los elementos a comparar.
   :type fro: array_like
   :param to: Valor de referencia para encontrar el elemento más cercano.
   :type to: array_like
   :param which: Especifica si se busca el elemento real o complejo más cercano.
   :type which: {'real', 'complex'}

   :returns: Índice del elemento más cercano en el arreglo de partida.
   :rtype: int

   :raises AssertionError: Si el argumento 'which' no es 'real' o 'complex'.

   .. seealso:: :func:`zpk2sos_analog`

   .. rubric:: Example

   >>> import numpy as np
   >>> from pytc2.sistemas_lineales import _nearest_real_complex_idx
   >>> fro = np.array([1, 2, 3, 4])
   >>> to = 2.5
   >>> nearest_idx = _nearest_real_complex_idx(fro, to, 'real')
   >>> print("El índice del elemento real más cercano a", to, "es:", nearest_idx)


.. py:function:: _cplxreal(z, tol=None)

   Separa en partes complejas y reales, combinando pares conjugados.

   El vector de entrada unidimensional `z` se divide en sus elementos complejos (`zc`) y reales (`zr`).
   Cada elemento complejo debe ser parte de un par conjugado complejo, que se combinan en un solo número
   (con parte imaginaria positiva) en la salida. Dos números complejos se consideran un par conjugado si sus partes
   real e imaginaria difieren en magnitud por menos de ``tol * abs(z)``.

   :param z: Vector de números complejos para ordenar y dividir.
   :type z: array_like
   :param tol: Tolerancia relativa para probar la realidad e igualdad conjugada.
               El valor predeterminado es ``100 * espaciado(1)`` del tipo de datos de `z`
               (es decir, 2e-14 para float64).
   :type tol: float, opcional

   :returns: * **zc** (*ndarray*) -- Elementos complejos de `z`, donde cada par se representa por un solo valor
               con parte imaginaria positiva, ordenada primero por parte real y luego
               por magnitud de la parte imaginaria. Los pares se promedian cuando se combinan
               para reducir el error.
             * **zr** (*ndarray*) -- Elementos reales de `z` (aquellos que tienen parte imaginaria menor que
               `tol` veces su magnitud), ordenados por valor.

   :raises ValueError: Si hay números complejos en `z` para los cuales no se puede encontrar un conjugado.

   .. seealso:: :func:`zpk2sos_analog`

   Exampless
   ---------
   >>> import numpy as np
   >>> from pytc2.sistemas_lineales import _cplxreal
   >>> a = [4, 3, 1, 2-2j, 2+2j, 2-1j, 2+1j, 2-1j, 2+1j, 1+1j, 1-1j]
   >>> zc, zr = _cplxreal(a)
   >>> print(zc)
   [ 1.+1.j  2.+1.j  2.+1.j  2.+2.j]
   >>> print(zr)
   [ 1.  3.  4.]




.. py:function:: _one_sos2tf(mySOS)

   Convierte una sección de segundo orden (SOS) en coeficientes de función de transferencia.

   :param mySOS: Vector que define una sección de segundo orden (SOS) del sistema.
   :type mySOS: array_like

   :returns: * **num** (*ndarray*) -- Coeficientes del numerador de la función de transferencia.
             * **den** (*ndarray*) -- Coeficientes del denominador de la función de transferencia.

   :raises ValueError: Si la entrada no es un vector con al menos 6 elementos.

   .. seealso:: :func:`sos2tf_analog`

   .. rubric:: Examples

   >>> import numpy as np
   >>> from pytc2.sistemas_lineales import _one_sos2tf
   >>> mySOS = [1, -1.9, 1, 1, -1.6, 0.64]
   >>> num, den = _one_sos2tf(mySOS)
   >>> print(num)
   [1, -1.9, 1]
   >>> print(den)
   [1, -1.6, 0.64]


.. py:function:: _build_poly_str(this_poly)

   Construye una cadena de caracteres que representa un polinomio.

   :param this_poly: Coeficientes del polinomio.
   :type this_poly: ndarray

   :returns: Cadena de caracteres que representa el polinomio.
   :rtype: str

   :raises ValueError: Si `this_poly` no es un array de numpy.

   .. seealso:: :func:`pretty_print_lti`, :func:`pretty_print_bicuad_omegayq`

   .. rubric:: Examples

   >>> import numpy as np
   >>> from pytc2.sistemas_lineales import _build_poly_str
   >>> this_poly = np.array([1, -2, 3])
   >>> poly_str = _build_poly_str(this_poly)
   >>> print(poly_str)
   's^2 - 2 s + 3'


.. py:function:: _build_omegayq_str(this_quad_poly, den=np.array([]))

   Construye una cadena de caracteres que representa un polinomio parametrizado
   mediante :math:`\omega_0` y Q.

   :param this_quad_poly: Coeficientes del polinomio cuadrático.
   :type this_quad_poly: ndarray
   :param den: Coeficientes del denominador. El valor predeterminado es np.array([]).
   :type den: ndarray, opcional

   :returns: Cadena de caracteres que representa el polinomio parametrizado.
   :rtype: str

   :raises ValueError: Si `this_poly` no es un array de numpy.

   .. seealso:: :func:`pretty_print_lti`, :func:`pretty_print_bicuad_omegayq`

   .. rubric:: Examples

   >>> import numpy as np
   >>> from pytc2.sistemas_lineales import _build_omegayq_str
   >>> this_quad_poly = np.array([1, 2, 3])
   >>> den = np.array([4, 5, 6])
   >>> omegaq_str = _build_omegayq_str(this_quad_poly, den)
   >>> print(omegaq_str)
   r'$s\,0.08333\,\frac{2}{0.1667}$'


.. py:function:: _complementaryColor(my_hex)

   Returns el color RGB complementario.


   :param my_hex: Código hexadecimal del color.
   :type my_hex: str

   :returns: Código hexadecimal del color complementario.
   :rtype: str

   :raises ValueError: Si `my_hex` no es una cadena de caracteres válida o no tiene la longitud correcta.

   .. seealso:: :func:`pzmap`

   .. rubric:: Examples

   >>> from pytc2.sistemas_lineales import _complementaryColor
   >>> _complementaryColor('FFFFFF')
   '000000'


