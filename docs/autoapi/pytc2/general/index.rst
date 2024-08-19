pytc2.general
=============

.. py:module:: pytc2.general

.. autoapi-nested-parse::

   Created on Thu Mar  2 11:22:31 2023

   Originally based on the work of Combination of 2011 Christopher Felton
   Further modifications were added for didactic purposes
   by Mariano Llamedo llamedom _at_ frba_utn_edu_ar

   @author: marianux



Attributes
----------

.. autoapisummary::

   pytc2.general.pytc2_full_path
   pytc2.general.small_val
   pytc2.general.s
   pytc2.general.w


Functions
---------

.. autoapisummary::

   pytc2.general.get_home_directory
   pytc2.general.pp
   pytc2.general.factorSOS
   pytc2.general.symbfunc2tf
   pytc2.general.simplify_n_monic
   pytc2.general.Chebyshev_polynomials
   pytc2.general.a_equal_b_latex_s
   pytc2.general.expr_simb_expr
   pytc2.general.to_latex
   pytc2.general.print_latex
   pytc2.general.print_console_alert
   pytc2.general.print_console_subtitle
   pytc2.general.print_subtitle
   pytc2.general.db2nepper
   pytc2.general.nepper2db


Module Contents
---------------

.. py:data:: pytc2_full_path

   Path a donde se encuentra pyTC2 localmente.

.. py:data:: small_val

   Es un valor muy pequeño para que las funciones que tienen restringido su evaluación
   en 0 no arrojen warnings ni errores. e.g. los cálculos de los logaritmos

.. py:data:: s

   Variable compleja de Laplace s = σ + j.ω
   En caso de necesitar usarla, importar el símbolo desde este módulo.

.. py:data:: w

   Fourier real variable ω
   En caso de necesitar usarla, importar el símbolo desde este módulo.

.. py:function:: get_home_directory()

.. py:function:: pp(z1, z2)

   Asocia en paralelo dos impedancias o en serie dos admitancias.


   :param z1: Inmitancia 1.
   :type z1: Symbolic o float
   :param z2: Inmitancia 2.
   :type z2: Symbolic o float

   :returns: **zp** -- Inmitancia resultante.
   :rtype: Symbolic o float

   :raises ValueError: Si alguno de los argumentos no es de tipo `Symbolic`.:

   .. seealso:: :func:`print_latex`, :func:`to_latex`, :func:`a_equal_b_latex_s`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import pp
   >>> # Asociación en paralelo de dos impedancias
   >>> z1 = sp.symbols('z1')
   >>> z2 = sp.symbols('z2')
   >>> zp = pp(z1, z2)
   >>> print(zp)
   z1*z2/(z1 + z2)
   >>> # Asociación en serie de dos admitancias
   >>> y1 = 1/z1
   >>> y2 = 1/z2
   >>> yp = pp(y1, y2)
   >>> print(yp)
   1/(z1*z2*(1/z2 + 1/z1))


.. py:function:: factorSOS(ratfunc, decimals=4)

   Factoriza una función racional simbólica, en polinomios de segundo y primer
   orden.


   :param ratfunc: Función racional simbólica.
   :type ratfunc: Expr. simbólica
   :param decimals: Cantidad de decimales para la evaluación simbólica.
   :type decimals: entero

   :returns: Función racional simbólica factorizada.
   :rtype: Expr. simbólica

   :raises ValueError: Si la entrada no es una expresión simbólica.

   .. seealso:: :func:`symbfunc2tf`, :func:`simplify_n_monic`, :func:`a_equal_b_latex_s`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s, factorSOS
   >>> tt = (s**4 + 8*s**3 + 18*s**2 + 11*s + 2)/(s**3 + 16*s**2 + 65*s + 14)
   >>> factorized_tt, _, _ = factorSOS(tt)
   >>> print(factorized_tt)
   (s + 0.382)*(s + 0.438)*(s + 2.62)*(s + 4.56)/((s + 0.228)*(s + 7.0)*(s + 8.77))


.. py:function:: symbfunc2tf(tt)

   Convierte una función racional simbólica, con coeficientes numéricos
   (convertibles a flotante), en un objeto transfer function.


   :param tt: Función racional simbólica.
   :type tt: Expr. simbólica

   :returns: TransferFunction que representa numéricamente la función.
   :rtype: TransferFunction

   :raises ValueError: Si la entrada no es una expresión simbólica.

   .. seealso:: :func:`simplify_n_monic`, :func:`to_latex`, :func:`factorSOS`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s, symbfunc2tf
   >>> tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
   >>> simplified_tt = symbfunc2tf(tt)
   >>> print(simplified_tt)
   TransferFunctionContinuous(
   array([0.5, 1.5, 1. ]),
   array([1. , 2.5, 1.5]),
   dt: None
   )


.. py:function:: simplify_n_monic(tt)

   Simplifica un polinomio de fracciones en forma mónica.


   :param tt: Polinomio de fracciones a simplificar.
   :type tt: Expr

   :returns: Polinomio simplificado en forma monica.
   :rtype: Expr

   :raises ValueError: Si la entrada no es una expresión simbólica.

   .. seealso:: :func:`print_latex`, :func:`to_latex`, :func:`a_equal_b_latex_s`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s, simplify_n_monic
   >>> tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
   >>> simplified_tt = simplify_n_monic(tt)
   >>> print(simplified_tt)
   (s + 2)/(2*s + 3)


.. py:function:: Chebyshev_polynomials(nn)

   Calcula el polinomio de Chebyshev de grado nn.


   :param nn: Grado del polinomio de Chebyshev.
   :type nn: int

   :returns: **Ts** -- Matriz de parámetros de transferencia scattering.
   :rtype: Symbolic Matrix

   :raises ValueError: Si nn no es un entero positivo.

   .. seealso:: :func:`print_latex`, :func:`to_latex`, :func:`a_equal_b_latex_s`

   .. rubric:: Examples

   >>> from pytc2.general import Chebyshev_polynomials
   >>> Ts = Chebyshev_polynomials(3)
   >>> print(Ts)
   w*(4*w**2 - 3)


.. py:function:: a_equal_b_latex_s(a, b)

   A partir de un string o expresión de SymPy (a), y otra expresión de SymPy (b):

   .. math:: a = b

   en un nuevo string formateado para visualizarse en LaTeX.


   :param a: Símbolo o cadena para el lado izquierdo de la igualdad.
   :type a: Symbolic or str
   :param b: Símbolo o cadena para el lado derecho de la igualdad.
   :type b: Symbolic, str o lista de ambas

   :returns: **str** -- String formateado en LaTeX representando la igualdad.
   :rtype: string

   :raises ValueError: Si a no es un símbolo ni una cadena.
       Si b no es un símbolo.

   .. seealso:: :func:`expr_simb_expr`, :func:`print_latex`, :func:`to_latex`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import a_equal_b_latex_s, print_latex
   >>> s = sp.symbols('s')
   >>> tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
   >>> print(a_equal_b_latex_s(sp.symbols('tt'), tt))
   'tt=\frac{s^{2} + 3 s + 2}{2 s^{2} + 5 s + 3}$'
   >>> print_latex(a_equal_b_latex_s(sp.symbols('tt'), tt))
   [LaTex formated equation]


.. py:function:: expr_simb_expr(a, b, symbol='=')

   A partir de un string o expresión de SymPy (a), y otra expresión de SymPy (b):

   a symbol b

   en un nuevo string formateado para visualizarse en LaTeX.


   :param a: Símbolo o cadena para el lado izquierdo de la expresión.
   :type a: Symbolic or str
   :param b: Símbolo o cadena para el lado derecho de la expresión.
   :type b: Symbolic or str
   :param symbol: Símbolo de operación entre a y b (por defecto es '=').
   :type symbol: str, optional

   :returns: String formateado en LaTeX representando la expresión.
   :rtype: str

   :raises ValueError: Si a no es un símbolo ni una cadena.
       Si b no es un símbolo.

   .. seealso:: :func:`a_equal_b_latex_s`, :func:`print_latex`, :func:`to_latex`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import expr_simb_expr, print_latex
   >>> s = sp.symbols('s')
   >>> tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
   >>> tt1 = (s**2 + 4*s + 7) / (2*s**2 + 5*s + 3)
   >>> print_latex(expr_simb_expr('tt', tt1, '\neq'))
   [LaTex formated equation]
   >>> print_latex(expr_simb_expr('tt', tt))
   [LaTex formated equation]


.. py:function:: to_latex(unsimbolo)

   Convierte un símbolo en un string formateado para visualizarse en LaTeX.


   :param unsimbolo: Símbolo o cadena a convertir a formato LaTeX.
   :type unsimbolo: Symbolic or str

   :returns: String formateado en LaTeX.
   :rtype: str

   :raises ValueError: Si unsimbolo no es un símbolo ni una cadena.

   .. seealso:: :func:`print_latex`, :func:`to_latex`, :func:`to_latex`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import to_latex, print_latex
   >>> print(to_latex(sp.symbols('x')))
   $x$
   >>> print_latex(to_latex(sp.symbols('x')))
   [LaTex formated equation]


.. py:function:: print_latex(unstr)

   Muestra una expresión LaTeX en formato matemático.


   :param unstr: Cadena que representa la expresión LaTeX.
   :type unstr: str

   :returns: Esta función no devuelve nada, simplemente muestra la expresión en formato LaTeX.
   :rtype: None

   :raises ValueError: Si unstr no es una cadena.

   .. seealso:: :func:`print_subtitle`, :func:`to_latex`, :func:`to_latex`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import to_latex, print_latex
   >>> print(to_latex('x'))
   $x$
   >>> print_latex(to_latex('x'))
   [LaTex formated equation]


.. py:function:: print_console_alert(unstr)

   Imprime una cadena rodeada por símbolos de alerta en la consola.


   :param unstr: Cadena a imprimir.
   :type unstr: str

   :returns: Esta función no devuelve nada, simplemente imprime la cadena en la consola.
   :rtype: None

   :raises ValueError: Si unstr no es una cadena.

   .. seealso:: :func:`print_subtitle`, :func:`print_latex`, :func:`print_alert`

   .. rubric:: Examples

   >>> from pytc2.general import print_console_alert
   >>> print_console_alert('Advertencia: Datos incompletos')
   ##################################
   # Advertencia: Datos incompletos #
   ##################################


.. py:function:: print_console_subtitle(unstr)

   Imprime un subtítulo en la consola.


   :param unstr: Cadena que representa el subtítulo.
   :type unstr: str

   :returns: Esta función no devuelve nada, simplemente imprime el subtítulo en la consola.
   :rtype: None

   :raises ValueError: Si unstr no es una cadena.

   .. seealso:: :func:`print_subtitle`, :func:`print_latex`, :func:`print_console_alert`

   .. rubric:: Examples

   >>> from pytc2.general import print_console_subtitle
   >>> print_console_subtitle('Subtítulo')
   Subtítulo
   ---------


.. py:function:: print_subtitle(unstr)

   Imprime un subtítulo.


   :param unstr: Cadena que representa el subtítulo.
   :type unstr: str

   :returns: Esta función no devuelve nada, simplemente imprime el subtítulo.
   :rtype: None

   :raises ValueError: Si unstr no es una cadena.

   .. seealso:: :func:`print_latex`, :func:`print_console_alert`, :func:`print_console_subtitle`

   .. rubric:: Examples

   >>> from pytc2.general import print_subtitle
   >>> print_subtitle('Subtítulo')
   <IPython.core.display.Markdown object>


.. py:function:: db2nepper(at_en_db)

   Convierte una magnitud en decibels a su equivalente en nepers.


   :param at_en_db: Magnitud en decibelios a convertir.
   :type at_en_db: float or numpy.ndarray

   :returns: Equivalente en nepers.
   :rtype: float or numpy.ndarray

   :raises ValueError: Si at_en_db no es de tipo `float`.:

   .. seealso:: :func:`nepper2db`

   .. rubric:: Examples

   >>> from pytc2.general import db2nepper
   >>> db2nepper(20.)
   2.3025850929940455
   >>> db2nepper(1.)
   0.11512925464970228


.. py:function:: nepper2db(at_en_np)

   Convierte una magnitud en neperios a su equivalente en decibelios.


   :param at_en_np: Magnitud en neperios a convertir.
   :type at_en_np: float or numpy.ndarray

   :returns: Equivalente en decibelios.
   :rtype: float or numpy.ndarray

   :raises ValueError: Si at_en_db no es de tipo `float`.:

   .. seealso:: :func:`db2nepper`

   .. rubric:: Examples

   >>> from pytc2.general import nepper2db
   >>> nepper2db(1.)
   8.685889638065037
   >>> nepper2db(2.3025850929940455)
   20.


