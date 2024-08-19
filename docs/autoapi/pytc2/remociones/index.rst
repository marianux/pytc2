pytc2.remociones
================

.. py:module:: pytc2.remociones

.. autoapi-nested-parse::

   Created on Thu Mar  2 14:12:53 2023

   @author: mariano



Attributes
----------

.. autoapisummary::

   pytc2.remociones.sig
   pytc2.remociones.sig_pos


Functions
---------

.. autoapisummary::

   pytc2.remociones.tanque_z
   pytc2.remociones.tanque_y
   pytc2.remociones.trim_poly_s
   pytc2.remociones.trim_func_s
   pytc2.remociones.modsq2mod_s
   pytc2.remociones.isFRP
   pytc2.remociones.remover_polo_sigma
   pytc2.remociones.remover_polo_jw
   pytc2.remociones.remover_polo_dc
   pytc2.remociones.remover_polo_infinito
   pytc2.remociones.remover_valor_en_infinito
   pytc2.remociones.remover_valor_en_dc


Module Contents
---------------

.. py:data:: sig

   versión simbólica de sigma, parte real de la variable compleja de Laplace
   s = σ + j.ω
   En caso de necesitar usarla, importar el símbolo desde este módulo.

.. py:data:: sig_pos

   versión simbólica de sigma, parte real positiva de la variable compleja
   de Laplace s = σ + j.ω
   En caso de necesitar usarla, importar el símbolo desde este módulo.

.. py:function:: tanque_z(doska, omegasq)

   Calcula los valores de L y C que componen un tanque resonante LC
   (tanque Z), a partir del valor del residuo (:math:`2.k`) y la omega al cuadrado
   (:math:`\omega^2`) de la expresión de impedancia dada por:

   .. math:: Z_{LC} = \frac{2.k_i.s}{(s^2+\omega^2_i)} = \frac{1}{(s.\frac{1}{2.k_i} + \frac{1}{s \frac{2.k_i}{\omega^2_i} })}

   .. math:: C = \frac{1}{2.k_i}

   .. math:: L = \frac{2.k_i}{\omega^2_i}


   :param doska: Dos veces el residuo.
   :type doska: Symbolic
   :param omegasq: Cuadrado de la omega a la que el tanque resuena.
   :type omegasq: Symbolic

   :returns: * **L** (*Symbolic*) -- Valor del inductor
             * **C** (*Symbolic*) -- Valor del capacitor

   :raises ValueError: Si doska u omegasq no son una instancia de sympy.Expr.

   .. seealso:: :func:`tanque_y`, :func:`trim_func_s`, :func:`isFRP`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import a_equal_b_latex_s, print_latex
   >>> from pytc2.remociones import tanque_z
   >>> k, o = sp.symbols('k, o')
   >>> # Sea la siguiente función de excitación
   >>> L, C = tanque_z( k, o )
   >>> print_latex(a_equal_b_latex_s(sp.symbols('L'), L))
   [LaTex formated equation] '$L=\frac{k}{o}$'
   >>> print_latex(a_equal_b_latex_s(sp.symbols('C'), C))
   [LaTex formated equation] '$C=\frac{1}{k}$'


.. py:function:: tanque_y(doska, omegasq)

   Calcula los valores de L y C que componen un tanque resonante LC
   (tanque Y), a partir del valor del residuo (:math:`2.k`) y la omega al cuadrado
   (:math:`\omega^2`) de la expresión de admitancia dada por:

   .. math:: Y_{LC} = \frac{2.k_i.s}{(s^2+\omega^2_i)} = \frac{1}{(s.\frac{1}{2.k_i} + \frac{1}{s \frac{2.k_i}{\omega^2_i} })}

   .. math:: L = \frac{1}{2.k_i}

   .. math:: C = \frac{2.k_i}{\omega^2_i}


   :param doska: Dos veces el residuo.
   :type doska: Symbolic
   :param omegasq: Cuadrado de la omega a la que el tanque resuena.
   :type omegasq: Symbolic

   :returns: * **L** (*Symbolic*) -- Valor del inductor
             * **C** (*Symbolic*) -- Valor del capacitor

   :raises ValueError: Si doska u omegasq no son una instancia de sympy.Expr.

   .. seealso:: :func:`tanque_z`, :func:`trim_func_s`, :func:`isFRP`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import a_equal_b_latex_s, print_latex
   >>> from pytc2.remociones import tanque_y
   >>> k, o = sp.symbols('k, o')
   >>> # Sea la siguiente función de excitación
   >>> L, C = tanque_y( k, o )
   >>> print_latex(a_equal_b_latex_s(sp.symbols('L'), L))
   [LaTex formated equation] '$C=\frac{1}{k}$'
   >>> print_latex(a_equal_b_latex_s(sp.symbols('C'), C))
   [LaTex formated equation] '$L=\frac{k}{o}$'


.. py:function:: trim_poly_s(this_poly, tol=10**(-6))

   Descarta los coeficientes de un polinomio *this_poly* cuyos valores estén por debajo de
   *tol*.


   :param this_poly: Expresión simbólica del polinomio a ajustar.
   :type this_poly: Symbolic polynomial
   :param tol: Mínimo valor permitido para un coeficiente.
   :type tol: float

   :returns: **poly_acc** -- Polinomio ajustado.
   :rtype: Symbolic

   :raises ValueError: Si this_poly no es una instancia de sympy.Expr polinomial.
       Si tol no es un flotante.

   .. seealso:: :func:`trim_func_s`, :func:`modsq2mod_s`, :func:`isFRP`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s
   >>> from pytc2.remociones import trim_poly_s
   >>> this_poly = sp.poly( 1e-10*s**3 + 2*s**2 + s + 1 , s)
   >>> trim_poly = trim_poly_s( this_poly )
   >>> print(trim_poly)
   2.0*s**2 + 1.0*s + 1.0


.. py:function:: trim_func_s(rat_func, tol=10**(-6))

   Descarta los coeficientes de una función racional *rat_func* cuyos valores estén por debajo de
   *tol*.


   :param rat_func: Expresión simbólica de la función racional a ajustar.
   :type rat_func: Symbolic expresion
   :param tol: Mínimo valor permitido para un coeficiente.
   :type tol: float

   :returns: **trim_func** -- Función racional ajustada.
   :rtype: Symbolic

   :raises ValueError: Si rat_func no es una instancia de sympy.Expr.
       Si tol no es un flotante.

   .. seealso:: :func:`trim_poly_s`, :func:`isFRP`, :func:`trim_poly_s`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s
   >>> from pytc2.remociones import trim_func_s
   >>> rat_func = ( 1e-10*s**3 + 2*s**2 + s + 1)/( 4.3e-10*s**2 + 2*s + 5)
   >>> trim_func = trim_func_s( rat_func )
   >>> print(trim_func)
   (2.0*s**2 + 1.0*s + 1.0)/(2.0*s + 5.0)


.. py:function:: modsq2mod_s(this_func)

   Esta función halla una función de variable compleja T(s), cuyo módulo se
   expresa como la factorización:

   .. math:: \vert T(j\omega) \vert^2 = T(j\omega).T(-j\omega)

   .. math:: T(s) = T(j\omega)\Big\vert_{\omega = s/j}

   Es decir que de todas la singularidades presentes en :math:`\vert T(j\omega) \vert^2`,
   el factor :math:`T(s)` sólo contendrá aquellas que se encuentren en el semiplano izquierdo.

   :param this_func: Expresión simbólica de la función :math:`\vert T(j\omega) \vert^2` a factorizar.
   :type this_func: Symbolic expresion

   :returns: **trim_func** -- Función :math:`T(s)` factorizada.
   :rtype: Symbolic

   :raises ValueError: Si this_func no es una instancia de sympy.Expr.

   .. seealso:: :func:`isFRP`, :func:`trim_func_s`, :func:`trim_poly_s`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s
   >>> from pytc2.remociones import modsq2mod_s
   >>> this_func = (  s**4 + 6*s**2 + 9)/( s**4 - 2*s**2 + 1)
   >>> factor_func = modsq2mod_s( this_func )
   >>> print(factor_func)
   (s**2 + 3)/(s**2 + 2*s + 1)


.. py:function:: isFRP(Imm)

   Chequear si la expresión simbólica Imm es una Función Real y Positiva (FRP).


   :param Imm: La inmitancia a chequear si es FRP.
   :type Imm: symbolic rational function

   :returns: **isFRP** -- A boolean with TRUE value if ff is FRP.
   :rtype: boolean

   :raises ValueError: Si this_func no es una instancia de sympy.Expr.

   .. seealso:: :func:`remover_polo_dc`, :func:`remover_polo_infinito`, :func:`remover_polo_jw`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s
   >>> from pytc2.remociones import isFRP
   >>> Imm = (s**2 + 4*s + 3)/(s**2 + 2*s)
   >>> print(isFRP( Imm ))
   True
   >>> Imm = (s**2 - 4*s + 3)/(s**2 - 2*s)
   >>> print(isFRP( Imm ))
   False


.. py:function:: remover_polo_sigma(imm, sigma, isImpedance=True, isRC=True, sigma_zero=None)

   Se removerá el residuo en sobre el eje :math:`\sigma` (sigma) de la impedancia
   o admitancia (imm) de forma completa o parcial.
   Como resultado de la remoción total, quedará otra función racional definida
   como:

   .. math:: Z_{R} = Z - \frac{k_i}{s + \sigma_i}

   siendo

   .. math:: k_i = \lim\limits _{s\to -\sigma_i} Z (s + \sigma_i)

   Cabe destacar que :math:`Z_{R}` ya no tiene un polo en :math:`\sigma_i`.

   Sin embargo, en cuanto se especifique :math:`\sigma_z`, la remoción parcial
   estará definida como:

   .. math:: Z_{R}\biggr\rfloor_{s=-\sigma_z}= 0 = Z - \frac{k_i}{s + \sigma_i}\biggr\rfloor_{s=-\sigma_z}

   siendo

   .. math:: k_i = Z.(s + \sigma_i)\biggr\rfloor_{s=-\sigma_z}

   Cabe destacar que, para la remoción parcial, :math:`Z_{R}` tendra un cero en
   :math:`\sigma_z` y un polo en :math:`\sigma_i`.


   :param imm: Inmitancia o función que se utilizará para la remoción. Es una función racional
               simbólica que tendrá un polo de orden 1 en :math:`\sigma_i`.
   :type imm: Symbolic
   :param sigma: Frecuencia :math:`\sigma_i` a la que la inmitancia deberá tener un polo.
   :type sigma: float
   :param isImpedance: Booleano que indica si la función imm es una impedancia o admitancia.
   :type isImpedance: bool
   :param isRC: Booleano que indica si la función imm es RC o RL.
   :type isRC: bool
   :param sigma_zero: Frecuencia :math:`\sigma_z` a la que la inmitancia tendrá un cero luego
                      de la remoción.
   :type sigma_zero: float

   :returns: * **imit_r** (*Symbolic*) -- Imitancia luego de la remoción
             * **kk** (*Symbolic*) -- Expresión completa del término removido :math:`\frac{k_i}{s + \sigma_i}`.
             * **R** (*Symbolic*) -- Valor del componente resistivo en la remoción.
             * **CoL** (*Symbolic*) -- Valor del componente capacitivo o inductivo en la remoción.

   :raises ValueError: Si Imm no es una instancia de sympy.Expr.
       Si sigma o sigma_zero no son flotantes.
       Si isImpedance o isRC no son booleanos.

   .. seealso:: :func:`remover_polo_dc`, :func:`remover_polo_infinito`, :func:`remover_polo_jw`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s, a_equal_b_latex_s, print_latex
   >>> from pytc2.remociones import remover_polo_sigma
   >>> # Sea la siguiente función de excitación
   >>> ZZ = (s**2 + 13*s + 32)/(2*(s+1)*(s+6))
   >>> # removemos R1-C1
   >>> sigma_R1C1 = -1
   >>> Z4, ZR1C1, R1, C1 = remover_polo_sigma(ZZ, sigma = sigma_R1C1, isImpedance = True, isRC = True )
   >>> print_latex(a_equal_b_latex_s('Z_3', ZR1C1))
   '$Z_3=\frac{2}{s + 1}$'
   >>> print_latex(a_equal_b_latex_s('Z_4', Z4))
   '$Z_4=\frac{s + 8}{2 \left(s + 6\right)}$'


.. py:function:: remover_polo_jw(imit, omega=None, isImpedance=True, omega_zero=None)

   Se removerá el residuo en sobre el eje :math:`j.\omega` (jota-omega) de la
   impedancia o admitancia (imm) de forma completa o parcial.
   Como resultado de la remoción total, quedará otra función racional definida
   como:

   .. math:: I_{R}=I-\frac{2.k.s}{s^{2}+\omega^{2}}

   siendo

   .. math:: 2.k=\lim\limits _{s^2\to-\omega^2}I\frac{s^{2}+\omega^{2}}{s}

   Cabe destacar que :math:`I_{R}` ya no tendrá sendos polos complejos conjugados en en :math:`\pm\omega`.

   Sin embargo, en cuanto se especifique :math:`\omega_z`, la remoción parcial
   estará definida como:

   .. math:: I_{R}\biggr\rfloor_{s^{2}=-\omega_{z}^{2}}=0=I-\frac{2.k^{'}.s}{s^{2}+\omega^{2}}\biggr\rfloor_{s^{2}=-\omega_{z}^{2}}

   siendo

   .. math:: 2.k^{'}=I.\frac{s^{2}+\omega^{2}}{s}\biggr\rfloor_{s^{2}=-\omega_z^{2}}

   Cabe destacar que, para la remoción parcial, :math:`I_{R}` tendra sendos ceros en
   :math:`\pm j.\omega_z` y sendos polos en :math:`\pm j.\omega`.


   :param imit: Inmitancia o función que se utilizará para la remoción. Es una función racional
                simbólica que tendrá un polo de orden 1 en :math:`j\omega`.
   :type imit: Symbolic
   :param omega: Frecuencia :math:`\sigma_i` a la que la inmitancia deberá tener un polo.
   :type omega: float
   :param isImpedance: Booleano que indica si la función imm es una impedancia o admitancia.
   :type isImpedance: bool
   :param omega_zero: Frecuencia :math:`\sigma_z` a la que la inmitancia tendrá un cero luego
                      de la remoción.
   :type omega_zero: float

   :returns: * **imit_r** (*Symbolic*) -- Imitancia luego de la remoción
             * **kk** (*Symbolic*) -- Expresión completa del término removido :math:`\frac{2.k.s}{s^{2}+\omega^{2}}`.
             * **R** (*Symbolic*) -- Valor del componente resistivo en la remoción.
             * **CoL** (*Symbolic*) -- Valor del componente capacitivo o inductivo en la remoción.

   :raises ValueError: Si Imm no es una instancia de sympy.Expr.
       Si sigma o sigma_zero no son flotantes.
       Si isImpedance o isRC no son booleanos.

   .. seealso:: :func:`remover_polo_dc`, :func:`remover_polo_infinito`, :func:`remover_polo_sigma`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s, a_equal_b_latex_s, print_latex
   >>> from pytc2.remociones import remover_polo_jw
   >>> # Sea la siguiente función de excitación
   >>> YY = (s * (3*s**2+7) )/((s**2+1)*(s**2+3))
   >>> # removemos R1-C1
   >>> omega_L2C2 = 1
   >>> Y4, Yt2, L2, C2 = remover_polo_jw(YY, isImpedance = False, omega = omega_L2C2 )
   >>> print_latex(a_equal_b_latex_s('Y_3(s)', Yt2))
   '$Y_3(s)=\frac{2 s}{s^{2} + 1}$'
   >>> print_latex(a_equal_b_latex_s('Y_4(s)', Y4))
   '$Y_4(s)=\frac{s}{s^{2} + 3}$'


.. py:function:: remover_polo_dc(imit, omega_zero=None, isSigma=False)

   Se removerá el residuo en continua (:math:`j.0`) de la
   impedancia o admitancia (inmitancia o imit) de forma completa o parcial.
   Como resultado de la remoción total, quedará otra función racional definida
   como:

   .. math:: I_{R}=I-\frac{k_0}{s}

   siendo

   .. math:: k_0=\lim\limits _{s\to0}I.s

   Cabe destacar que :math:`I_{R}` ya no tendrá polo en :math:`j.0`.

   Sin embargo, en cuanto se especifique :math:`\omega_z`, la remoción parcial
   estará definida como:

   .. math:: I_{R}\biggr\rfloor_{s^{2}=-\omega_z^{2}}=0=I-\frac{k_{0}^{'}}{s}\biggr\rfloor_{s^{2}=-\omega_z^{2}}

   siendo

   .. math:: k_{0}^{'}=I.s\biggr\rfloor_{s^{2}=-\omega_z^{2}}

   Cabe destacar que, para la remoción parcial, :math:`I_{R}` tendra sendos ceros en
   :math:`\pm j.\omega_z` y un polo en :math:`j.0`.


   :param imit: Inmitancia o función que se utilizará para la remoción. Es una función racional
                simbólica que tendrá un polo de orden 1 en :math:`j\omega`.
   :type imit: Symbolic
   :param isSigma: Booleano que indica si la función imm es una impedancia o admitancia.
   :type isSigma: bool
   :param omega_zero: Frecuencia :math:`\sigma_z` a la que la inmitancia tendrá un cero luego
                      de la remoción.
   :type omega_zero: float

   :returns: * **imit_r** (*Symbolic*) -- Imitancia luego de la remoción
             * **k_cero** (*Symbolic*) -- Expresión completa del término removido :math:`\frac{2.k.s}{s^{2}+\omega^{2}}`.

   :raises ValueError: Si imit no es una instancia de sympy.Expr.
       Si omega_zero no es flotante.
       Si isSigma o isRC no son booleanos.

   .. seealso:: :func:`remover_polo_jw`, :func:`remover_polo_infinito`, :func:`remover_polo_sigma`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s, a_equal_b_latex_s, print_latex
   >>> from pytc2.remociones import remover_polo_dc
   >>> # Sea la siguiente función de excitación
   >>> YY = 3*s*(s**2+sp.Rational(7,3))/(s**2+2)/(s**2+5)
   >>> omega_L2C2 = 1
   >>> Z2, Zc1 = remover_polo_dc(1/YY, omega_zero = omega_L2C2 )
   >>> # Zc1 es la admitancia removida
   >>> # extraigo C1
   >>> C1 = 1/(s*Zc1)
   >>> print_latex(a_equal_b_latex_s('Z_1(s)', Zc1))
   $Z_1(s)=\frac{1}{s}$
   >>> print_latex(a_equal_b_latex_s('Z_2(s)', Z2))
   $Z_2(s)=\frac{\left(s^{2} + 1\right) \left(s^{2} + 3\right)}{s \left(3 s^{2} + 7\right)}$'


.. py:function:: remover_polo_infinito(imit, omega_zero=None, isSigma=False)

   Se removerá el residuo en infinito  de la impedancia o admitancia (inmitancia
   o imit) de forma completa o parcial. Como resultado de la remoción total,
   quedará otra función racional definida como:

   .. math:: I_R = I - s.k_\infty

   siendo

   .. math:: k_{\infty}=\lim\limits _{s\to\infty}I.\frac{1}{s}

   Cabe destacar que :math:`I_{R}` ya no tendrá polo en :math:`j.\infty`.

   En cuanto se especifique :math:`\omega_z`, la remoción parcial estará definida
   como:

   .. math:: I_{R}\biggr\rfloor_{s^{2}=-\omega_z^{2}}=0=I-s.k_{\infty}^{'}\biggr\rfloor_{s^{2}=-\omega_z^{2}}

   siendo

   .. math:: k_{\infty}^{'}=I.\frac{1}{s}\biggr\rfloor_{s^{2}=-\omega_z^{2}}

   Cabe destacar que, para la remoción parcial, :math:`I_{R}` tendra sendos ceros en
   :math:`\pm j.\omega_z` y un polo en :math:`j.\infty`. Lo anterior se cumple
   siempre que isSigma = False, de lo contrario

   .. math:: I_{R}\biggr\rfloor_{s=-\omega_z}=0=I-s.k_{\infty}^{'}\biggr\rfloor_{s=-\omega_z}

   siendo

   .. math:: k_{\infty}^{'}=I.\frac{1}{s}\biggr\rfloor_{s=-\omega_z}

   Al igual que antes, destacar que para la remoción parcial, :math:`I_{R}` tendrá
   un cero en :math:`-\sigma_z = \omega_z` y un polo en :math:`j.\infty`.


   :param imit: Inmitancia o función que se utilizará para la remoción. Es una función racional
                simbólica que tendrá un polo de orden 1 en :math:`j\omega`.
   :type imit: Symbolic
   :param isSigma: Booleano que indica si la función imm es una impedancia o admitancia.
   :type isSigma: bool
   :param omega_zero: Frecuencia :math:`\sigma_z` a la que la inmitancia tendrá un cero luego
                      de la remoción.
   :type omega_zero: float

   :returns: * **imit_r** (*Symbolic*) -- Imitancia luego de la remoción
             * **k_inf** (*Symbolic*) -- Expresión completa del término removido :math:`s.k_{\infty}`.

   :raises ValueError: Si imit no es una instancia de sympy.Expr.
       Si omega_zero no es flotante.
       Si isSigma o isRC no son booleanos.

   .. seealso:: :func:`remover_polo_dc`, :func:`remover_polo_jw`, :func:`remover_polo_sigma`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s, a_equal_b_latex_s, print_latex
   >>> from pytc2.remociones import remover_polo_infinito
   >>> # Sea la siguiente función de excitación
   >>> YY = 3*s*(s**2+sp.Rational(7,3))/(s**2+2)/(s**2+5)
   >>> Z2, Z1 = remover_polo_infinito(1/YY)
   >>> # Z1 es la admitancia removida
   >>> # extraigo L1
   >>> L1 = Z1/s
   >>> print_latex(a_equal_b_latex_s('Z_1(s)', Z1))
   '$Z_1(s)=\frac{s}{3}$'
   >>> print_latex(a_equal_b_latex_s('Z_2(s)', Z2))
   '$Z_2(s)=\frac{2 \cdot \left(7 s^{2} + 15\right)}{3 s \left(3 s^{2} + 7\right)}$'


.. py:function:: remover_valor_en_infinito(imit, sigma_zero=None)

   Se removerá un valor real de la impedancia o admitancia (inmitancia
   o imit) de forma completa o parcial. Como resultado de la remoción total,
   quedará otra función racional definida como:

   .. math:: I_R = I - k_\infty

   siendo

   .. math:: k_{\infty}=\lim\limits _{s\to\infty}I

   En cuanto se especifique :math:`\sigma_z`, la remoción parcial estará definida
   como:

   .. math:: I_{R}\biggr\rfloor_{s=-\sigma_z}=0=I-k_{\infty}^{'}\biggr\rfloor_{s=-\sigma_z}

   siendo

   .. math:: k_{\infty}^{'}=I\biggr\rfloor_{s=-\sigma_z}

   Cabe destacar que, para la remoción parcial, :math:`I_{R}` tendra un cero en
   :math:`-\sigma_z` y un valor real en :math:`\infty`.

   :param imit: Inmitancia o función que se utilizará para la remoción. Es una función racional
                simbólica que tendrá un polo de orden 1 en :math:`j\omega`.
   :type imit: Symbolic
   :param sigma_zero: Frecuencia :math:`\sigma_z` a la que la inmitancia tendrá un cero luego
                      de la remoción.
   :type sigma_zero: float

   :returns: * **imit_r** (*Symbolic*) -- Imitancia luego de la remoción
             * **k_inf** (*Symbolic*) -- Expresión completa del término removido :math:`s.k_{\infty}`.

   :raises ValueError: Si imit no es una instancia de sympy.Expr.
       Si sigma_zero no es flotante.

   .. seealso:: :func:`remover_valor_en_dc`, :func:`remover_polo_en_infinito`, :func:`remover_polo_en_dc`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s, a_equal_b_latex_s, print_latex
   >>> from pytc2.remociones import remover_valor_en_infinito
   >>> # Sea la siguiente función de excitación
   >>> ZZ = (s**2 + 13*s + 32)/(3*s**2 + 27*s+ 44)
   >>> Z2, Z1 = remover_valor_en_infinito(ZZ)
   >>> print_latex(a_equal_b_latex_s('Z_1(s)', Z1))
   '$Z_1(s)=\frac{1}{3}$'
   >>> print_latex(a_equal_b_latex_s('Z_2(s)', Z2))
   '$Z_2(s)=\frac{4 \cdot \left(3 s + 13\right)}{3 \cdot \left(3 s^{2} + 27 s + 44\right)}$'


.. py:function:: remover_valor_en_dc(imit, sigma_zero=None)

   Se removerá un valor constante en continua (s=0) de la imitancia (imit) de forma
   completa. Como resultado de la remoción, quedará otra función racional definida
   como:

   .. math:: I_R = I - k_0

   siendo

   .. math:: k_0 = \lim\limits _{s \to 0}I

   En cuanto se especifique :math:`\sigma_z`, la remoción parcial estará definida
   como:

   .. math:: I_{R}\biggr\rfloor_{s=-\sigma_z}=0=I-k_{0}^{'}\biggr\rfloor_{s=-\sigma_z}

   siendo

   .. math:: k_{0}^{'}=I\biggr\rfloor_{s=-\sigma_z}

   Cabe destacar que, para la remoción parcial, :math:`I_{R}` tendra un cero en
   :math:`-\sigma_z` y un valor real en 0.

   :param imit: Inmitancia o función que se utilizará para la remoción. Es una función racional
                simbólica que tendrá un polo de orden 1 en :math:`j\omega`.
   :type imit: Symbolic
   :param sigma_zero: Frecuencia :math:`\sigma_z` a la que la inmitancia tendrá un cero luego
                      de la remoción.
   :type sigma_zero: float

   :returns: * **imit_r** (*Symbolic*) -- Imitancia luego de la remoción
             * **k_0** (*Symbolic*) -- Expresión completa del término removido :math:`k_0`.

   :raises ValueError: Si imit no es una instancia de sympy.Expr.
       Si sigma_zero no es flotante.

   .. seealso:: :func:`remover_valor_en_infinito`, :func:`remover_polo_en_infinito`, :func:`remover_polo_en_dc`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.general import s, a_equal_b_latex_s, print_latex
   >>> from pytc2.remociones import remover_valor_en_dc
   >>> # Sea la siguiente función de excitación
   >>> ZZ = (s**2 + 13*s + 32)/(3*s**2 + 27*s+ 44)
   >>> Z2, Z1 = remover_valor_en_dc(1/ZZ)
   >>> print_latex(a_equal_b_latex_s('Z_1(s)', Z1))
   :math:`$Z_1(s)=\frac{11}{8}$`
   >>> print_latex(a_equal_b_latex_s('Z_2(s)', Z2))
   :math:`$Z_2(s)=\frac{s \left(13 s + 73\right)}{8 \left(s^{2} + 13 s + 32\right)}$`


