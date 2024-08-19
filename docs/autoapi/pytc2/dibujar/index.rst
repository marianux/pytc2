pytc2.dibujar
=============

.. py:module:: pytc2.dibujar

.. autoapi-nested-parse::

   Created on Thu Mar  2 11:24:17 2023

   @author: mariano



Attributes
----------

.. autoapisummary::

   pytc2.dibujar.elementos_dic
   pytc2.dibujar.elementos_keys
   pytc2.dibujar.elementos_keys_str
   pytc2.dibujar.elementos_keys_str


Functions
---------

.. autoapisummary::

   pytc2.dibujar.dibujar_Tee
   pytc2.dibujar.dibujar_Pi
   pytc2.dibujar.dibujar_lattice
   pytc2.dibujar.dibujar_cauer_RC_RL
   pytc2.dibujar.dibujar_cauer_LC
   pytc2.dibujar.dibujar_foster_derivacion
   pytc2.dibujar.dibujar_foster_serie
   pytc2.dibujar.dibujar_puerto_entrada
   pytc2.dibujar.dibujar_puerto_salida
   pytc2.dibujar.dibujar_espaciador
   pytc2.dibujar.dibujar_funcion_exc_abajo
   pytc2.dibujar.dibujar_funcion_exc_arriba
   pytc2.dibujar.dibujar_elemento_serie
   pytc2.dibujar.dibujar_espacio_derivacion
   pytc2.dibujar.dibujar_cierre
   pytc2.dibujar.dibujar_elemento_derivacion
   pytc2.dibujar.dibujar_tanque_RC_serie
   pytc2.dibujar.dibujar_tanque_RC_derivacion
   pytc2.dibujar.dibujar_tanque_RL_serie
   pytc2.dibujar.dibujar_tanque_RL_derivacion
   pytc2.dibujar.dibujar_tanque_serie
   pytc2.dibujar.dibujar_tanque_derivacion


Module Contents
---------------

.. py:data:: elementos_dic

.. py:data:: elementos_keys

.. py:data:: elementos_keys_str

.. py:data:: elementos_keys_str

.. py:function:: dibujar_Tee(ZZ, return_components=False)

   Dibuja una red Tee a partir de la matriz Z.

   :param ZZ: Matriz de impedancia Z.
   :type ZZ: sympy.Matrix
   :param return_components: Indica si se deben devolver los componentes individuales de la red (Za, Zb, Zc). Por defecto es False.
   :type return_components: bool, optional

   :returns: Si return_components es True, devuelve una lista con los componentes individuales de la red (Za, Zb, Zc).
             Si return_components es False, no devuelve nada.
   :rtype: list or None

   :raises ValueError: Si ZZ no es una instancia de sympy.Matrix.

   .. seealso:: :func:`dibujar_Pi`, :func:`dibujar_lattice`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.dibujar import dibujar_Tee
   >>> dibujar_Tee(sp.Matrix([[5, 2], [2, 6]]))
   [dibujo de la red]

   Ver el `tutorial de cuadripolos elementales <https://pytc2.readthedocs.io/en/latest/cuadripolos_elementales.html>`__ para
   observar el resultado de ésta y otras funciones.



.. py:function:: dibujar_Pi(YY, return_components=False)

   Dibuja una red Pi a partir de la matriz Y.


   :param YY: Matriz de admitancia Y.
   :type YY: Symbolic Matrix
   :param return_components: Indica si se deben devolver los componentes individuales de la red (Ya, Yb, Yc). Por defecto es False.
   :type return_components: bool, optional

   :returns: Si return_components es True, devuelve una lista con los componentes individuales de la red (Ya, Yb, Yc).
             Si return_components es False, no devuelve nada.
   :rtype: None or list

   :raises ValueError: Si YY no es una instancia de sympy.Matrix.

   .. seealso:: :func:`dibujar_Tee`, :func:`dibujar_lattice`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.dibujar import dibujar_Pi
   >>> Ya, Yb, Yc = dibujar_Pi(sp.Matrix([[5, -2], [-2, 6]]), return_components=True)
   [dibujo de la red]

   Ver el `tutorial de cuadripolos elementales <https://pytc2.readthedocs.io/en/latest/cuadripolos_elementales.html>`__ para
   observar el resultado de ésta y otras funciones.



.. py:function:: dibujar_lattice(ZZ, return_components=False)

   Dibuja una red Lattice a partir de una matriz de parámetros Z.


   :param ZZ: Parámetros Z de la red. Si no se proporciona, solo se genera el dibujo. El valor predeterminado es None.
   :type ZZ: Matriz simbólica, opcional
   :param return_components: Indica si se deben devolver los componentes de la red Lattice simétrica (Za y Zb). El valor predeterminado es False.
   :type return_components: bool, opcional

   :returns: Si return_components es True, devuelve una lista con los componentes Za y Zb de la red Lattice simétrica.
             Si return_components es False, devuelve None.
   :rtype: list or None

   :raises ValueError: Si ZZ no es una instancia de sympy.Matrix.
       Si ZZ no es de 2x2

   .. seealso:: :func:`dibujar_Pi`, :func:`dibujar_Tee`

   Ejemplos
   --------
   >>> import sympy as sp
   >>> from pytc2.dibujar import dibujar_lattice
   >>> Za, Zb = dibujar_lattice(sp.Matrix([[5, 2], [2, 6]]), return_components=True)


   Ver el `tutorial de cuadripolos elementales <https://pytc2.readthedocs.io/en/latest/cuadripolos_elementales.html>`__ para
   observar el resultado de ésta y otras funciones.



.. py:function:: dibujar_cauer_RC_RL(ki=None, y_exc=None, z_exc=None)

   Dibuja una red disipativa escalera (RC-RL) a partir de una expansión en
   fracciones continuas (Método de Cauer). Dependiendo se especifique `z_exc`
   o `y_exc` y el tipo de residuos de `ki` se dibujará la red correspondiente.
   En caso que se trate de redes RC, la forma matemática será:

   .. math:: Z_{RC}(s)= \frac{1}{s.C_1} + \frac{1}{ \frac{1}{R_1} + \frac{1}{ \frac{1}{s.C_2} + \cdots } } =
        R_1 + \frac{1}{ s.C_1 + \frac{1}{ R_2 + \cdots } }

   .. math:: Y_{RC}(s)= s.C_1 + \frac{1}{ R_1 + \frac{1}{ s.C_2 + \cdots } } =
        \frac{1}{R_1} + \frac{1}{ s.C_1 + \frac{1}{ \frac{1}{R_2} + \cdots } }

   :param ki: Será una lista que contenga los residuos [k0, ki, koo ] como expresiones
              simbólicas. Esta lista la provee la función :func:`cauer_RC`. El valor
              predeterminado es None. Siendo:

              * k0  : Residuo de la función en DC o :math:`\sigma \to 0`.
              * koo : Residuo de la función en infinito o :math:`\sigma \to \infty`.
              * ki  : Residuo de la función en :math:`\sigma_i` o :math:`\sigma \to -\sigma_i`
   :type ki: lista con expresiones simbólicas

   :rtype: None

   :raises ValueError: Si y_exc y z_exc no son una instancia de sympy.Expr.

   .. seealso:: :func:`cauer_RC`, :func:`foster_zRC2yRC`, :func:`dibujar_cauer_LC`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.sintesis_dipolo import cauer_RC
   >>> from pytc2.dibujar import dibujar_cauer_RC_RL
   >>> s = sp.symbols('s ', complex=True)
   >>> # Sea la siguiente función de excitación
   >>> ZRC = (s**2 + 4*s + 3)/(s**2 + 2*s)
   >>> # Implementaremos FF mediante Cauer 1 o remociones continuas en infinito
   >>> koo, ZRC_cauer_oo, rem = cauer_RC(ZRC, remover_en_inf=True)
   >>> # Tratamos a nuestra función inmitancia como una Z
   >>> dibujar_cauer_RC_RL(koo, z_exc = ZRC_cauer_oo)
   >>> # Tratamos a nuestra función inmitancia como una Y
   >>> dibujar_cauer_RC_RL(koo, y_exc = ZRC_cauer_oo)


.. py:function:: dibujar_cauer_LC(ki=None, y_exc=None, z_exc=None)

   Dibuja una red escalera no disipativa, a partir de la expansión en fracciones
   continuas (Método de Cauer). Dependiendo se especifique `z_exc`
   o `y_exc` y el tipo de residuos de `ki` se dibujará la red correspondiente.
   La forma matemática será:

   .. math:: Z(s)= \frac{1}{s.C_1} + \frac{1}{ \frac{1}{s.L_1} + \frac{1}{ \frac{1}{s.C_2} + \cdots } } =
            s.L_1 + \frac{1}{ s.C_1 + \frac{1}{ s.L_2 + \cdots } }

   .. math:: Y(s)= \frac{1}{s.L_1} + \frac{1}{ \frac{1}{s.C_1} + \frac{1}{ \frac{1}{s.L_2} + \cdots } } =
            s.C_1 + \frac{1}{ s.L_1 + \frac{1}{ s.C_2 + \cdots } }


   :param ki: Será una lista que contenga los residuos [k0, ki, koo ] como expresiones
              simbólicas. Esta lista la provee la función :func:`cauer`.
              El valor predeterminado es None. Siendo:

              * k0  : Residuo de la función en DC o :math:`s \to 0`.
              * koo : Residuo de la función en infinito o :math:`s \to \infty`.
              * ki  : Residuo de la función en :math:`\omega_i` o :math:`s^2 \to -\omega^2_i`
   :type ki: lista con expresiones simbólicas

   :rtype: None

   :raises ValueError: Si y_exc y z_exc no son una instancia de sympy.Expr.

   .. seealso:: :func:`cauer_LC`, :func:`foster_zRC2yRC`, :func:`dibujar_cauer_LC`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.sintesis_dipolo import cauer_LC
   >>> from pytc2.dibujar import dibujar_cauer_LC
   >>> s = sp.symbols('s ', complex=True)
   >>> # Sea la siguiente función de excitación
   >>> FF = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
   >>> # Implementaremos FF mediante Cauer 1 o remociones continuas en infinito
   >>> koo, F_cauer_oo, rem = cauer_LC(FF, remover_en_inf=True)
   >>> # Tratamos a nuestra función inmitancia como una Z
   >>> dibujar_cauer_LC(koo, z_exc = F_cauer_oo)
   >>> # Tratamos a nuestra función inmitancia como una Y
   >>> dibujar_cauer_LC(koo, y_exc = F_cauer_oo)


.. py:function:: dibujar_foster_derivacion(k0=None, koo=None, ki=None, kk=None, y_exc=None)

   Dibuja una red no disipativa a partir de una expansión en fracciones simples
   (Método de Foster). La forma matemática es:

   .. math:: Y(s)= \frac{k_0}{s} + k_\infty.s + \sum_{i=1}^N\frac{2.k_i.s}{s^2+\omega_i^2}

   Esta función provee una interpretación circuital al resultado de la función
   :func:`foster`.


   :param k0: Residuo de la función en DC o :math:`s \to 0`. El valor predeterminado es None.
   :type k0: simbólica, opcional
   :param koo: Residuo de la función en infinito o :math:`s \to \infty`. El valor predeterminado es None.
   :type koo: simbólica, opcional
   :param ki: Residuo de la función en :math:`\omega_i` o :math:`s^2 \to -\omega^2_i`. El valor predeterminado es None.
   :type ki: simbólica, list o tuple opcional
   :param kk: Residuo de la función en :math:`\sigma_i` o :math:`\omega \to -\omega_i`. El valor predeterminado es None.
   :type kk: simbólica, opcional

   :rtype: None

   :raises ValueError: Si cualquiera de los argumentos no son una instancia de sympy.Expr.

   .. seealso:: :func:`foster`, :func:`foster_zRC2yRC`, :func:`dibujar_foster_serie`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.sintesis_dipolo import foster
   >>> from pytc2.dibujar import dibujar_foster_derivacion
   >>> s = sp.symbols('s ', complex=True)
   >>> # Sea la siguiente función de excitación
   >>> FF = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
   >>> # Se expande FF a la Foster
   >>> k0, koo, ki_wi, _, FF_foster = foster(FF)
   >>> # Tratamos a nuestra función imitancia como una Z
   >>> dibujar_foster_derivacion(k0 = k0, koo = koo, ki = ki_wi, y_exc = FF)


.. py:function:: dibujar_foster_serie(k0=None, koo=None, ki=None, kk=None, z_exc=None)

   Dibuja una red no disipativa a partir de una expansión en fracciones simples
   (Método de Foster). La forma matemática es:

   .. math:: Z(s)= \frac{k_0}{s} + k_\infty.s + \sum_{i=1}^N\frac{2.k_i.s}{s^2+\omega_i^2}

   Esta función provee una interpretación circuital al resultado de la función
   :func:`foster`.


   :param k0: Residuo de la función en DC o :math:`s \to 0`. El valor predeterminado es None.
   :type k0: simbólica, opcional
   :param koo: Residuo de la función en infinito o :math:`s \to \infty`. El valor predeterminado es None.
   :type koo: simbólica, opcional
   :param ki: Residuo de la función en :math:`\omega_i` o :math:`s^2 \to -\omega^2_i`. El valor predeterminado es None.
   :type ki: simbólica, list o tuple opcional
   :param kk: Residuo de la función en :math:`\sigma_i` o :math:`\omega \to -\omega_i`. El valor predeterminado es None.
   :type kk: simbólica, opcional

   :rtype: None

   :raises ValueError: Si cualquiera de los argumentos no son una instancia de sympy.Expr.

   .. seealso:: :func:`foster`, :func:`foster_zRC2yRC`, :func:`dibujar_foster_paralelo`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.sintesis_dipolo import foster
   >>> from pytc2.dibujar import dibujar_foster_serie
   >>> s = sp.symbols('s ', complex=True)
   >>> # Sea la siguiente función de excitación
   >>> FF = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)
   >>> # Se expande FF a la Foster
   >>> k0, koo, ki_wi, _, FF_foster = foster(FF)
   >>> # Tratamos a nuestra función imitancia como una Z
   >>> dibujar_foster_serie(k0 = k0, koo = koo, ki = ki_wi, z_exc = FF)


.. py:function:: dibujar_puerto_entrada(d, port_name=None, voltage_lbl=None, current_lbl=None)

   Dibuja un puerto de entrada a una red eléctrica diagramada mediante
   :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param port_name: Nombre del puerto. El valor predeterminado es None.
   :type port_name: string, opcional
   :param voltage_lbl: Etiqueta o nombre para la tensión del puerto. El valor predeterminado es None.
   :type voltage_lbl: string, tuple o list opcional
   :param current_lbl: Etiqueta o nombre para la corrientedel puerto. El valor predeterminado es None.
   :type current_lbl: string, opcional

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_funcion_exc_abajo`, :func:`dibujar_elemento_serie`, :func:`dibujar_puerto_salida`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_puerto_entrada, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Za")
   >>> d = dibujar_elemento_derivacion(d, "Z", sym_label="Zb")
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_puerto_salida(d, port_name=None, voltage_lbl=None, current_lbl=None)

   Dibuja un puerto de salida a una red eléctrica diagramada mediante
   :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param port_name: Nombre del puerto. El valor predeterminado es None.
   :type port_name: string, opcional
   :param voltage_lbl: Etiqueta o nombre para la tensión del puerto. El valor predeterminado es None.
   :type voltage_lbl: string, tuple o list opcional
   :param current_lbl: Etiqueta o nombre para la corrientedel puerto. El valor predeterminado es None.
   :type current_lbl: string, opcional

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_funcion_exc_abajo`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_puerto_entrada`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_puerto_entrada, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Za")
   >>> d = dibujar_elemento_derivacion(d, "Z", sym_label="Zb")
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_espaciador(d)

   Dibuja un espacio horizontal en un esquema dibujado mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_funcion_exc_abajo`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_puerto_entrada`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_espaciador, dibujar_puerto_entrada, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Za")
   >>> d = dibujar_espaciador(d)
   >>> d = dibujar_elemento_derivacion(d, "Z", sym_label="Zb")
   >>> d = dibujar_espaciador(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_funcion_exc_abajo(d, func_label, sym_func, k_gap_width=0.5, hacia_salida=False, hacia_entrada=False)

   Dibuja una ecuación correspondiente a la función de excitación definida en
   un dipolo de una red eléctrica diagramada mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param func_label: Etiqueta o nombre de la función de excitación.
   :type func_label: string
   :param sym_func: Un valor o expresión simbólica de la función `func_label` a indicar.
   :type sym_func: string, Real, symbolic expr.
   :param k_gap_width: Anchura del espacio destinado para la expresión proporcional a la escala del esquemático.
                       El valor predeterminado es `0.5*d.unit`.
   :type k_gap_width: Real, opcional
   :param hacia_salida: Booleano para indicar si la función se mide hacia la salida. El valor predeterminado es False.
   :type hacia_salida: boolean, opcional
   :param hacia_entrada: Booleano para indicar si la función se mide hacia la entrada. El valor predeterminado es False.
   :type hacia_entrada: string, opcional

   :returns: * **d** (*schemdraw.Drawing*) -- Objeto Drawing del módulo :mod:`schemdraw`.
             * **lbl** (*schemdraw.label*) -- Handle a la etiqueta visualizado.

   :raises None:

   .. seealso:: :func:`dibujar_funcion_exc_arriba`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_tanque_RC_serie`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> Za, Zb = sp.symbols('Za, Zb', complex=True)
   >>> # Sea la siguiente función de excitación
   >>> ZZ = Za+Zb
   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_funcion_exc_abajo, dibujar_puerto_entrada, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_funcion_exc_abajo(d,
   >>>                                  'Z',
   >>>                                  ZZ,
   >>>                                  hacia_salida = True)
   >>> d = dibujar_elemento_serie(d, "Z", Za)
   >>> d = dibujar_elemento_derivacion(d, "Z", Zb)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_funcion_exc_arriba(d, func_label, sym_func, k_gap_width=0.5, hacia_salida=False, hacia_entrada=False)

   Dibuja una ecuación correspondiente a la función de excitación definida en
   un dipolo de una red eléctrica diagramada mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param func_label: Etiqueta o nombre de la función de excitación.
   :type func_label: string
   :param sym_func: Un valor o expresión simbólica de la función `func_label` a indicar.
   :type sym_func: string, Real, symbolic expr.
   :param k_gap_width: Anchura del espacio destinado para la expresión proporcional a la escala del esquemático.
                       El valor predeterminado es `0.5*d.unit`.
   :type k_gap_width: Real, opcional
   :param hacia_salida: Booleano para indicar si la función se mide hacia la salida. El valor predeterminado es False.
   :type hacia_salida: boolean, opcional
   :param hacia_entrada: Booleano para indicar si la función se mide hacia la entrada. El valor predeterminado es False.
   :type hacia_entrada: string, opcional

   :returns: * **d** (*schemdraw.Drawing*) -- Objeto Drawing del módulo :mod:`schemdraw`.
             * **lbl** (*schemdraw.label*) -- Handle a la etiqueta visualizado.

   :raises None:

   .. seealso:: :func:`dibujar_funcion_exc_arriba`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_tanque_RC_serie`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> Za, Zb = sp.symbols('Za, Zb', complex=True)
   >>> # Sea la siguiente función de excitación
   >>> ZZ = Za+Zb
   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_funcion_exc_arriba, dibujar_puerto_entrada, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_funcion_exc_arriba(d,
   >>>                                  'Z',
   >>>                                  ZZ,
   >>>                                  hacia_salida = True)
   >>> d = dibujar_elemento_serie(d, "Z", Za)
   >>> d = dibujar_elemento_derivacion(d, "Z", Zb)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_elemento_serie(d, elemento, sym_label='')

   Dibuja un elemento en serie para una red eléctrica diagramada mediante
   :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param elemento: Un elemento a dibujar implementado en :mod:`schemdraw.elements` o un
                    string que apunte al elemento. Ej. 'R': Resistor,
                    'Z' o 'Y': ResistorIEC, 'C': Capacitor, 'L': Inductor, Line, Dot, Gap,
                    Arrow.
   :type elemento: str o elemento en schemdraw.elements
   :param sym_label: Un valor o expresión simbólica del elemento a dibujar.
   :type sym_label: string, Real, symbolic expr.

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_funcion_exc_arriba`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_tanque_RC_derivacion`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_puerto_entrada, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Za")
   >>> d = dibujar_elemento_derivacion(d, "Z", sym_label="Zb")
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_espacio_derivacion(d)

   Dibuja un espacio enb una red eléctrica diagramada mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_cierre`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_tanque_RC_derivacion`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_espacio_derivacion, dibujar_puerto_entrada, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Za")
   >>> d = dibujar_espacio_derivacion(d)
   >>> d = dibujar_elemento_derivacion(d, "Z", sym_label="Zb")
   >>> d = dibujar_espacio_derivacion(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_cierre(d)

   Dibuja un cierre entre el conductor superior e inferior en una red eléctrica
   diagramada mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_espacio_derivacion`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_tanque_RC_derivacion`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_cierre, dibujar_puerto_entrada, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Za")
   >>> d = dibujar_elemento_derivacion(d, "Z", sym_label="Zb")
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_cierre(d)
   >>> display(d)


.. py:function:: dibujar_elemento_derivacion(d, elemento, sym_label='', with_nodes=True)

   Dibuja un elemento en derivación para una red eléctrica diagramada mediante
   :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param elemento: Un elemento a dibujar implementado en :mod:`schemdraw`. Ej. Resistor,
                    ResistorIEC, Capacitor, Inductor, Line, Dot, Gap, Arrow.
   :type elemento: schemdraw.elements
   :param sym_label: Un valor o expresión simbólica del elemento a dibujar.
   :type sym_label: string, Real, symbolic expr.
   :param with_nodes = bool: Este booleano controla si la rama dibujada tendrá nodos o no. Es útil
                             al dibujar el primer elemento de una red, donde el nodo no suele ser
                             necesario.
   :param opcional: Este booleano controla si la rama dibujada tendrá nodos o no. Es útil
                    al dibujar el primer elemento de una red, donde el nodo no suele ser
                    necesario.

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_funcion_exc_arriba`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_tanque_RC_derivacion`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_espacio_derivacion, dibujar_puerto_entrada, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Za")
   >>> d = dibujar_espacio_derivacion(d)
   >>> d = dibujar_elemento_derivacion(d, "Z", sym_label="Zb")
   >>> d = dibujar_espacio_derivacion(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_tanque_RC_serie(d, resistor_label='', capacitor_lbl='')

   Dibuja un tanque RC (resistor y capacitor en paralelo) conectado en serie
   a una red eléctrica diagramada mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param resistor_label: Un valor o expresión simbólica del resistor a dibujar.
   :type resistor_label: string o symbolic expr.
   :param capacitor_lbl: Un valor o expresión simbólica del capacitor a dibujar.
   :type capacitor_lbl: string o symbolic expr.

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_funcion_exc_arriba`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_tanque_RC_derivacion`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_puerto_entrada, dibujar_tanque_RC_serie, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_tanque_RC_serie(d, "R_a", "C_a")
   >>> d = dibujar_elemento_derivacion(d, "Z", sym_label="Zb")
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_tanque_RC_derivacion(d, resistor_label='', capacitor_lbl='')

   Dibuja un tanque RC (resistor y capacitor en serie) conectado en derivación
   a una red eléctrica diagramada mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param resistor_label: Un valor o expresión simbólica del resistor a dibujar.
   :type resistor_label: string o symbolic expr.
   :param capacitor_lbl: Un valor o expresión simbólica del capacitor a dibujar.
   :type capacitor_lbl: string o symbolic expr.

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_tanque_RC_serie`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_funcion_exc_arriba`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_puerto_entrada, dibujar_tanque_RC_derivacion, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Za")
   >>> d = dibujar_tanque_RC_derivacion(d, "R_b", "C_b")
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_tanque_RL_serie(d, resistor_label='', inductor_label='')

   Dibuja un tanque RL (resistor e inductor en paralelo) conectado en serie
   a una red eléctrica diagramada mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param resistor_label: Un valor o expresión simbólica del resistor a dibujar.
   :type resistor_label: string o symbolic expr.
   :param inductor_label: Un valor o expresión simbólica del inductor a dibujar.
   :type inductor_label: string o symbolic expr.

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_funcion_exc_arriba`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_tanque_RL_derivacion`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_puerto_entrada, dibujar_tanque_RL_serie, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_tanque_RL_serie(d, "R_a", "L_a")
   >>> d = dibujar_elemento_derivacion(d, "Z", sym_label="Zb")
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_tanque_RL_derivacion(d, resistor_label='', inductor_label='')

   Dibuja un tanque RL (resistor e inductor en serie) conectado en derivación
   a una red eléctrica diagramada mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param resistor_label: Un valor o expresión simbólica del resistor a dibujar.
   :type resistor_label: string o symbolic expr.
   :param inductor_label: Un valor o expresión simbólica del inductor a dibujar.
   :type inductor_label: string o symbolic expr.

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_tanque_RL_serie`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_funcion_exc_arriba`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_puerto_entrada, dibujar_tanque_RL_derivacion, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Za")
   >>> d = dibujar_tanque_RL_derivacion(d, "R_b", "L_b")
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_tanque_serie(d, inductor_label='', capacitor_label='')

   Dibuja un tanque LC (inductor y capacitor en paralelo) conectado en serie
   a una red eléctrica diagramada mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param inductor_label: Un valor o expresión simbólica del inductor a dibujar.
   :type inductor_label: string o symbolic expr.
   :param capacitor_label: Un valor o expresión simbólica del capacitor a dibujar.
   :type capacitor_label: string o symbolic expr.

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_funcion_exc_arriba`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_tanque_RL_derivacion`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_puerto_entrada, dibujar_tanque_serie, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_tanque_serie(d, "L_a", "C_a")
   >>> d = dibujar_elemento_derivacion(d, "Z", sym_label="Zb")
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


.. py:function:: dibujar_tanque_derivacion(d, inductor_label='', capacitor_label='')

   Dibuja un tanque LC (inductor y capacitor en serie) conectado en derivación
   a una red eléctrica diagramada mediante :mod:`schemdraw`.


   :param d: Objeto Drawing del módulo :mod:`schemdraw`.
   :type d: schemdraw.Drawing
   :param inductor_label: Un valor o expresión simbólica del inductor a dibujar.
   :type inductor_label: string o symbolic expr.
   :param capacitor_label: Un valor o expresión simbólica del capacitor a dibujar.
   :type capacitor_label: string o symbolic expr.

   :returns: **d** -- Objeto Drawing del módulo :mod:`schemdraw`.
   :rtype: schemdraw.Drawing

   :raises None:

   .. seealso:: :func:`dibujar_tanque_serie`, :func:`dibujar_elemento_derivacion`, :func:`dibujar_tanque_RL_derivacion`

   .. rubric:: Examples

   >>> from schemdraw import Drawing
   >>> from pytc2.dibujar import dibujar_puerto_entrada, dibujar_tanque_derivacion, dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_puerto_salida
   >>> d = Drawing(unit=4)
   >>> d = dibujar_puerto_entrada(d)
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Za")
   >>> d = dibujar_tanque_derivacion(d, "L_a", "C_a")
   >>> d = dibujar_elemento_serie(d, "Z", sym_label="Zc")
   >>> d = dibujar_puerto_salida(d)
   >>> display(d)


