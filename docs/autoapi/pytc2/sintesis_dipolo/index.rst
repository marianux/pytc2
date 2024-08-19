pytc2.sintesis_dipolo
=====================

.. py:module:: pytc2.sintesis_dipolo

.. autoapi-nested-parse::

   Created on Thu Mar  2 11:26:00 2023

   @author: mariano



Functions
---------

.. autoapisummary::

   pytc2.sintesis_dipolo.cauer_RC
   pytc2.sintesis_dipolo.cauer_LC
   pytc2.sintesis_dipolo.foster_zRC2yRC
   pytc2.sintesis_dipolo.foster


Module Contents
---------------

.. py:function:: cauer_RC(imm, remover_en_inf=True)

   Realiza una expansión en fracciones continuas sobre una inmitancia (imm),
   removiendo en DC o :math:`\infty` dependiendo de (remover_en_inf). Este
   procedimiento se conoce como métodos de Cauer I y II. En el ejemplo de
   :math:`Z_{RC}` se remueve en DC y para el caso de :math:`Y_{RC}`
   en :math:`\infty`.

   .. math:: Z_{RC}(s)= \frac{1}{s.C_1} + \frac{1}{ \frac{1}{R_1} + \frac{1}{ \frac{1}{s.C_2} + \cdots } } =
        R_1 + \frac{1}{ s.C_1 + \frac{1}{ R_2 + \cdots } }

   .. math:: Y_{RC}(s)= s.C_1 + \frac{1}{ R_1 + \frac{1}{ s.C_2 + \cdots } } =
        \frac{1}{R_1} + \frac{1}{ s.C_1 + \frac{1}{ \frac{1}{R_2} + \cdots } }

   :param imm: La inmitancia a expandir en fracciones continuas..
   :type imm: symbolic rational function
   :param remover_en_inf: Determina en qué extremo se realiza la remoción.
   :type remover_en_inf: boolean

   :rtype: A list k0 with the i-th k0_i resulted from continued fraction expansion.

   :raises ValueError: Si y_exc y z_exc no son una instancia de sympy.Expr.

   .. seealso:: :func:`cauer_LC`, :func:`dibujar_cauer_RC_RL`, :func:`dibujar_cauer_LC`

   .. rubric:: Example

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


.. py:function:: cauer_LC(imm, remover_en_inf=True)

   Dibuja una red escalera no disipativa, a partir de la expansión en fracciones
   continuas (Método de Cauer). Se remueve en DC o :math:`\infty` dependiendo
   de *remover_en_inf*. En los siguientes ejemplos se expande tanto :math:`Z(s)`
   como :math:`Y(s)`, y se remueve a la izquierda en DC y a la derecha en :math:`\infty`.
   La forma matemática será:

   .. math:: Z(s)= \frac{1}{s.C_1} + \frac{1}{ \frac{1}{s.L_1} + \frac{1}{ \frac{1}{s.C_2} + \cdots } } =
            s.L_1 + \frac{1}{ s.C_1 + \frac{1}{ s.L_2 + \cdots } }

   .. math:: Y(s)= \frac{1}{s.L_1} + \frac{1}{ \frac{1}{s.C_1} + \frac{1}{ \frac{1}{s.L_2} + \cdots } } =
            s.C_1 + \frac{1}{ s.L_1 + \frac{1}{ s.C_2 + \cdots } }


   :param imm: La inmitancia a expandir en fracciones continuas..
   :type imm: symbolic rational function
   :param remover_en_inf: Determina en qué extremo se realiza la remoción.
   :type remover_en_inf: boolean

   :returns: * **ko** (*lista de expresiones simbólicas*) -- Conjunto de términos con los residuos de forma :math:`\frac{k_0}{s}` y :math:`s.k_{\infty}`
             * **imm_as_cauer** (*symbolic rational function*) -- La función inmitancia expandida en fracciones contínuas.
             * **rem** (*symbolic rational function*) -- 0 en caso que la expansión sea exitosa, ó una función remanente que no
               puede ser expresada en formato Cauer.

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


.. py:function:: foster_zRC2yRC(k0=sp.Rational(0), koo=sp.Rational(0), ki_wi=sp.Rational(0), kk=sp.Rational(0), ZRC_foster=sp.Rational(0))

   Permite llegar a la forma foster de una inmitancia :math:`I(s)` (YRC - ZRL),
   a partir de la propia función :func:`foster` de expansión en fracciones
   simples, y una conversión término-a-término de cada residuo obtenido.

   De esa manera se comienza con la expansión foster( I(s)/s ), para luego
   realizarel siguiente mapeo de residuos:

   + :math:`k_\infty = kk`
   + :math:`k_k = k_0`
   + :math:`k_i = ki_wi`
   + :math:`I_F(s) = I(s)*s`

   :param k0: Residuo de la función en DC o :math:`s \to 0`. El valor predeterminado es 0.
   :type k0: simbólica, opcional
   :param koo: Residuo de la función en infinito o :math:`s \to \infty`. El valor predeterminado es 0.
   :type koo: simbólica, opcional
   :param ki_wi: Residuo de la función en :math:`\omega_i` o :math:`s^2 \to -\omega^2_i`. El valor predeterminado es 0.
   :type ki_wi: simbólica, list o tuple opcional
   :param kk: Residuo de la función en :math:`\sigma_i` o :math:`\omega \to -\omega_i`. El valor predeterminado es 0.
   :type kk: simbólica, opcional
   :param ZRC_foster: Función inmitancia :math:`I(s)` a expresar como :math:`I_F(s)`
   :type ZRC_foster: simbólica

   :returns: * **k0** (*simbólica, opcional*) -- No está permitido para esta forma el residuo en 0.
             * **koo** (*simbólica, opcional*) -- Residuo de la función en infinito o :math:`s \to \infty`.
             * **ki** (*simbólica, list o tuple opcional*) -- Residuo de la función en :math:`\omega_i` o :math:`s \to -\sigma_i`.
             * **kk** (*simbólica, opcional*) -- Residuo de la función en :math:`\sigma = 0`.
             * **YRC_foster** (*simbólica*) -- Función YRC expresada como :math:`I_F(s) = I(s)*s`

   :raises ValueError: Si cualquiera de los argumentos no son una instancia de sympy.Expr.

   .. seealso:: :func:`foster`, :func:`foster_zRC2yRC`, :func:`dibujar_foster_serie`

   .. rubric:: Examples

   >>> import sympy as sp
   >>> from pytc2.sintesis_dipolo import foster, foster_zRC2yRC
   >>> from pytc2.dibujar import dibujar_foster_derivacion
   >>> s = sp.symbols('s ', complex=True)
   >>> # Sea la siguiente función de excitación
   >>> YRC = 2*(s**2 + 4*s + 3)/(s**2 + 8*s + 12)
   >>> k0, koo, ki_wi, kk, YRC_foster = foster(YRC/s)
   >>> k0, koo, ki_wi, kk, YRC_foster = foster_zRC2yRC(k0, koo, ki_wi, kk, YRC_foster)
   >>> dibujar_foster_derivacion(k0 = k0, koo = koo, ki = ki_wi, y_exc = YRC_foster)


.. py:function:: foster(imm)

   Expande una función inmitancia :math:`I(s)` en fracciones simples, de acuerdo al método
   de Foster. La forma matemática es:

   .. math:: I(s)= \frac{k_0}{s} + k_\infty.s + \sum_{i=1}^N\frac{2.k_i.s}{s^2+\omega_i^2}

   Dependiendo la naturaleza de :math:`I(s)` como impedancia o admitancia,
   resultará en los métodos de Foster serie, o paralelo. También existen 3
   variantes 1) en caso que se trate de redes no disipativas (LC), y redes
   disipativas compuestos solo por dos elementos circuitales: RC - RL. 2) Las
   expresiones matemáticas para :math:`Z_{RC}` son las mismas que :math:`Y_{RL}`,
   mientras que 3) las de :math:`Z_{RL}` iguales a las de :math:`Y_{RC}`.


   :param k0: Residuo de la función en DC o :math:`s \to 0`. El valor predeterminado es 0.
   :type k0: simbólica, opcional
   :param koo: Residuo de la función en infinito o :math:`s \to \infty`. El valor predeterminado es 0.
   :type koo: simbólica, opcional
   :param ki: Residuo de la función en :math:`\omega_i` o :math:`s^2 \to -\omega^2_i`. El valor predeterminado es 0.
   :type ki: simbólica, list o tuple opcional
   :param kk: Residuo de la función en :math:`\sigma_i` o :math:`\omega \to -\omega_i`. El valor predeterminado es 0.
   :type kk: simbólica, opcional

   :returns: * **k0** (*simbólica, opcional*) -- El residuo en 0, expresado matemáticamente como :math:`\frac{k_0}{s}`.
             * **koo** (*simbólica, opcional*) -- Residuo de la función en infinito o :math:`s \to \infty`, que se
               corresponde al término :math:`k_\infty*s`.
             * **ki** (*simbólica, list o tuple opcional*) -- Residuo de la función en :math:`s^2 \to -\omega_i^2`, matemáticamente
               :math:`\frac{2.k_i.s}{s^2+\omega_i^2}`.
             * **kk** (*simbólica, opcional*) -- Residuo de la función en :math:`\sigma = 0`, para funciones disipativas.
             * **foster_form** (*simbólica*) -- Función YRC expresada como :math:`I_F(s) = I(s)*s`

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


