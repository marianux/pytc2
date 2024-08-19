pytc2.ltspice
=============

.. py:module:: pytc2.ltspice

.. autoapi-nested-parse::

   Created on Thu Mar  2 11:22:31 2023

   @author: mariano



Attributes
----------

.. autoapisummary::

   pytc2.ltspice.ltux
   pytc2.ltspice.ltuy
   pytc2.ltspice.filename_eq_base
   pytc2.ltspice.cap_num
   pytc2.ltspice.res_num
   pytc2.ltspice.ind_num
   pytc2.ltspice.node_num
   pytc2.ltspice.cur_x
   pytc2.ltspice.cur_y
   pytc2.ltspice.lt_wire_length
   pytc2.ltspice.res_der_str
   pytc2.ltspice.ind_der_str
   pytc2.ltspice.cap_der_str
   pytc2.ltspice.res_ser_str
   pytc2.ltspice.ind_ser_str
   pytc2.ltspice.cap_ser_str


Functions
---------

.. autoapisummary::

   pytc2.ltspice.ltsp_nuevo_circuito
   pytc2.ltspice.ltsp_capa_derivacion
   pytc2.ltspice.ltsp_ind_serie
   pytc2.ltspice.ltsp_etiquetar_nodo


Module Contents
---------------

.. py:data:: ltux
   :value: 16


   Unidades X para dibujar en la hoja de LTspice

.. py:data:: ltuy
   :value: 16


   Unidades Y para dibujar en la hoja de LTspice

.. py:data:: filename_eq_base
   :value: 'ltspice_equalizador_base'


   Archivo marco contenedor de las redes sintetizadas como ecualizadores/filtros

.. py:data:: cap_num
   :value: 1


   cuenta de capacitores

.. py:data:: res_num
   :value: 1


   cuenta de resistores

.. py:data:: ind_num
   :value: 1


   cuenta de inductores

.. py:data:: node_num
   :value: 1


   cuenta de nodos

.. py:data:: cur_x
   :value: 0


   cursor X para la localización de componentes

.. py:data:: cur_y
   :value: 0


   cursor Y para la localización de componentes

.. py:data:: lt_wire_length
   :value: 4


   tamaño estandard del cable

.. py:data:: res_der_str

   resistor en derivacion

.. py:data:: ind_der_str

   inductor en derivacion

.. py:data:: cap_der_str

   capacitor en derivacion

.. py:data:: res_ser_str

   resistor en serie

.. py:data:: ind_ser_str

   inductor en serie

.. py:data:: cap_ser_str

   capacitor en serie

.. py:function:: ltsp_nuevo_circuito(circ_name=None, circ_folder=None)

   Se genera un circuito nuevo en LTspice de nombre *circ_name*.

   :param circ_name: Nombre del circuito.
   :type circ_name: string
   :param circ_folder: Path a la carpeta donde se creará el archivo ASC y PLT de LTspice.
   :type circ_folder: str, opcional

   :returns: **circ_hdl** -- Handle al archivo de texto de LTspice para continuar construyendo el
             circuito.
   :rtype: archivo de texto

   :raises TypeError: Si ZZ no es una instancia de sympy.Matrix.

   .. seealso:: :func:`ltsp_capa_derivacion`, :func:`ltsp_ind_serie`

   .. rubric:: Examples

   >>> from pytc2.ltspice import ltsp_nuevo_circuito, ltsp_etiquetar_nodo, ltsp_ind_serie, ltsp_capa_derivacion, ltsp_etiquetar_nodo
   >>> circ_hdl = ltsp_nuevo_circuito('prueba1')
   >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vi')
   >>> ltsp_ind_serie(circ_hdl, 1.0)
   >>> ltsp_capa_derivacion(circ_hdl, 2.0)
   >>> ltsp_ind_serie(circ_hdl, 1.0)
   >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vo')
   >>> R01 = 1.0
   >>> R02 = 1.0
   >>> circ_hdl.writelines('TEXT -48 304 Left 2 !.param RG={:3.3f} RL={:3.3f}'.format(R01, R02))
   >>> circ_hdl.close()
   [ Buscar el archivo "ltsp_prueba.asc" en LTspice ]


.. py:function:: ltsp_capa_derivacion(circ_hdl, cap_value, cap_label=None)

   Incorpora un capacitor en derivación a un circuito en LTspice.

   :param circ_hdl: Handle al archivo LTspice.
   :type circ_hdl: archivo de texto LTspice
   :param cap_value: Valor del capacitor.
   :type cap_value: float o numéro simbólico
   :param cap_label: Etiqueta para identificar al capacitor en el circuito.
   :type cap_label: string o None

   :rtype: None

   :raises ValueError: Si cap_value no es numérico o el valor no es positivo.

   .. seealso:: :func:`ltsp_capa_derivacion`, :func:`ltsp_ind_serie`

   .. rubric:: Examples

   >>> from pytc2.ltspice import ltsp_nuevo_circuito, ltsp_etiquetar_nodo, ltsp_ind_serie, ltsp_capa_derivacion, ltsp_etiquetar_nodo
   >>> circ_hdl = ltsp_nuevo_circuito('prueba1')
   >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vi')
   >>> ltsp_ind_serie(circ_hdl, 1.0)
   >>> ltsp_capa_derivacion(circ_hdl, 2.0)
   >>> ltsp_ind_serie(circ_hdl, 1.0)
   >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vo')
   >>> R01 = 1.0
   >>> R02 = 1.0
   >>> circ_hdl.writelines('TEXT -48 304 Left 2 !.param RG={:3.3f} RL={:3.3f}'.format(R01, R02))
   >>> circ_hdl.close()
   [ Buscar el archivo "ltsp_prueba.asc" en LTspice ]


.. py:function:: ltsp_ind_serie(circ_hdl, ind_value, ind_label=None)

   Incorpora un inductor en serie a un circuito en LTspice.

   :param circ_hdl: Handle al archivo LTspice.
   :type circ_hdl: archivo de texto LTspice
   :param ind_value: Valor del inductor.
   :type ind_value: float o numéro simbólico
   :param ind_label: Etiqueta para identificar al inductor en el circuito.
   :type ind_label: string o None

   :rtype: None

   :raises ValueError: Si cap_value no es numérico o el valor no es positivo.

   .. seealso:: :func:`ltsp_capa_derivacion`, :func:`ltsp_ind_serie`

   .. rubric:: Examples

   >>> from pytc2.ltspice import ltsp_nuevo_circuito, ltsp_etiquetar_nodo, ltsp_ind_serie, ltsp_capa_derivacion, ltsp_etiquetar_nodo
   >>> circ_hdl = ltsp_nuevo_circuito('prueba1')
   >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vi')
   >>> ltsp_ind_serie(circ_hdl, 1.0)
   >>> ltsp_capa_derivacion(circ_hdl, 2.0)
   >>> ltsp_ind_serie(circ_hdl, 1.0)
   >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vo')
   >>> R01 = 1.0
   >>> R02 = 1.0
   >>> circ_hdl.writelines('TEXT -48 304 Left 2 !.param RG={:3.3f} RL={:3.3f}'.format(R01, R02))
   >>> circ_hdl.close()
   [ Buscar el archivo "ltsp_prueba.asc" en LTspice ]


.. py:function:: ltsp_etiquetar_nodo(circ_hdl, node_label=None)

   Asigna una etiqueta a un nodo de un circuito en LTspice.

   :param circ_hdl: Handle al archivo LTspice.
   :type circ_hdl: archivo de texto LTspice
   :param node_label: Etiqueta para identificar al nodo en el circuito.
   :type node_label: string o None

   :rtype: None

   :raises ValueError: Si cap_value no es numérico o el valor no es positivo.

   .. seealso:: :func:`ltsp_capa_derivacion`, :func:`ltsp_ind_serie`

   .. rubric:: Examples

   >>> from pytc2.ltspice import ltsp_nuevo_circuito, ltsp_etiquetar_nodo, ltsp_ind_serie, ltsp_capa_derivacion, ltsp_etiquetar_nodo
   >>> circ_hdl = ltsp_nuevo_circuito('prueba1')
   >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vi')
   >>> ltsp_ind_serie(circ_hdl, 1.0)
   >>> ltsp_capa_derivacion(circ_hdl, 2.0)
   >>> ltsp_ind_serie(circ_hdl, 1.0)
   >>> ltsp_etiquetar_nodo(circ_hdl, node_label='vo')
   >>> R01 = 1.0
   >>> R02 = 1.0
   >>> circ_hdl.writelines('TEXT -48 304 Left 2 !.param RG={:3.3f} RL={:3.3f}'.format(R01, R02))
   >>> circ_hdl.close()
   [ Buscar el archivo "ltsp_prueba.asc" en LTspice ]


