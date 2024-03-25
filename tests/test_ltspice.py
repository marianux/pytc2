#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
import io as io
from pytc2 import ltspice as test_module



def test_nuevo_circuito_valid():
        
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        circ_hdl = test_module.ltsp_nuevo_circuito('prueba1')
        assert isinstance(circ_hdl, io.TextIOBase)
        
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
    finally:
        if circ_hdl:
            circ_hdl.close()

def test_nuevo_circuito_invalid():

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = test_module.ltsp_nuevo_circuito(1.)

def test_ltsp_dibujar_valid():

    try:
        circ_hdl = test_module.ltsp_nuevo_circuito('prueba1')
        test_module.ltsp_etiquetar_nodo(circ_hdl)
        test_module.ltsp_ind_serie(circ_hdl, 1.0, ind_label='L1')
        test_module.ltsp_capa_derivacion(circ_hdl, 2.0, cap_label='C1')
        test_module.ltsp_capa_derivacion(circ_hdl, 2.0)
        test_module.ltsp_etiquetar_nodo(circ_hdl, node_label='vo')
        circ_hdl.writelines('TEXT -48 304 Left 2 !.param RG=1.0 RL=1.0')
        
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
    finally:
        if circ_hdl:
            circ_hdl.close()

def test_ltsp_dibujar_invalid():

    circ_hdl = test_module.ltsp_nuevo_circuito('prueba1')

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        test_module.ltsp_capa_derivacion(circ_hdl, '2.0')

    with pytest.raises(ValueError):
        test_module.ltsp_ind_serie(circ_hdl, '2.0')

    with pytest.raises(ValueError):
        test_module.ltsp_capa_derivacion(circ_hdl, -2.0)

    with pytest.raises(ValueError):
        test_module.ltsp_ind_serie(circ_hdl, -2.0)
    
    with pytest.raises(ValueError):
        test_module.ltsp_ind_serie(circ_hdl, 2.0, ind_label=1.0)
    
    with pytest.raises(ValueError):
        test_module.ltsp_capa_derivacion(circ_hdl, 2.0, cap_label=1.0)

    with pytest.raises(ValueError):
        test_module.ltsp_etiquetar_nodo(circ_hdl, node_label=1.0)

    circ_hdl.close()
