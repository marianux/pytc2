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
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
    
    assert isinstance(circ_hdl, io.TextIOBase)

    circ_hdl.close()

def test_nuevo_circuito_invalid():

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = test_module.ltsp_nuevo_circuito(1.)


def test_ltsp_capa_derivacion_invalid():

    circ_hdl = test_module.ltsp_nuevo_circuito('prueba1')

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = test_module.ltsp_capa_derivacion(circ_hdl, 2.0)

    circ_hdl.close()
