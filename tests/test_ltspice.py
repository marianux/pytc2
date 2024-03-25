#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
import sympy as sp
from pytc2 import dibujar as test_module
import matplotlib.pyplot as plt
from pytc2.sintesis_dipolo import cauer_RC, cauer_LC, foster
from schemdraw import Drawing

k, o = sp.symbols('k, o')
s = sp.symbols('s ', complex=True)
z1, z2 = sp.symbols('z1, z2')

@pytest.mark.parametrize(
    "func_ptr, args, true_val",
    [
       (   
           test_module.dibujar_Tee,
           (sp.Matrix([[5, 2], [2, 6]]), True),
           [3, 2, 4]
        ),
       (   
           test_module.dibujar_Pi,
           (sp.Matrix([[5, -2], [-2, 6]]), True),
           [3, 2, 4]
        ),
       (   
           test_module.dibujar_lattice,
           (sp.Matrix([[5, 2], [2, 6]]), True),
           [3, 7]
        )
       
   ]
)   
def test_dibujar_valid(func_ptr, args, true_val):
    
    
    p1, p2 = args

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        out = func_ptr(p1, p2)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    for result, truev in zip(out, true_val):
        assert sp.simplify(result - truev).is_zero 


@pytest.mark.parametrize(
    "func_ptr",
    [
        test_module.dibujar_Tee,
        test_module.dibujar_Pi,
        test_module.dibujar_Pi
    ]
)   
def test_dibujar_invalid_input(func_ptr):

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr(sp.Matrix([[5, 2], [2, 6]]), 'a')
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a')




# cerramos todas las figuras que se hayan abierto, para no perjudicar otras pruebas
def test_last_call():

    plt.close('all')

ZRC = (s**2 + 4*s + 3)/(s**2 + 2*s)

# Implementaremos FF mediante Cauer 1 o remociones continuas en infinito
koo_RC, ZRC_cauer_oo_RC, _ = cauer_RC(ZRC, remover_en_inf=True)

FF = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)

# Implementaremos FF mediante Cauer 1 o remociones continuas en infinito
koo_LC, F_cauer_oo_LC, _ = cauer_LC(FF, remover_en_inf=True)

# Se expande FF a la Foster
k0_fd, koo_fd, ki_wi_fd, _, FF_foster_fd = foster(FF)


@pytest.mark.parametrize(
    "funciones_parametros",
    [
        {
            test_module.dibujar_cauer_RC_RL: {'ki' : koo_RC, 'y_exc' : ZRC_cauer_oo_RC, 'z_exc' : None},
        },
        {
            test_module.dibujar_cauer_RC_RL: {'ki' : koo_RC, 'z_exc' : ZRC_cauer_oo_RC, 'y_exc' : None},
        },
        {
            test_module.dibujar_cauer_LC:    {'ki' : koo_LC, 'y_exc' : F_cauer_oo_LC, 'z_exc' : None},
        },
        {
            test_module.dibujar_cauer_LC:    {'ki' : koo_LC, 'z_exc' : F_cauer_oo_LC, 'y_exc' : None},
        },
        {
            test_module.dibujar_foster_derivacion:    {'k0' : k0_fd, 'koo' : koo_fd, 'ki' : ki_wi_fd, 'kk' : None, 'y_exc' : FF},
        },
        {
            test_module.dibujar_foster_serie:    {'k0' : k0_fd, 'koo' : koo_fd, 'ki' : ki_wi_fd, 'kk' : None, 'z_exc' : FF},
        },
        
   ]
)   
def test_dibujar_canonicas_valid(funciones_parametros):

    for funcion, parametros in funciones_parametros.items():

        # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
        try:
            funcion(**parametros)
        except ValueError:
            pytest.fail("Se levantó un ValueError incorrectamente.")
        
        
@pytest.mark.parametrize(
    "funciones_parametros",
    [
        {
            test_module.dibujar_cauer_RC_RL: {'ki' : 'koo_RC', 'y_exc' : ZRC_cauer_oo_RC, 'z_exc' : None},
        },
        {
            test_module.dibujar_cauer_RC_RL: {'ki' : koo_RC, 'y_exc' : 'ZRC_cauer_oo_RC'},
        },
        {
            test_module.dibujar_cauer_RC_RL: {'ki' : koo_RC, 'z_exc' : 'ZRC_cauer_oo_RC'},
        },
        {
            test_module.dibujar_cauer_LC:    {'ki' : 'koo_LC', 'y_exc' : F_cauer_oo_LC, 'z_exc' : None},
        },
        {
            test_module.dibujar_cauer_LC:    {'ki' : koo_LC, 'y_exc' : 'F_cauer_oo_LC', 'z_exc' : None},
        },
        {
            test_module.dibujar_cauer_LC:    {'ki' : koo_LC, 'z_exc' : 'F_cauer_oo_LC'},
        },
        {
            test_module.dibujar_foster_derivacion:    {'k0' : 'k0_fd', 'koo' : koo_fd, 'ki' : ki_wi_fd, 'kk' : None, 'y_exc' : FF},
        },
        {
            test_module.dibujar_foster_derivacion:    {'k0' : k0_fd, 'koo' : 'koo_fd', 'ki' : ki_wi_fd, 'kk' : None, 'y_exc' : FF},
        },
        {
            test_module.dibujar_foster_derivacion:    {'k0' : k0_fd, 'koo' : koo_fd, 'ki' : 'ki_wi_fd', 'kk' : None, 'y_exc' : FF},
        },
        {
            test_module.dibujar_foster_derivacion:    {'k0' : k0_fd, 'koo' : koo_fd, 'ki' : ki_wi_fd, 'kk' : 'None', 'y_exc' : FF},
        },
        {
            test_module.dibujar_foster_derivacion:    {'k0' : k0_fd, 'koo' : koo_fd, 'ki' : ki_wi_fd, 'kk' : None, 'y_exc' : 'FF'},
        },
        {
            test_module.dibujar_foster_serie:    {'k0' : 'k0_fd', 'koo' : koo_fd, 'ki' : ki_wi_fd, 'kk' : None, 'z_exc' : FF},
        },
        {
            test_module.dibujar_foster_serie:    {'k0' : k0_fd, 'koo' : 'koo_fd', 'ki' : ki_wi_fd, 'kk' : None, 'z_exc' : FF},
        },
        {
            test_module.dibujar_foster_serie:    {'k0' : k0_fd, 'koo' : koo_fd, 'ki' : 'ki_wi_fd', 'kk' : None, 'z_exc' : FF},
        },
        {
            test_module.dibujar_foster_serie:    {'k0' : k0_fd, 'koo' : koo_fd, 'ki' : ki_wi_fd, 'kk' : 'None', 'z_exc' : FF},
        },
        {
            test_module.dibujar_foster_serie:    {'k0' : k0_fd, 'koo' : koo_fd, 'ki' : ki_wi_fd, 'kk' : None, 'z_exc' : 'FF'},
        },
        
   ]
)   
def test_dibujar_canonicas_invalid(funciones_parametros):

    for funcion, parametros in funciones_parametros.items():

        # Verificar que se levante un ValueError 
        with pytest.raises(ValueError):
            funcion(**parametros)


d = Drawing(unit=4)

@pytest.mark.parametrize(
    "funciones_parametros",
    [
        {
            test_module.dibujar_puerto_entrada: {'d': d},
        },
        {
            test_module.dibujar_puerto_entrada: {'d': d, 'port_name':'port_name', 'voltage_lbl' : 'v_name', 'current_lbl' : 'i_name'},
        },
        {
            test_module.dibujar_puerto_salida: {'d': d, },
        },
        {
            test_module.dibujar_puerto_salida: {'d': d,  'port_name':'port_name', 'voltage_lbl' : 'v_name', 'current_lbl' : 'i_name'},
        },
        {
            test_module.dibujar_espaciador: {'d': d },
        },
        {
            test_module.dibujar_funcion_exc_abajo: {'d': d, 'func_label':'Z_{in}', 'sym_func': ZRC, 'k_gap_width': 1., 'hacia_salida': True, 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_arriba: {'d': d, 'func_label':'Z_{in}', 'sym_func': ZRC, 'k_gap_width': 1., 'hacia_salida': True, 'hacia_entrada': True},
        },
        {
            test_module.dibujar_elemento_serie: {'d': d, 'elemento': 'Z', 'sym_label': 'Za' },
        },
        {
            test_module.dibujar_espacio_derivacion: {'d': d },
        },
        {
            test_module.dibujar_cierre: {'d': d },
        },
        {
            test_module.dibujar_elemento_derivacion: {'d': d, 'elemento': 'Z', 'sym_label': 'Za' },
        },
        {
            test_module.dibujar_tanque_RC_serie: {'d': d, 'resistor_label': 'R', 'capacitor_lbl': 'C' },
        },
        {
            test_module.dibujar_tanque_RC_derivacion: {'d': d, 'resistor_label': 'R', 'capacitor_lbl': 'C' },
        },
        {
            test_module.dibujar_tanque_RL_serie: {'d': d, 'resistor_label': 'R', 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_RL_derivacion: {'d': d, 'resistor_label': 'R', 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_serie: {'d': d, 'capacitor_label': 'C', 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_derivacion: {'d': d, 'capacitor_label': 'C', 'inductor_label': 'L' },
        },
    ]
)   
def test_dibujar_redes_valid(funciones_parametros):

    for funcion, parametros in funciones_parametros.items():

        # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
        try:
            d = funcion(**parametros)
        except ValueError:
            pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(d, Drawing) 
        
        
@pytest.mark.parametrize(
    "funciones_parametros",
    [
        {
            test_module.dibujar_puerto_entrada: {'d': 'd', 'port_name':'port_name', 'voltage_lbl' : 'v_name', 'current_lbl' : 'i_name'},
        },
        {
            test_module.dibujar_puerto_entrada: {'d': d, 'port_name': 1., 'voltage_lbl' : 'v_name', 'current_lbl' : 'i_name'},
        },
        {
            test_module.dibujar_puerto_entrada: {'d': d, 'port_name':'port_name', 'voltage_lbl' : 1., 'current_lbl' : 'i_name'},
        },
        {
            test_module.dibujar_puerto_entrada: {'d': d, 'port_name':'port_name', 'voltage_lbl' : 'v_name', 'current_lbl' : 1.},
        },
        {
            test_module.dibujar_puerto_entrada: {'d': 'd'},
        },
        {
            test_module.dibujar_puerto_salida: {'d': 'd', 'port_name':'port_name', 'voltage_lbl' : 'v_name', 'current_lbl' : 'i_name'},
        },
        {
            test_module.dibujar_puerto_salida: {'d': d, 'port_name': 1., 'voltage_lbl' : 'v_name', 'current_lbl' : 'i_name'},
        },
        {
            test_module.dibujar_puerto_salida: {'d': d, 'port_name':'port_name', 'voltage_lbl' : 1., 'current_lbl' : 'i_name'},
        },
        {
            test_module.dibujar_puerto_salida: {'d': d, 'port_name':'port_name', 'voltage_lbl' : 'v_name', 'current_lbl' : 1.},
        },
        {
            test_module.dibujar_puerto_salida: {'d': 'd'},
        },
        {
            test_module.dibujar_espaciador: {'d': 'd' },
        },
        {
            test_module.dibujar_funcion_exc_abajo: {'d': 'd', 'func_label':'Z_{in}', 'sym_func': ZRC, 'k_gap_width': 1., 'hacia_salida': True, 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_abajo: {'d': d, 'func_label':1., 'sym_func': ZRC, 'k_gap_width': 1., 'hacia_salida': True, 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_abajo: {'d': d, 'func_label':'Z_{in}',                  'k_gap_width': 1., 'hacia_salida': True, 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_abajo: {'d': d, 'func_label':'Z_{in}', 'sym_func': ZRC, 'k_gap_width': '1.', 'hacia_salida': True, 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_abajo: {'d': d, 'func_label':'Z_{in}', 'sym_func': ZRC, 'k_gap_width': 1., 'hacia_salida': 'True', 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_abajo: {'d': d, 'func_label':'Z_{in}', 'sym_func': ZRC, 'k_gap_width': 1., 'hacia_salida': True, 'hacia_entrada': 'True'},
        },
        {
            test_module.dibujar_funcion_exc_abajo: {'d': d },
        },
        {
            test_module.dibujar_funcion_exc_arriba: {'d': 'd', 'func_label':'Z_{in}', 'sym_func': ZRC, 'k_gap_width': 1., 'hacia_salida': True, 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_arriba: {'d': d, 'func_label':1., 'sym_func': ZRC, 'k_gap_width': 1., 'hacia_salida': True, 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_arriba: {'d': d, 'func_label':'Z_{in}',                  'k_gap_width': 1., 'hacia_salida': True, 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_arriba: {'d': d, 'func_label':'Z_{in}', 'sym_func': ZRC, 'k_gap_width': '1.', 'hacia_salida': True, 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_arriba: {'d': d, 'func_label':'Z_{in}', 'sym_func': ZRC, 'k_gap_width': 1., 'hacia_salida': 'True', 'hacia_entrada': True},
        },
        {
            test_module.dibujar_funcion_exc_arriba: {'d': d, 'func_label':'Z_{in}', 'sym_func': ZRC, 'k_gap_width': 1., 'hacia_salida': True, 'hacia_entrada': 'True'},
        },
        {
            test_module.dibujar_funcion_exc_arriba: {'d': d },
        },
        {
            test_module.dibujar_elemento_serie: {'d': 'd', 'elemento': 'Z', 'sym_label': 'Za' },
        },
        {
            test_module.dibujar_elemento_serie: {'d': d, 'elemento': 'X', 'sym_label': 'Za' },
        },
        {
            test_module.dibujar_espacio_derivacion: {'d': 'd' },
        },
        {
            test_module.dibujar_cierre: {'d': 'd' },
        },
        {
            test_module.dibujar_elemento_derivacion: {'d': 'd', 'elemento': 'Y', 'sym_label': 'Ya' },
        },
        {
            test_module.dibujar_elemento_derivacion: {'d': d, 'elemento': 'X', 'sym_label': 'Ya' },
        },
        {
            test_module.dibujar_tanque_RC_serie: {'d': 'd', 'resistor_label': 'R', 'capacitor_lbl': 'C' },
        },
        {
            test_module.dibujar_tanque_RC_serie: {'d': d, 'resistor_label': 1., 'capacitor_lbl': 'C' },
        },
        {
            test_module.dibujar_tanque_RC_serie: {'d': d, 'resistor_label': 'R', 'capacitor_lbl': 1. },
        },
        {
            test_module.dibujar_tanque_RC_derivacion: {'d': 'd', 'resistor_label': 'R', 'capacitor_lbl': 'C' },
        },
        {
            test_module.dibujar_tanque_RC_derivacion: {'d': d, 'resistor_label': 1., 'capacitor_lbl': 'C' },
        },
        {
            test_module.dibujar_tanque_RC_derivacion: {'d': d, 'resistor_label': 'R', 'capacitor_lbl': 1. },
        },
        {
            test_module.dibujar_tanque_RL_serie: {'d': 'd', 'resistor_label': 'R', 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_RL_serie: {'d': d, 'resistor_label': 1., 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_RL_serie: {'d': d, 'resistor_label': 'R', 'inductor_label': 1. },
        },
        {
            test_module.dibujar_tanque_RL_derivacion: {'d': 'd', 'resistor_label': 'R', 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_RL_derivacion: {'d': d, 'resistor_label': 1., 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_RL_derivacion: {'d': d, 'resistor_label': 'R', 'inductor_label': 1. },
        },
        {
            test_module.dibujar_tanque_serie: {'d': 'd', 'capacitor_label': 'C', 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_serie: {'d': d, 'capacitor_label': 1., 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_serie: {'d': d, 'capacitor_label': 'C', 'inductor_label': 1. },
        },
        {
            test_module.dibujar_tanque_derivacion: {'d': 'd', 'capacitor_label': 'C', 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_derivacion: {'d': d, 'capacitor_label': 1., 'inductor_label': 'L' },
        },
        {
            test_module.dibujar_tanque_derivacion: {'d': d, 'capacitor_label': 'C', 'inductor_label': 1. },
        },
    ]
)   

def test_dibujar_redes_invalid(funciones_parametros):

    for funcion, parametros in funciones_parametros.items():

        # Verificar que se levante un ValueError 
        with pytest.raises((ValueError, TypeError)):
            funcion(**parametros)


# def test__valid():
    
# def test__invalid_input():
