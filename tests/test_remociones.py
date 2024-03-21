#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
import sympy as sp
from pytc2 import remociones as test_module


k, o = sp.symbols('k, o')

@pytest.mark.parametrize(
    "func_ptr, true_val",
    [
       (   
           test_module.tanque_z,
           ( k/o, 1/k )
        ), 
       (   
           test_module.tanque_y,
           ( 1/k, k/o )
        ), 
      
    ]
)
def test_tanques_valid(func_ptr, true_val):
    
    L_t, C_t = true_val

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        L, C = func_ptr(k,o)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(L, sp.Expr)
    assert isinstance(C, sp.Expr)

    # Verificar el correcto resultado
    assert sp.simplify(L_t - L) == sp.Rational(0)
    
    assert sp.simplify(C_t - C) == sp.Rational(0)


@pytest.mark.parametrize(
    "func_ptr",
    [
       # S -> Ts
        test_module.tanque_z,
        test_module.tanque_y,
    ]
)
def test_tanques_invalid_input(func_ptr):

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        L, C = func_ptr(k, 1.)
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        L, C = func_ptr(1., k)
        
s = sp.symbols('s ', complex=True)

@pytest.mark.parametrize(
    "func_ptr, args, true_val",
    [
       (   
           test_module.trim_poly_s,
           sp.poly( 1e-10*s**3 + 2*s**2 + s + 1 , s),
           2.0*s**2 + 1.0*s + 1.0
        ), 
       (   
           test_module.trim_func_s,
           ( 1e-10*s**3 + 2*s**2 + s + 1)/( 4.3e-10*s**2 + 2*s + 5),
           (2.0*s**2 + 1.0*s + 1.0)/(2.0*s + 5.0)

        ), 
      
    ]
)   
def test_trims_valid(func_ptr, args, true_val):
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        out = func_ptr(args)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(out, sp.Expr)

    # Verificar el correcto resultado
    assert sp.simplify(out - true_val) == sp.Rational(0)
    
    
@pytest.mark.parametrize(
    "func_ptr, args",
    [
       (   
           test_module.trim_poly_s,
           sp.poly( 1e-10*s**3 + 2*s**2 + s + 1 , s)
        ), 
       (   
           test_module.trim_func_s,
           ( 1e-10*s**3 + 2*s**2 + s + 1)/( 4.3e-10*s**2 + 2*s + 5)

        ), 
      
    ]
)   
def test_trims_invalid_input(func_ptr, args):
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr(1., tol = 1e-6)
    
    with pytest.raises(ValueError):
        out = func_ptr(args, tol = '1e-6')


# TODO: todavía la función modsq2mod_s no funciona correctamente. Probar el ejemplo.
# def test_modsq2mod_s_valid():
    
#     this_func = (  s**4 + 6*s**2 + 9)/( s**4 - 2*s**2 + 1)
    
#     true_val = (s**2 + 3)/(s**2 + 2*s + 1)
    
#     # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
#     try:
#         factor_func = test_module.modsq2mod_s( this_func )
#     except ValueError:
#         pytest.fail("Se levantó un ValueError incorrectamente.")

#     # Verificar que el tipo de resultado sea sp.Matrix
#     assert isinstance(factor_func, sp.Expr)

#     # Verificar el correcto resultado
#     assert sp.simplify(true_val - factor_func) == sp.Rational(0)
    
    
# def test_modsq2mod_s_invalid_input():

#     # Verificar que se levante un ValueError 
#     with pytest.raises(ValueError):
#         factor_func = test_module.modsq2mod_s( 1.0 )


def test_isFRP_valid():

    Imm = (s**2 + 4*s + 3)/(s**2 + 2*s)    

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        test_module.isFRP(Imm)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")

    # Imm si es FRP
    assert test_module.isFRP(Imm)

    Imm = (s**2 - 4*s + 3)/(s**2 - 2*s)
    
    # Imm no es FRP
    assert not test_module.isFRP(Imm)

    
def test_isFRP_invalid_input():
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        test_module.isFRP(1.)

    

@pytest.mark.parametrize(
    "func_ptr, args, true_val",
    [
       (   
           test_module.remover_polo_sigma,
           ((s**2 + 13*s + 32)/(2*(s+1)*(s+6)), -1.),
           ((s + 8)/(2*(s + 6)), 2/(s + 1))
        ), 
       (   
           test_module.remover_polo_jw,
           ((s * (3*s**2+7) )/((s**2+1)*(s**2+3)), 1.),
           (s/(s**2 + 3), 2*s/(s**2 + 1))
        ), 
       (   
           test_module.remover_polo_dc,
           ( ((s**2+2)*(s**2+5))/(3*s*(s**2+sp.Rational(7,3))), 1.),
           ((s**2 + 1)*(s**2 + 3)/(s*(3*s**2 + 7)), 1/s)
        ), 
       (   
           test_module.remover_polo_infinito,
           ( ((s**2+2)*(s**2+5))/(3*s*(s**2+sp.Rational(7,3))), None),
           (2*(7*s**2 + 15)/(3*s*(3*s**2 + 7)), s/3)
        ), 
       (   
           test_module.remover_valor_en_infinito,
           ( (s**2 + 13*s + 32)/(3*s**2 + 27*s+ 44), None),
           (4*(3*s + 13)/(3*(3*s**2 + 27*s + 44)), 1/3)
        ), 
       (   
           test_module.remover_valor_en_dc,
           ( (3*s**2 + 27*s+ 44)/(s**2 + 13*s + 32), None),
           (s*(13*s + 73)/(8*(s**2 + 13*s + 32)), 11/8)
        ), 
     
    ]
)   
def test_remociones_valid(func_ptr, args, true_val):
    
    
    imm, sigma_omega = args

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        out = func_ptr(imm, sigma_omega)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(out[0], sp.Expr)
    assert isinstance(out[1], sp.Expr)

    # Verificar el correcto resultado
    assert sp.simplify(out[0] - true_val[0]) == sp.Rational(0)
    assert sp.simplify(out[1] - true_val[1]) == sp.Rational(0)
    
    
def test_remover_polo_sigma_invalid_input():
    
    func_ptr = test_module.remover_polo_sigma
    
    imm = (s**2 + 13*s + 32)/(2*(s+1)*(s+6))

    sigma_R1C1 = -1

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a', sigma_R1C1)

    with pytest.raises(ValueError):
        out = func_ptr(imm, 'a')

    with pytest.raises(ValueError):
        out = func_ptr(imm, sigma_R1C1, isImpedance = 2)

    with pytest.raises(ValueError):
        out = func_ptr(imm, sigma_R1C1, isRC = 2)

    with pytest.raises(ValueError):
        out = func_ptr(imm, sigma_R1C1, sigma_zero = 'a')

def test_remover_polo_jw_invalid_input():
    
    func_ptr = test_module.remover_polo_jw
    
    imm = (s * (3*s**2+7) )/((s**2+1)*(s**2+3))

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a')

    with pytest.raises(ValueError):
        out = func_ptr(imm, omega = 'a')

    with pytest.raises(ValueError):
        out = func_ptr(imm, isImpedance = 2)

    with pytest.raises(ValueError):
        out = func_ptr(imm, omega_zero = 'a')


def test_remover_polo_jw_invalid_input():
    
    func_ptr = test_module.remover_polo_jw
    
    imm = (s * (3*s**2+7) )/((s**2+1)*(s**2+3))

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a')

    with pytest.raises(ValueError):
        out = func_ptr(imm, omega = 'a')

    with pytest.raises(ValueError):
        out = func_ptr(imm, isImpedance = 2)


    with pytest.raises(ValueError):
        out = func_ptr(imm, omega_zero = 'a')

@pytest.mark.parametrize(
    "func_ptr, imm",
    [
       (   
           test_module.remover_polo_dc,
           ((s**2+2)*(s**2+5))/(3*s*(s**2+sp.Rational(7,3))),
        ), 
       (   
           test_module.remover_polo_infinito,
           ((s**2+2)*(s**2+5))/(3*s*(s**2+sp.Rational(7,3))),
        ), 
     
    ]
)   
def test_remover_polo_dc_invalid_input(func_ptr, imm):

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a')

    with pytest.raises(ValueError):
        out = func_ptr(imm, omega_zero = 'a')

    with pytest.raises(ValueError):
        out = func_ptr(imm, isSigma = 2)

@pytest.mark.parametrize(
    "func_ptr, imm",
    [
       (   
           test_module.remover_valor_en_infinito,
           (s**2 + 13*s + 32)/(3*s**2 + 27*s+ 44),
        ), 
       (   
           test_module.remover_valor_en_dc,
           (3*s**2 + 27*s+ 44)/(s**2 + 13*s + 32),
        ), 
     
    ]
)   
def test_remover_valor_invalid_input(func_ptr, imm):

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a')

    with pytest.raises(ValueError):
        out = func_ptr(imm, sigma_zero = 'a')


# def test__valid():
    
# def test__invalid_input():
