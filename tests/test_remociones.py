#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
import numpy as np
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



# def test__valid():
    
# def test__invalid_input():
