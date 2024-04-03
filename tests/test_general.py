#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
import sympy as sp
from pytc2 import general as test_module
from numbers import Real


k, o = sp.symbols('k, o')
s = sp.symbols('s ', complex=True)
z1, z2 = sp.symbols('z1, z2')

@pytest.mark.parametrize(
    "func_ptr, args, true_val",
    [
       (   
           test_module.pp,
           (2., 2.),
           1.
        ), 
       (   
           test_module.pp,
           (z1, z2),
           z1*z2/(z1+z2)
        )     
    ]
)   
def test_pp_valid(func_ptr, args, true_val):
    
    
    zz1, zz2 = args

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        out = func_ptr(zz1, zz2)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        

    if isinstance(out, sp.Expr):

        # Verificar el correcto resultado
        assert sp.simplify(out - true_val).is_zero 
    
    else:
        
        # Verificar que el tipo de resultado sea sp.Matrix
        assert isinstance(out, Real) 
    
        # Verificar el correcto resultado
        assert (out - true_val) == 0


def test_cauer_invalid_input():

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = test_module.pp(z1, 'a')
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = test_module.pp('a', z2)
        
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = test_module.pp(2., 'a')
        
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = test_module.pp('a', 2.)


def test_simplify_n_monic_valid():
    
    tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
    truev = (s + 2)/(2*s + 3)
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        simplified_tt = test_module.simplify_n_monic(tt)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(simplified_tt, sp.Expr) 

    # Verificar el correcto resultado
    assert sp.simplify(simplified_tt - truev).is_zero 
    
    
def test_simplify_n_monic_invalid_input():
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        simplified_tt = test_module.simplify_n_monic('tt')
    
    
w = test_module.w    

@pytest.mark.parametrize(
    "args, true_val",
    [
       (   
           0,
           sp.Rational(1)
        ), 
       (   
           1,
           w
        ),    
       (   
           2,
           2*w**2-1
        ),     
       (   
           3,
           4*w**3-3*w
        )     
    ]
)    
def test_Chebyshev_polynomials_valid(args, true_val):
    
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        simplified_tt = test_module.Chebyshev_polynomials(args)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(true_val, sp.Expr) 

    # Verificar el correcto resultado
    assert sp.simplify(true_val - true_val).is_zero 
    
    
def test_Chebyshev_polynomials_invalid_input():
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        simplified_tt = test_module.Chebyshev_polynomials('tt')
    
    with pytest.raises(ValueError):
        simplified_tt = test_module.Chebyshev_polynomials(2.5)
    

def test_a_equal_b_latex_s_valid():
    
    tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
    truev = '$z_{1}=\\frac{s^{2} + 3 s + 2}{2 s^{2} + 5 s + 3}$'
    truev2 = '$z_{1}=\\left[ \\frac{s^{2} + 3 s + 2}{2 s^{2} + 5 s + 3}, \\  z_{1}\\right]$'
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        outv = test_module.a_equal_b_latex_s(z1, tt)
        outv1 = test_module.a_equal_b_latex_s('z_{1}', tt)
        outv2 = test_module.a_equal_b_latex_s('z_{1}', [tt, z1])
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(outv, str) 

    # Verificar el correcto resultado
    assert outv == truev
    assert outv1 == truev
    assert outv2 == truev2
    
    
def test_a_equal_b_latex_s_invalid_input():
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        simplified_tt = test_module.a_equal_b_latex_s('tt', 1.)

    with pytest.raises(ValueError):
        simplified_tt = test_module.a_equal_b_latex_s(['tt', '1.'], '1.')

    with pytest.raises(ValueError):
        simplified_tt = test_module.a_equal_b_latex_s( 1., ['tt', '1.'])

    with pytest.raises(ValueError):
        simplified_tt = test_module.a_equal_b_latex_s(1., z1)

def test_expr_simb_expr_valid():
    
    tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
    truev = '$z_{1}\\neq\\frac{s^{2} + 3 s + 2}{2 s^{2} + 5 s + 3}$'
    truev2 = '$z_{1}\\neq\\left[ \\frac{s^{2} + 3 s + 2}{2 s^{2} + 5 s + 3}, \\  z_{1}\\right]$'
    asym = '\\neq'
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        outv = test_module.expr_simb_expr(z1, tt, symbol = asym)
        outv1 = test_module.expr_simb_expr('z_{1}', tt, symbol = asym)
        outv2 = test_module.expr_simb_expr('z_{1}', [tt, z1], symbol = asym)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(outv, str) 

    # Verificar el correcto resultado
    assert outv == truev
    assert outv1 == truev
    assert outv2 == truev2
    
    
def test_expr_simb_expr_invalid_input():
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        simplified_tt = test_module.expr_simb_expr('tt', 1.)

    with pytest.raises(ValueError):
        simplified_tt = test_module.expr_simb_expr('tt', '1.')

    with pytest.raises(ValueError):
        simplified_tt = test_module.expr_simb_expr( 1., ['tt', '1.'])

    with pytest.raises(ValueError):
        simplified_tt = test_module.expr_simb_expr(1., z1)

    with pytest.raises(ValueError):
        simplified_tt = test_module.expr_simb_expr(1., z1, symbol = 1.)

def test_to_latex_valid():
    
    tt = '\\alpha'
    truev = '$\\alpha$'
    truev1 = '$z_{1}$'
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        outv = test_module.to_latex(tt)
        outv1 = test_module.to_latex(z1)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(outv, str) 
    assert isinstance(outv1, str) 

    # Verificar el correcto resultado
    assert outv == truev
    assert outv1 == truev1
    
    
def test_db2nepper_invalid_input():

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        outv = test_module.to_latex(1.)

def test_db2nepper_latex_valid():
    
    at_db = 8.685889638065037
    truev = 1.
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        outv = test_module.nepper2db(truev)
        outv1 = test_module.db2nepper(at_db)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(outv, Real) 
    assert isinstance(outv1, Real) 

    # Verificar el correcto resultado
    assert outv == at_db
    assert outv1 == truev
    
    
def test_str_to_latex_invalid_input():

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        outv = test_module.nepper2db('1.')

    with pytest.raises(ValueError):
        outv = test_module.db2nepper('1.')

# def test__valid():
    
# def test__invalid_input():
