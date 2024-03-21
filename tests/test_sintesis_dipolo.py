#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
import numpy as np
import sympy as sp
from pytc2 import sintesis_dipolo as test_module


k, o = sp.symbols('k, o')
s = sp.symbols('s ', complex=True)

@pytest.mark.parametrize(
    "func_ptr, args, true_val",
    [
       (   
           test_module.cauer_RC,
           ((s**2 + 4*s + 3)/(s**2 + 2*s), True),
           (1, s/2, 4, s/6)
        ), 
       (   
           test_module.cauer_RC,
           ((s**2 + 4*s + 3)/(s**2 + 2*s), False),
           (3/(2*s), 4/5, 25/(2*s), 1/5)
        ), 
       (   
           test_module.cauer_LC,
           ((2*s**4 + 20*s**2 + 18)/(s**3 + 4*s), True),
           (2*s, s/12, 24*s/5, 5*s/36)
        ), 
       (   
           test_module.cauer_LC,
           ((2*s**4 + 20*s**2 + 18)/(s**3 + 4*s), False),
           (9/(2*s), 8/(31*s), 961/(30*s), 15/(62*s))
        ), 
     
    ]
)   
def test_cauer_valid(func_ptr, args, true_val):
    
    
    imm, sigma_omega = args

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        out = func_ptr(imm, sigma_omega)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    for result in out[0]:
        assert isinstance(result, sp.Expr) 

    # Verificar el correcto resultado
    for result, truev in zip(out[0], true_val):
        assert sp.simplify(result - truev).is_zero 
    


@pytest.mark.parametrize(
    "func_ptr",
    [
        # S -> Ts
        test_module.cauer_RC,
        test_module.cauer_LC,
    ]
)
def test_cauer_invalid_input(func_ptr):

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr(s+1/s, 'a')
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('s+1/s', remover_en_inf = True)
        
        
def test_foster_valid():
    
    
    imm = (2*s**4 + 20*s**2 + 18)/(s**3 + 4*s)

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        k0, koo, ki_wi, _, FF_foster = test_module.foster(imm)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Expr
    assert isinstance(k0, sp.Expr) 
    assert isinstance(koo, sp.Expr) 
    for result in ki_wi[0]:
        assert isinstance(result, sp.Expr) 
    
    # Verificar el correcto resultado
    assert sp.simplify(k0 - sp.Rational('9/2')).is_zero 
    assert sp.simplify(koo - sp.Rational('2')).is_zero 
    assert sp.simplify(ki_wi[0][0] - sp.Rational('8/15')).is_zero 
    assert sp.simplify(ki_wi[0][1] - sp.Rational('2/15')).is_zero 
        
def test_foster_invalid():
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        k0, koo, ki_wi, _, FF_foster = test_module.foster('imm')


def test_foster_conversion_valid():
    
    YRC = 2*(s**2 + 4*s + 3)/(s**2 + 8*s + 12)
    k0, koo, ki_wi, kk, YRC_foster = test_module.foster(YRC/s)

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        k0, koo, ki_wi, kk, YRC_foster = test_module.foster_zRC2yRC(k0, koo, ki_wi, kk, YRC_foster)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Expr
    assert isinstance(k0, type(None)) 
    assert isinstance(koo, type(None)) 
    
    for ii in ki_wi:
        for result in ii:
            assert isinstance(result, sp.Expr) 
    
    # Verificar el correcto resultado
    assert sp.simplify(ki_wi[0][0] - sp.Rational('24/5')).is_zero 
    assert sp.simplify(ki_wi[0][1] - sp.Rational('4/5')).is_zero 
    assert sp.simplify(ki_wi[1][0] - sp.Rational('8')).is_zero 
    assert sp.simplify(ki_wi[1][1] - sp.Rational('4')).is_zero 
        

# def test__valid():
    
# def test__invalid_input():
