#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
import numpy as np
import sympy as sp
import scipy.signal as sig
from pytc2 import cuadripolos as test_module

# símbolos genéricos 
S11, S12, S21, S22 = sp.symbols('S11, S12, S21, S22')


@pytest.mark.parametrize(
    "func_ptr, true_val",
    [
       # S -> Ts
       (   
           test_module.S2Ts_s,
           sp.Matrix([[1/S21, -S22/S21], [S11/S21, -S11*S22/S21 + S12]])
        ), 
       # Ts -> S 
       (   
           test_module.Ts2S_s,
           sp.Matrix([[S21/S11, S22 - S12*S21/S11], [1/S11, -S12/S11]])
        ), 
       # Ts -> Tabcd 
       (   
           test_module.Ts2Tabcd_s,
           sp.Matrix([[S11/2 + S12/2 + S21/2 + S22/2, S11/2 - S12/2 + S21/2 - S22/2], [S11/2 + S12/2 - S21/2 - S22/2, S11/2 - S12/2 - S21/2 + S22/2]])
        ), 
       # Tabcd -> S 
       (   
           test_module.Tabcd2S_s,
           sp.Matrix([[(S11 + S12 - S21 - S22)/(S11 + S12 + S21 + S22), 2*(S11*S22 - S12*S21)/(S11 + S12 + S21 + S22)], [2/(S11 + S12 + S21 + S22), (-S11 + S12 - S21 + S22)/(S11 + S12 + S21 + S22)]])
        ), 
       # S -> Tabcd 
       (   
           test_module.S2Tabcd_s,
           sp.Matrix([[(-S11*S22 + S11 + S12*S21 - S22 + 1)/(2*S21), (S11*S22 + S11 - S12*S21 + S22 + 1)/(2*S21)], [(S11*S22 - S11 - S12*S21 - S22 + 1)/(2*S21), (-S11*S22 - S11 + S12*S21 + S22 + 1)/(2*S21)]])
        ), 
       # Y -> Tabcd 
       (   
           test_module.Y2Tabcd_s,
           sp.Matrix([[-S22/S21, -1/S21], [-(S11*S22 - S12*S21)/S21, -S22/S21]])
        ), 
       # Tabcd -> Y
       (   
           test_module.Tabcd2Y_s,
           sp.Matrix([[S22/S12, -(S11*S22 - S12*S21)/S12], [-1/S12, S11/S12]])
        ), 
       # Z -> Tabcd 
       (   
           test_module.Z2Tabcd_s,
           sp.Matrix([[S11/S21, (S11*S22 - S12*S21)/S21], [1/S21, S22/S21]])
        ), 
       # Tabcd -> Z
       (   
           test_module.Tabcd2Z_s,
           sp.Matrix([[S11/S21, (S11*S22 - S12*S21)/S21], [1/S21, S22/S21]])
        ), 
      
    ]
)
def test_matrix_conversion_simbolic(func_ptr, true_val):
    
    Spar = sp.Matrix([[S11, S12],
                      [S21, S22]])

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Ts = func_ptr(Spar)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Ts, sp.MatrixBase)

    # Verificar el correcto resultado
    assert sp.simplify(Ts - true_val) == sp.Matrix([[0, 0], [0, 0]])


@pytest.mark.parametrize(
    "func_ptr, Spar_invalid",
    [
       # S -> Ts
       (   
           test_module.S2Ts_s,
           sp.Matrix([[S11, S12], [0, S22]])
        ), 
       # Ts -> S 
       (   
           test_module.Ts2S_s,
           sp.Matrix([[0, S12], [S21, S22]])
        ), 
       # Y -> Tabcd 
       (   
           test_module.Y2Tabcd_s,
           sp.Matrix([[S11, S12], [0, S22]])
        ), 
       # Tabcd -> Y
       (   
           test_module.Tabcd2Y_s,
           sp.Matrix([[S11, 0], [S21, S22]])
        ), 
       # Z -> Tabcd 
       (   
           test_module.Z2Tabcd_s,
           sp.Matrix([[S11, S12], [0, S22]])
        ), 
       # Tabcd -> Z
       (   
           test_module.Tabcd2Z_s,
           sp.Matrix([[S11, S12], [0, S22]])
        ), 
      
    ]
)
def test_matrix_conversion_simbolic_invalid_input(func_ptr, Spar_invalid):
    
    Spar_mal_tipo = np.random.randn(2,2)
    
    Spar_mala_dim = sp.Matrix([[S11, S12, S22],
                               [S11, S12, S22],
                               [S11, S12, S22]])

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(Spar_mal_tipo)
    
    # Verificar que se levante un ValueError al pasar una matriz mal dimensionada
    with pytest.raises(ValueError):
        _ = func_ptr(Spar_mala_dim)
    
    # Verificar que se levante un ValueError al pasar una matriz no válida
    with pytest.raises(ValueError):
        _ = func_ptr(Spar_invalid)
    

@pytest.mark.parametrize(
    "func_ptr, Spar_invalid",
    [
       # Ts -> Tabcd 
       (   
           test_module.Ts2Tabcd_s,
           sp.Matrix([[0, S12], [S21, S22]])
        ), 
       # Tabcd -> S 
       (   
           test_module.Tabcd2S_s,
           sp.Matrix([[S11, S12, S22], [S11, S12, S22], [S11, S12, S22]])
        ), 
       # S -> Tabcd 
       (   
           test_module.S2Tabcd_s,
           sp.Matrix([[S11, S12], [0, S22]])
        ), 
    ]
)
def test_matrix_conversion_simbolic_invalid_input_2(func_ptr, Spar_invalid):
    
    Spar_mal_tipo = np.random.randn(2,2)
    
    Spar_mal_tipo2 = np.random.randn(1)
    
    Spar_mala_dim = sp.Matrix([[S11, S12, S22],
                               [S11, S12, S22],
                               [S11, S12, S22]])
    
    Spar_valida = sp.Matrix([[S11, S12],
                             [S21, S22]])

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        _ = func_ptr(Spar_mal_tipo)
    
    # Verificar que se levante un ValueError al pasar un Z01 numérica
    with pytest.raises(ValueError):
        _ = func_ptr(Spar_valida, Z01=Spar_mal_tipo2)
    
    # Verificar que se levante un ValueError al pasar un Z02 numérica
    with pytest.raises(ValueError):
        _ = func_ptr(Spar_valida, Z02=Spar_mal_tipo2)
    
    # Verificar que se levante un ValueError al pasar una matriz mal dimensionada
    with pytest.raises(ValueError):
        _ = func_ptr(Spar_mala_dim)
    
    # Verificar que se levante un ValueError al pasar una matriz no válida
    with pytest.raises(ValueError):
        _ = func_ptr(Spar_invalid)
    

def test_I2Tabcd_s_valid():
    
    true_val = sp.Matrix([[sp.sqrt(S12/S21)*sp.cosh(S11), sp.sqrt(S12*S21)*sp.sinh(S11)], [sp.sinh(S11)/sp.sqrt(S12*S21), sp.sqrt(S21/S12)*sp.cosh(S11)]])

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Ts = test_module.I2Tabcd_s(S11, S12, S21)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Ts, sp.MatrixBase)

    # Verificar el correcto resultado
    assert sp.simplify(Ts - true_val) == sp.Matrix([[0, 0], [0, 0]])

def test_I2Tabcd_s_invalid_input():

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        Ts = test_module.I2Tabcd_s(1, S12, S21)
    
    # Verificar que se levante un ValueError al pasar una matriz mal dimensionada
    with pytest.raises(ValueError):
        Ts = test_module.I2Tabcd_s(S12, 1, S21)
    
    # Verificar que se levante un ValueError al pasar una matriz no válida
    with pytest.raises(ValueError):
        Ts = test_module.I2Tabcd_s(S12, S21, 1)

    

y11, y12, y21, y22 = sp.symbols('y11, y12, y21, y22', complex=True)
z11, z12, z21, z22 = sp.symbols('z11, z12, z21, z22', complex=True)
A, B, C, D = sp.symbols('A, B, C, D', complex=True)
Ai, Bi, Ci, Di = sp.symbols('Ai, Bi, Ci, Di', complex=True)
h11, h12, h21, h22 = sp.symbols('h11, h12, h21, h22', complex=True)
g11, g12, g21, g22 = sp.symbols('g11, g12, g21, g22', complex=True)
v1, v2, i1, i2 = sp.symbols('v1, v2, i1, i2', complex=True)

# Parámetros Z (impedancia - circ. abierto)
ZZ = sp.Matrix([[z11, z12], [z21, z22]])
# vars. dependientes
vv = sp.Matrix([[v1], [v2]])
# vars. INdependientes
ii = sp.Matrix([[i1], [i2]])

# Parámetros Y (admitancia - corto circ.)
YY = sp.Matrix([[y11, y12], [y21, y22]])
# vars. dependientes
# ii = sp.Matrix([[i1], [i2]])
# vars. INdependientes
# vv = sp.Matrix([[v1], [v2]])

# Parámetros H (híbridos h)
HH = sp.Matrix([[h11, h12], [h21, h22]])
# vars. dependientes
h_dep = sp.Matrix([[v1], [i2]])
# vars. INdependientes
h_ind = sp.Matrix([[i1], [v2]])

# Parámetros G (híbridos g)
GG = sp.Matrix([[g11, g12], [g21, g22]])
# vars. dependientes
g_dep = sp.Matrix([[i1], [v2]])
# vars. INdependientes
g_ind = sp.Matrix([[v1], [i2]])

# Parámetros Tabcd (Transmisión, ABCD)
TT = sp.Matrix([[A, -B], [C, -D]])
# vars. dependientes
t_dep = sp.Matrix([[v1], [i1]])
# vars. INdependientes.  (Signo negativo de corriente)
t_ind = sp.Matrix([[v2], [i2]])

# Parámetros Tdcba (Transmisión inversos, DCBA)
TTi = sp.Matrix([[Ai, Bi], [-Ci, -Di]])
# vars. dependientes
ti_dep = sp.Matrix([[v2], [i2]])
# vars. INdependientes. (Signo negativo de corriente)
ti_ind = sp.Matrix([[v1], [i1]])

# Diccionario con la definición de cada modelo
model_dct = [ { 'model_name': 'Z', 'matrix': ZZ, 'dep_var': vv, 'indep_var':ii },
              { 'model_name': 'Y', 'matrix': YY, 'dep_var': ii, 'indep_var':vv },
              { 'model_name': 'H', 'matrix': HH, 'dep_var': h_dep, 'indep_var':h_ind },
              { 'model_name': 'G', 'matrix': GG, 'dep_var': g_dep, 'indep_var':g_ind },
              { 'model_name': 'T', 'matrix': TT, 'dep_var': t_dep, 'indep_var':t_ind, 'neg_i2_current': True },
              { 'model_name': 'Ti', 'matrix': TTi, 'dep_var': ti_dep, 'indep_var':ti_ind, 'neg_i2_current': True}
            ]


def test_Model_conversion_valid():
   
    for dst_model in model_dct:
        
        for src_model in model_dct:
            
            try:
                HH_z = test_module.Model_conversion( src_model, dst_model )
            except ValueError:
                pytest.fail("Se levantó un ValueError incorrectamente.")

malos_indices = np.random.randn(2)

# diccionarios inválidos
invalid_model_dct = [ {  'matrix': ZZ, 'dep_var': vv, 'indep_var':malos_indices },
                      { 'model_name': 'W', 'dep_var': ii, 'indep_var':malos_indices },
                      { 'model_name': 'W', 'matrix': ZZ, 'dep_var': vv, 'indep_var':malos_indices },                      
                      { 'model_name': 'W', 'matrix': ZZ, 'dep_var': vv, 'indep_var':malos_indices },                      
                      { 'model_name': 'H', 'matrix': HH, 'indep_var':h_ind },
                      { 'model_name': 'G', 'matrix': GG, 'dep_var': g_dep},
                      { 'dep_var': t_dep, 'indep_var':t_ind, 'neg_i2_current': True },
                      { 'model_name': 'Ti', 'matrix': TTi, 'dep_var': ti_dep}
                    ]
    
def test_Model_conversion_invalid_input():
    
    for dst_model in model_dct:
        
        for src_model in invalid_model_dct:
            
            with pytest.raises(ValueError):
                
                HH_z = test_module.Model_conversion( src_model, dst_model)
            
    for dst_model in invalid_model_dct:
        
        for src_model in model_dct:
            
            with pytest.raises(ValueError):
                
                HH_z = test_module.Model_conversion( src_model, dst_model)
            

def test_y2mai_valid():


    true_val = sp.Matrix([[y11, y12, -y11 - y12], [y21, y22, -y21 - y22], [-y11 - y21, -y12 - y22, y11 + y12 + y21 + y22]])

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Ts = test_module.y2mai(YY)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Ts, sp.MatrixBase)

    # Verificar el correcto resultado
    assert sp.simplify(Ts - true_val) == sp.Matrix([[0, 0, 0],[0, 0, 0],[0, 0, 0]])

    
def test_y2mai_invalid_input():

    with pytest.raises(ValueError):
        HH_z = test_module.y2mai(1)
    
    with pytest.raises(ValueError):
        HH_z = test_module.y2mai(np.random.randn(2,2))
    
    
def test_may2y_valid():

    true_val = sp.Matrix([[y11, y12], [y21, y22]])
    Ymai = sp.Matrix([[y11, y12, -y11 - y12], [y21, y22, -y21 - y22], [-y11 - y21, -y12 - y22, y11 + y12 + y21 + y22]])

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Ts = test_module.may2y(Ymai, 2)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Ts, sp.MatrixBase)

    # Verificar el correcto resultado
    assert sp.simplify(Ts - true_val) == sp.Matrix([[0, 0],[0, 0]])
    
    
def test_may2y_invalid_input():
    
    Ymai = sp.Matrix([[y11, y12, -y11 - y12], [y21, y22, -y21 - y22], [-y11 - y21, -y12 - y22, y11 + y12 + y21 + y22]])
    
    with pytest.raises(ValueError):
        YY = test_module.may2y(1, 2)
    
    with pytest.raises(ValueError):
        YY = test_module.may2y(np.random.randn(3,3), 3)

    with pytest.raises(ValueError):
        YY = test_module.may2y(Ymai, 3)
    
    with pytest.raises(ValueError):
        YY = test_module.may2y(Ymai, [2, 'a'])

    with pytest.raises(ValueError):
        YY = test_module.may2y(Ymai, 'a')


@pytest.mark.parametrize(
    "func_ptr, true_val",
    [
       # Y -> Tabcd
       (   
           test_module.Y2Tabcd,
           np.array([[-1.33333333, -0.33333333], [ 0.66666667, -0.33333333]])
        ), 
       # Tabcd -> Y
       (   
           test_module.Tabcd2Y,
           np.array([[ 2.,   1. ], [-0.5, 0.5]])
        ),
       # Z -> Tabcd
       (   
           test_module.Z2Tabcd,
           np.array([[ 0.33333333, -0.66666667], [ 0.33333333,  1.33333333]])
        ), 
       # Tabcd -> Z
       (   
           test_module.Tabcd2Z,
           np.array([[ 0.33333333, -0.66666667], [ 0.33333333,  1.33333333]])
        ), 
      
    ]
)
def test_matrix_conversion_numeric(func_ptr, true_val):
    
    matrix = np.array([[1., 2.],
                       [3., 4.]])

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Ts = func_ptr(matrix)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Ts, np.ndarray)

    # Verificar el correcto resultado
    assert np.median(Ts - true_val) < 1e-6


@pytest.mark.parametrize(
    "func_ptr, Spar_invalid",
    [
     
       # Y -> Tabcd
       (   
           test_module.Y2Tabcd,
           np.array([[-1.33333333, -0.33333333], [ 0., -0.33333333]])
        ), 
       # Z -> Tabcd
       (   
           test_module.Z2Tabcd,
           np.array([[ 0.33333333, -0.66666667], [ 0.,  1.33333333]])
        ), 
       # Tabcd -> Y
       (   
           test_module.Tabcd2Y,
           np.array([[ 2.,   0. ], [-0.5, 0.5]])
        ),
       # Tabcd -> Z
       (   
           test_module.Tabcd2Z,
           np.array([[ 0.33333333, -0.66666667], [ 0.,  1.33333333]])
        ) 

    ]
)
def test_matrix_conversion_numeric_invalid_input(func_ptr, Spar_invalid):
    
    Spar_mala_dim = np.random.randn(3,3)
    
    Spar_mal_tipo = sp.Matrix([[S11, S12, S22],
                               [S11, S12, S22],
                               [S11, S12, S22]])

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        Ts = func_ptr(Spar_mal_tipo)
    
    # Verificar que se levante un ValueError al pasar una matriz mal dimensionada
    with pytest.raises(ValueError):
        Ts = func_ptr(Spar_mala_dim)
    
    # Verificar que se levante un ValueError al pasar una matriz no válida
    with pytest.raises(ValueError):
        Ts = func_ptr(Spar_invalid)
    

def test_I2Tabcd_valid():
    
    true_val = np.array([[0.68073771+0.8074316j,  1.5553376 +3.18055853j],
                         [0.25922293+0.53009309j, 1.02110657+1.21114739j]])

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Ts = test_module.I2Tabcd(1.+1.j, 2, 3)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Ts, np.ndarray)

    # Verificar el correcto resultado
    assert np.median(Ts - true_val) < 1e-6

    
def test_I2Tabcd_invalid_input():

    # Verificar que se levante un ValueError al pasar una matriz numérica como gamma
    with pytest.raises(ValueError):
        Ts = test_module.I2Tabcd( np.array([1.0, 2.0]), 2, 3)
    
    # Verificar que se levante un ValueError al pasar un nivel de impedancia negativo
    with pytest.raises(ValueError):
        Ts = test_module.I2Tabcd( 1., -2., 3)
    
    # Verificar que se levante un ValueError al pasar un nivel de impedancia negativo
    with pytest.raises(ValueError):
        Ts = test_module.I2Tabcd( 1., 2., -3)

Zexc = sp.symbols('Zexc')
Z01 = sp.symbols('Z01')
Z02 = sp.symbols('Z02')
Y01 = sp.symbols('Y01')
Y02 = sp.symbols('Y02')
Yexc = sp.symbols('Yexc')

@pytest.mark.parametrize(
    "func_ptr, f_exc, p1, p2, true_val",
    [
       # Z_s
       (   
           test_module.SparZ_s,
           sp.symbols('Zexc'),
           sp.symbols('Z01'),
           sp.symbols('Z02'),
           sp.Matrix([[(-Z01 + Z02 + Zexc)/(Z01 + Z02 + Zexc), 2*Z01*sp.sqrt(Z02/Z01)/(Z01 + Z02 + Zexc)], [2*Z02*sp.sqrt(Z01/Z02)/(Z01 + Z02 + Zexc), (Z01 - Z02 + Zexc)/(Z01 + Z02 + Zexc)]])
        ), 
       (   
           test_module.SparZ_s,
           sp.symbols('Zexc'),
           sp.symbols('Z01'),
           sp.symbols('Z01'),
           sp.Matrix([[Zexc/(2*Z01 + Zexc), 2*Z01/(2*Z01 + Zexc)], [2*Z01/(2*Z01 + Zexc), Zexc/(2*Z01 + Zexc)]])
        ), 
       # Y_s
       (   
           test_module.SparY_s,
           sp.symbols('Yexc'),
           sp.symbols('Y01'),
           sp.symbols('Y02'),
           sp.Matrix([[(Y01 - Y02 - Yexc)/(Y01 + Y02 + Yexc), 2*Y01*sp.sqrt(Y01/Y02)/(Y01 + Y02 + Yexc)], [2*Y02*sp.sqrt(Y02/Y01)/(Y01 + Y02 + Yexc), (-Y01 + Y02 - Yexc)/(Y01 + Y02 + Yexc)]])
        ), 
       (   
           test_module.SparY_s,
           sp.symbols('Yexc'),
           sp.symbols('Y01'),
           sp.symbols('Y01'),
           sp.Matrix([[-Yexc/(2*Y01 + Yexc), 2*Y01/(2*Y01 + Yexc)], [2*Y01/(2*Y01 + Yexc), -Yexc/(2*Y01 + Yexc)]])
        ), 
      
    ]
)
def test_matrix_def_simbolic(func_ptr, f_exc, p1, p2, true_val):
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Ts = func_ptr(f_exc, p1, p2)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Ts, sp.MatrixBase)

    # Verificar el correcto resultado
    assert sp.simplify(Ts - true_val) == sp.Matrix([[0, 0], [0, 0]])


@pytest.mark.parametrize(
    "func_ptr",
    [
       # Z_s
           test_module.SparZ_s,
       # Y_s
           test_module.SparY_s,
      
    ]
)
def test_matrix_def_simbolic_invalid_input(func_ptr):
    
    Spar_mal_tipo = np.random.randn(2,2)

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        Ts = func_ptr(Spar_mal_tipo, y11)
    
    # Verificar que se levante un ValueError al pasar un float
    with pytest.raises(ValueError):
        Ts = func_ptr(y11, y11, 1.)
    

@pytest.mark.parametrize(
    "func_ptr, f_exc1, f_exc2, true_val",
    [
       # LYZ_s
       (   
           test_module.TabcdLYZ_s,
           sp.symbols('Yexc'),
           sp.symbols('Zexc'),
           sp.Matrix([[1, Zexc], [Yexc, Yexc*Zexc + 1]])
        ), 
       (   
           test_module.TabcdLZY_s,
           sp.symbols('Zexc'),
           sp.symbols('Yexc'),
           sp.Matrix([[Yexc*Zexc + 1, Zexc], [Yexc, 1]])
        ), 
      
    ]
)
def test_matrix_def1_simbolic(func_ptr, f_exc1, f_exc2, true_val):
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Ts = func_ptr(f_exc1, f_exc2)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Ts, sp.MatrixBase)

    # Verificar el correcto resultado
    assert sp.simplify(Ts - true_val) == sp.Matrix([[0, 0], [0, 0]])


@pytest.mark.parametrize(
    "func_ptr",
    [
       # Z_s
           test_module.TabcdLYZ_s,
       # Y_s
           test_module.TabcdLZY_s,
      
    ]
)
def test_matrix_def1_simbolic_invalid_input(func_ptr):

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        Ts = func_ptr(1., y11)
    
    # Verificar que se levante un ValueError al pasar un float
    with pytest.raises(ValueError):
        Ts = func_ptr(y11, 1.)
    

@pytest.mark.parametrize(
    "func_ptr, f_exc, true_val",
    [
       # LYZ_s
       (   
           test_module.TabcdZ_s,
           sp.symbols('Zexc'),
           sp.Matrix([[1, Zexc], [0, 1]])
        ), 
       (   
           test_module.TabcdY_s,
           sp.symbols('Yexc'),
           sp.Matrix([[1, 0], [Yexc, 1]])
        ), 
      
    ]
)
def test_matrix_def2_simbolic(func_ptr, f_exc, true_val):
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Ts = func_ptr(f_exc)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Ts, sp.MatrixBase)

    # Verificar el correcto resultado
    assert sp.simplify(Ts - true_val) == sp.Matrix([[0, 0], [0, 0]])


@pytest.mark.parametrize(
    "func_ptr",
    [
       # Z_s
           test_module.TabcdZ_s,
       # Y_s
           test_module.TabcdY_s,
      
    ]
)
def test_matrix_def2_simbolic_invalid_input(func_ptr):

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        Ts = func_ptr(1.)

Y1, Y2, Y3 = sp.symbols('Y1 Y2 Y3', complex=True)
G = sp.symbols('G', real=True, positive=True)
s = sp.symbols('s ', complex=True)

@pytest.mark.parametrize(
    "func_ptr, true_val",
    [
         (
             test_module.calc_MAI_vtransf_ij_mn,
             -1/(2*G*s + 2*s**2*(G*s + 1) + 1)
         ),
         (
             test_module.calc_MAI_ztransf_ij_mn,
             -1/(2*G*s**2 + G + 2*s)
         )
    ]
)
def test_MAI_transf( func_ptr, true_val):

    input_port = [0, 1]
    output_port = [3, 1]
    #      Nodos: 0      1        2        3
    Ymai = sp.Matrix([  
                    [ Y1,    0,      -Y1,      0],
                    [ 0,    Y2+G,    -Y2,     -G],
                    [ -Y1,  -Y2,    Y1+Y2+Y3, -Y3],
                    [ 0,    -G,      -Y3,      Y3+G ]
                    ])
    # Butter de 3er orden doblemente cargado
    Ymai = Ymai.subs(Y1, 1/s/sp.Rational('1'))
    Ymai = Ymai.subs(Y3, 1/s/sp.Rational('1'))
    Ymai = Ymai.subs(Y2, s*sp.Rational('2'))

    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Zmai = func_ptr(Ymai, output_port[0], output_port[1], input_port[0], input_port[1], verbose=False)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Zmai, sp.Expr)

    # Verificar el correcto resultado
    assert sp.simplify(Zmai - true_val) == sp.Rational('0')


@pytest.mark.parametrize(
    "func_ptr",
    [
     test_module.calc_MAI_vtransf_ij_mn,
     test_module.calc_MAI_ztransf_ij_mn      
    ]
)
def test_MAI_transf_invalid_input(func_ptr):

    Ymai_mal_tipo = np.random.randn(2,2)
    
    Ymai = sp.Matrix([[S11, S12, S22],
                               [S11, S12, S22],
                               [S11, S12, S22]])

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        Zmai = func_ptr(Ymai_mal_tipo, 0, 1, 3, 1, verbose=False)

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        Zmai = func_ptr(Ymai, 0.5, 1, 3, 1, verbose=False)
    
    with pytest.raises(ValueError):
        Zmai = func_ptr(Ymai, 0, 0.5, 3, 1, verbose=False)
    
    with pytest.raises(ValueError):
        Zmai = func_ptr(Ymai, 0, 1, 0.5, 1, verbose=False)
    
    with pytest.raises(ValueError):
        Zmai = func_ptr(Ymai, 0, 1, 3, 0.5, verbose=False)
    
    with pytest.raises(ValueError):
        Zmai = func_ptr(Ymai, 0, 1, 3, 1, verbose=1.0)
    

    
def test_MAI_impedance_valid():
    
    input_port = [0, 1]
    #      Nodos: 0      1        2        3
    Ymai = sp.Matrix([  
                    [ Y1,    0,      -Y1,      0],
                    [ 0,    Y2+G,    -Y2,     -G],
                    [ -Y1,  -Y2,    Y1+Y2+Y3, -Y3],
                    [ 0,    -G,      -Y3,      Y3+G ]
                    ])
    # Butter de 3er orden doblemente cargado
    Ymai = Ymai.subs(Y1, 1/s/sp.Rational('1'))
    Ymai = Ymai.subs(Y3, 1/s/sp.Rational('1'))
    Ymai = Ymai.subs(Y2, s*sp.Rational('2'))

    true_val = (2*G*s + 2*s**2*(G*s + 1) + 1)/(2*G*s**2 + G + 2*s)
    
    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        Zmai = test_module.calc_MAI_impedance_ij(Ymai, input_port[0], input_port[1], verbose=False)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")
        
    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(Zmai, sp.Expr)

    # Verificar el correcto resultado
    assert sp.simplify(Zmai - true_val) == sp.Rational('0')
    
    
def test_MAI_impedance_invalid_input():
    
    Ymai_mal_tipo = np.random.randn(2,2)
    
    Ymai = sp.Matrix([[S11, S12, S22],
                               [S11, S12, S22],
                               [S11, S12, S22]])

    func_ptr = test_module.calc_MAI_impedance_ij

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        Zmai = func_ptr(Ymai_mal_tipo, 0, 1, verbose=False)

    # Verificar que se levante un ValueError al pasar una matriz numérica
    with pytest.raises(ValueError):
        Zmai = func_ptr(Ymai, 0.5, 1, verbose=False)
    
    with pytest.raises(ValueError):
        Zmai = func_ptr(Ymai, 0, 0.5, verbose=False)
    
    with pytest.raises(ValueError):
        Zmai = func_ptr(Ymai, 0, 1, verbose=1.0)
    
    

# def test__valid():
    
# def test__invalid_input():

