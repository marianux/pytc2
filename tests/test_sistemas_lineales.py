#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
from scipy.signal import TransferFunction
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from pytc2 import sistemas_lineales as test_module
import scipy.signal as sig

def test_tfcascade():
    # Definir las funciones de transferencia de prueba
    tfa = TransferFunction([1, 2], [1, 3, 2])
    tfb = TransferFunction([1], [1, 4])

    # Calcular el resultado esperado de la cascada
    expected_num = [1, 2]
    expected_den = [1, 7, 14, 8]

    # Llamar a la función tfcascade con las funciones de transferencia de prueba
    result = test_module.tfcascade(tfa, tfb)

    # Verificar que el tipo de resultado sea TransferFunction
    assert isinstance(result, TransferFunction)

    # Verificar que los coeficientes del numerador y el denominador sean los esperados
    assert np.array_equal(result.num, expected_num)
    assert np.array_equal(result.den, expected_den)

def test_tfcascade_invalid_input():
    # Definir una función de transferencia incorrecta
    tfa = TransferFunction([1, 2], [1, 3, 2])
    invalid_tfb = "not_a_transfer_function"

    # Verificar que se levante un ValueError al pasar una función de transferencia no válida
    with pytest.raises(ValueError):
        test_module.tfcascade(tfa, invalid_tfb)

    # Definir funciones de transferencia válidas
    valid_tfa = TransferFunction([1, 2], [1, 3, 2])
    valid_tfb = TransferFunction([1], [1, 4])

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        test_module.tfcascade(valid_tfa, valid_tfb)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")

def test_tfadd_valid_input():
    # Definir funciones de transferencia de prueba
    tfa = TransferFunction([1, 2], [3, 4])
    tfb = TransferFunction([5, 6], [7, 8])

    # Calcular el resultado esperado de la suma
    expected_num = np.polyadd(np.polymul(tfa.num, tfb.den), np.polymul(tfa.den, tfb.num))
    expected_den = np.polymul(tfa.den, tfb.den)

    # Llamar a la función tfadd con las funciones de transferencia de prueba
    result = test_module.tfadd(tfa, tfb)

    # Verificar que el tipo de resultado sea TransferFunction
    assert isinstance(result, TransferFunction)

    # Verificar que los coeficientes del numerador y del denominador sean los esperados
    assert np.array_equal(result.num, expected_num)
    assert np.array_equal(result.den, expected_den)

def test_tfadd_invalid_input():
    # Definir una función de transferencia incorrecta
    tfa = TransferFunction([1, 2], [3, 4])
    invalid_tfb = "not_a_transfer_function"

    # Verificar que se levante un ValueError al pasar una función de transferencia no válida
    with pytest.raises(ValueError):
        test_module.tfadd(tfa, invalid_tfb)

    # Definir funciones de transferencia válidas
    valid_tfa = TransferFunction([1, 2], [3, 4])
    valid_tfb = TransferFunction([5, 6], [7, 8])

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        test_module.tfadd(valid_tfa, valid_tfb)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")

# varias SOS que prueben muchas posibilidades de sistemas válidos
@pytest.mark.parametrize(
    "mySOS",
    [
        np.array([ [1., 7., 1., 1., 9., 1.],
                   [1., 2., 1., 1., 3., 1.]]),
        np.array([ [0., 1., 1., 0., 1., 2.],
                   [1., 3./7., 9., 1., 4./5., 16.],                  
                   [1., 1./3., 1., 1., 5./2., 25.]]), 
        np.array([ [1., 1., 0., 1., 1., 0.],
                   [1., 3./7., 9., 1., 4./5., 16.],                  
                   [1., 1./3., 1., 1., 5./2., 25.]]), 
        np.array([ [1., 0., 0., 0., 1., 1.],
                   [0., 0., 1., 1., 4., 3.],                  
                   [1., 7./1., 10., 1., 1./5., 1.]]), 
    ])
def test_sos2tf_analog_valid_input(mySOS):
    # Definir una matriz SOS de prueba

    # tolerancia numérica
    tol = 1e-10

    # Calcular la función de transferencia analógica esperada
    # Calcular la matriz SOS esperada de otra manera
    sos_num, sos_den = test_module._one_sos2tf(mySOS[0,:])
    expected_tf = TransferFunction(sos_num, sos_den)
    for ii in range(1, mySOS.shape[0]):
        sos_num, sos_den = test_module._one_sos2tf(mySOS[ii,:])
        tfa = TransferFunction(sos_num, sos_den)
        expected_tf = test_module.tfcascade(expected_tf, tfa)

    # Llamar a la función sos2tf_analog con la matriz SOS de prueba
    result_tf = test_module.sos2tf_analog(mySOS)

    # Verificar que el tipo de resultado sea TransferFunction
    assert isinstance(result_tf, TransferFunction)

    # Verificar que los coeficientes de la función de transferencia sean los esperados
    assert np.max(np.abs( result_tf.num - expected_tf.num)) < tol
    assert np.max(np.abs( result_tf.den - expected_tf.den)) < tol

def test_sos2tf_analog_invalid_input():
    # Definir una matriz SOS no válida (no 2D)
    invalid_SOS = [1, 0.5, 1, 1, 0.2, 1]

    # Verificar que se levante un ValueError al pasar una matriz SOS no válida
    with pytest.raises(ValueError):
        test_module.sos2tf_analog(invalid_SOS)

    # Definir una matriz SOS no válida (cada fila no tiene exactamente 6 elementos)
    invalid_SOS = np.array([[1, 0.5, 1, 1, 0.2]])

    # Verificar que se levante un ValueError al pasar una matriz SOS no válida
    with pytest.raises(ValueError):
        test_module.sos2tf_analog(invalid_SOS)

    # Definir una matriz SOS válida pero no como instancia de ndarray
    invalid_SOS = [[1, 0.5, 1, 1, 0.2, 1], [1, 1, 1, 1, 1, 1]]

    # Verificar que se levante un ValueError al pasar una matriz SOS no válida
    with pytest.raises(ValueError):
        test_module.sos2tf_analog(invalid_SOS)


@pytest.mark.parametrize(
    "obj",
    [
        [np.array([1. , 37./21., 229./21., 95./7., 87./7., 9. ]), # num
         np.array([1. , 53./10., 248./5., 146., 520., 800. ])],   # den
        [np.array([1. , 9., 16., 9., 1. ]),    # num
         np.array([1. , 12., 29., 12., 1. ])], # den
        [np.array([1. , 10./7., 6./7., 9. ]), # num
         np.array([1. , 53./10., 248./5., 146., 520., 800. ])],   # den
        [np.array([1. , 4./3., 4./3., 1. ]), # num
         np.array([1. , 53./10., 248./5., 146., 520., 800. ])],   # den
        [np.array([1. , 1., 0. ]), # num
         np.array([1. , 53./10., 248./5., 146., 520., 800. ])],   # den
        [np.array([1. , 2., 0., 0. ]), # num
         np.array([1. , 53./10., 248./5., 146., 520., 800. ])],   # den
        [np.array([1. , 2., 0. ]), # num
         np.array([1. , 8., 22., 24., 9.])],   # den
        [np.array([1. , 1./3., 1.]), # num
         np.array([1. , 8., 22., 24., 9.])],   # den
    ]
)
def test_tf2sos_analog_valid_input(obj):
    
    num, den = obj
    # Coeficientes numéricos y denóminos de una función de transferencia de prueba

    # tolerancia numérica
    tol = 1e-10

    # Llamar a la función tf2sos_analog con los coeficientes de prueba
    result_sos = test_module.tf2sos_analog(num, den)

    # Verificar que el resultado sea una instancia de ndarray
    assert isinstance(result_sos, np.ndarray)

    # Calcular la matriz SOS esperada de otra manera
    sos_num, sos_den = test_module._one_sos2tf(result_sos[0,:])
    expected_tf = TransferFunction(sos_num, sos_den)
    for ii in range(1, result_sos.shape[0]):
        sos_num, sos_den = test_module._one_sos2tf(result_sos[ii,:])
        tfa = TransferFunction(sos_num, sos_den)
        expected_tf = test_module.tfcascade(expected_tf, tfa)

    # Verificar que la forma y los elementos de la matriz SOS sean los esperados
    assert np.max(np.abs( expected_tf.num - num)) < tol
    assert np.max(np.abs( expected_tf.den - den)) < tol

def test_tf2sos_analog_invalid_input():
    # Coeficientes numéricos y denóminos no válidos (no instancias de arrays de numpy)
    invalid_num = 2
    invalid_den = [4, 5, 6]

    # Verificar que se levante un ValueError al pasar coeficientes no válidos
    with pytest.raises(ValueError):
        test_module.tf2sos_analog(invalid_num, invalid_den)

    invalid_num = [4, 5, 6]
    invalid_den = 3

    # Verificar que se levante un ValueError al pasar coeficientes no válidos
    with pytest.raises(ValueError):
        test_module.tf2sos_analog(invalid_num, invalid_den)

    # Coeficientes numéricos y denóminos no válidos (más ceros que polos)
    invalid_num = np.array([1, 2, 3, 4])
    invalid_den = np.array([4, 5, 6])

    # Verificar que se levante un ValueError al pasar coeficientes no válidos
    with pytest.raises(ValueError):
        test_module.tf2sos_analog(invalid_num, invalid_den)

def test_pretty_print_lti_valid_input():
    # Coeficientes numéricos y denóminos de una función de transferencia de prueba
    num = [1, 2, 3]
    den = [4, 5, 6]

    # Llamar a la función pretty_print_lti con los coeficientes de prueba
    result = test_module.pretty_print_lti(num, den, displaystr=False)

    # Verificar que se devuelva una cadena de texto
    assert isinstance(result, str)

    # Verificar que la cadena generada sea la esperada
    expected_str = r'\frac{s^2 \,\, 0.25 + s \,\, 0.5 + 0.75 }{s^2 + s \,\, 1.25 + 1.5 }'
    assert result == expected_str

def test_pretty_print_lti_invalid_num_type():
    # Coeficientes numéricos no válidos (no list o ndarray)
    invalid_num = '1, 2, 3'

    # Verificar que se levante un ValueError al pasar coeficientes numéricos no válidos
    with pytest.raises(ValueError):
        test_module.pretty_print_lti(invalid_num)

def test_pretty_print_lti_invalid_den_type():
    # Coeficientes denóminos no válidos (no list o ndarray)
    invalid_den = '4, 5, 6'

    # Verificar que se levante un ValueError al pasar coeficientes denóminos no válidos
    with pytest.raises(ValueError):
        test_module.pretty_print_lti([1, 2, 3], den=invalid_den)

def test_pretty_print_lti_invalid_displaystr_type():
    # Tipo de displaystr no válido (no bool)
    invalid_displaystr = 123

    # Verificar que se levante un ValueError al pasar displaystr no válido
    with pytest.raises(ValueError):
        test_module.pretty_print_lti([1, 2, 3], displaystr=invalid_displaystr)

## analogicos

# arbitrary
tf_arb = TransferFunction( np.array( [6, 24, 5, 96]),   # num
  np.array( [1,  2, 3,  4]) ) # den

# lowpass
tf_lp = TransferFunction( np.array( [0., 0., 3.] ),      # num
  np.array( [1., 1./3., 1.] ) ) # den

# notch
tf_notch = TransferFunction( np.array( [2, 0, 8] ),     # num
  np.array( [1, 2, 3, 4]) ) # den

# highpass
tf_hp = TransferFunction( np.array( [7., 0., 0.] ),    # num
  np.array( [1., 5./4., 25.] ) ) # den

# hp-notch
tf_hpnotch = TransferFunction( np.array( [5, 0, 80]),    # num
  np.array( [1, 2, 3, 4]) ) # den

# lp-notch
tf_lpnotch = TransferFunction( np.array( [9, 0, 9]),    # num
  np.array( [1, 2, 3, 4]) ) # den

# bilineal
tf_bili = TransferFunction( np.array( [0, 3, 3] ),    # num
  np.array( [0, 1, 2] ) ) # den

# bilineal
tf_bili2 = TransferFunction( np.array( [0, 5, 35]),    # num
  np.array( [0, 1, 3 ] ) ) # den

# lp 1er orden
tf_lp1 = TransferFunction( np.array( [0, 0, 15] ),    # num
  np.array( [0, 1, 3] ) ) # den

# hp 1er orden
tf_hp1 = TransferFunction( np.array( [0, 5, 0] ),   # num
  np.array( [0, 1, 3] ) ) # den
    

# digitales
fs = 1e3

# arbitrario
numz, denz = sig.bilinear(tf_arb.num, tf_arb.den, fs = fs)
tf_arb_dig = TransferFunction( numz, denz, dt=1/fs ) # den

# lowpass
numz, denz = sig.bilinear(tf_lp.num, tf_lp.den, fs = fs)
tf_lp_dig = TransferFunction( numz, denz, dt=1/fs ) # den

# notch
numz, denz = sig.bilinear(tf_notch.num, tf_notch.den, fs = fs)
tf_notch_dig = TransferFunction( numz, denz, dt=1/fs ) # den

# highpass
numz, denz = sig.bilinear(tf_hp.num, tf_hp.den, fs = fs)
tf_hp_dig = TransferFunction( numz, denz, dt=1/fs ) # den

# hp-notch
numz, denz = sig.bilinear(tf_hpnotch.num, tf_hpnotch.den, fs = fs)
tf_hpnotch_dig = TransferFunction( numz, denz, dt=1/fs ) # den

# lp-notch
numz, denz = sig.bilinear(tf_lpnotch.num, tf_lpnotch.den, fs = fs)
tf_lpnotch_dig = TransferFunction( numz, denz, dt=1/fs ) # den

# bilineal
numz, denz = sig.bilinear(tf_bili.num, tf_bili.den, fs = fs)
tf_bili_dig = TransferFunction( numz, denz, dt=1/fs ) # den

# bilineal
numz, denz = sig.bilinear(tf_bili2.num, tf_bili2.den, fs = fs)
tf_bili2_dig = TransferFunction( numz, denz, dt=1/fs ) # den

# lp 1er orden
numz, denz = sig.bilinear(tf_lp1.num, tf_lp1.den, fs = fs)
tf_lp1_dig = TransferFunction( numz, denz, dt=1/fs ) # den

# hp 1er orden
numz, denz = sig.bilinear(tf_hp1.num, tf_hp1.den, fs = fs)
tf_hp1_dig = TransferFunction( numz, denz, dt=1/fs ) # den

@pytest.mark.parametrize(
    "true_parameters, num, den",
    [
       # arbitrary
       ([sp.Rational('2'),  # w_od
          sp.Rational('3'),  # Q_d
          sp.Rational('4'),  # w_on
          sp.Rational('5'),  # Q_n
          sp.Rational('6')], # K
          np.array( [sp.Rational(6), sp.Rational(24,5), sp.Rational(96)] ),    # num
          np.array( [sp.Rational(1), sp.Rational(2,3), sp.Rational(4)] ) ), # den
       
       # lowpass
       ([sp.Rational('1'),  # w_od
          sp.Rational('3'),  # Q_d
          sp.Rational('1'),  # w_on
          sp.oo,  # Q_n
          sp.Rational('3')], # K
          np.array( [0., 0., 3.] ),    # num
          np.array( [1., 1./3., 1.] ) ), # den
       
       # notch
       ([sp.Rational('2'),  # w_od
          sp.Rational('3'),  # Q_d
          sp.Rational('2'),  # w_on
          sp.oo,  # Q_n
          sp.Rational('2')], # K
          np.array( [sp.Rational(2), sp.Rational(0), sp.Rational(8)] ),    # num
          np.array( [sp.Rational(1), sp.Rational(2,3), sp.Rational(4)] ) ), # den

       # highpass
       ([sp.Rational('5'),  # w_od
          sp.Rational('4'),  # Q_d
          sp.Rational('0'),  # w_on
          sp.Rational('0'),  # Q_n
          sp.Rational('7')], # K
          np.array( [7., 0., 0.] ),    # num
          np.array( [1., 5./4., 25.] ) ), # den
       
       # hp-notch
       ([sp.Rational('2'),  # w_od
          sp.Rational('3'),  # Q_d
          sp.Rational('4'),  # w_on
          sp.oo,  # Q_n
          sp.Rational('5')], # K
          np.array( [sp.Rational(5), sp.Rational(0), sp.Rational(80)] ),    # num
          np.array( [sp.Rational(1), sp.Rational(2,3), sp.Rational(4)] ) ), # den
       
       # lp-notch
       ([sp.Rational('2'),  # w_od
          sp.Rational('3'),  # Q_d
          sp.Rational('1'),  # w_on
          sp.oo,  # Q_n
          sp.Rational('9')], # K
          np.array( [sp.Rational(9), sp.Rational(0), sp.Rational(9)] ),    # num
          np.array( [sp.Rational(1), sp.Rational(2,3), sp.Rational(4)] ) ), # den
       
       # bilineal
       ([sp.Rational('2'),  # w_od
          sp.nan,  # Q_d
          sp.Rational('1'),  # w_on
          sp.nan,  # Q_n
          sp.Rational('3')], # K
          np.array( [sp.Rational(0), sp.Rational(3), sp.Rational(3)] ),    # num
          np.array( [sp.Rational(0), sp.Rational(1), sp.Rational(2)] ) ), # den
       
       # bilineal
       ([sp.Rational('3'),  # w_od
          sp.nan,  # Q_d
          sp.Rational('7'),  # w_on
          sp.nan,  # Q_n
          sp.Rational('5')], # K
          np.array( [sp.Rational(0), sp.Rational(5), sp.Rational(35)] ),    # num
          np.array( [sp.Rational(0), sp.Rational(1), sp.Rational(3)] ) ), # den
       
       # lp 1er orden
       ([sp.Rational('3'),  # w_od
          sp.nan,  # Q_d
          sp.Rational('3'),  # w_on
          sp.nan,  # Q_n
          sp.Rational('5')], # K
          np.array( [sp.Rational(0), sp.Rational(0), sp.Rational(15)] ),    # num
          np.array( [sp.Rational(0), sp.Rational(1), sp.Rational(3)] ) ), # den
       
       # hp 1er orden
       ([sp.Rational('3'),  # w_od
          sp.nan,  # Q_d
          sp.Rational('0'),  # w_on
          sp.nan,  # Q_n
          sp.Rational('5')], # K
          np.array( [sp.Rational(0), sp.Rational(5), sp.Rational(0)] ),    # num
          np.array( [sp.Rational(0), sp.Rational(1), sp.Rational(3)] ) ), # den
       
    ]
)
def test_parametrize_sos_valid_input_bicuad( true_parameters, num, den ):
    # Coeficientes de un sistema de segundo orden
    s = test_module.s

    # num = [K, K*w_on/Q_n,  K*w_on**2]
    # den = [1, w_od/Q_d, w_od**2]
    
    w_od, Q_d, w_on, Q_n, K = true_parameters
    
    num_poly = sp.Poly(num[0]*s**2 + num[1] * s + num[2], s)
    den_poly = sp.Poly(den[0]*s**2 + den[1] * s + den[2], s)

    # Llamar a la función parametrize_sos con los coeficientes de prueba
    result_num, result_den, result_w_on, result_Q_n, result_w_od, result_Q_d, result_K  = test_module.parametrize_sos(num_poly, den_poly)

    # Verificar que la tupla devuelta tenga los elementos esperados
    assert (result_K*result_num).expr == num_poly.expr
    assert result_den.expr == den_poly.expr
    assert result_w_od == w_od
    assert result_Q_d == Q_d
    assert result_w_on == w_on
    assert result_Q_n == Q_n
    assert result_K == K
    
def test_parametrize_sos_invalid_num_type():

    s = test_module.s
    
    # Coeficientes numéricos no válidos (no Poly)
    invalid_num = s

    # Verificar que se levante un ValueError al pasar coeficientes numéricos no válidos
    with pytest.raises(ValueError):
        test_module.parametrize_sos(invalid_num, sp.Poly(s + 1))

def test_parametrize_sos_invalid_den_type():
    s = test_module.s
    
    # Coeficientes del denominador no válidos (no Poly)
    invalid_den = s

    # Verificar que se levante un ValueError al pasar coeficientes del denominador no válidos
    with pytest.raises(ValueError):
        test_module.parametrize_sos(sp.Poly(s), invalid_den)

def test_pretty_print_bicuad_omegayq_complete_second_order():
    num = [1, 2/3, 4]
    den = [1, 5/4, 25]
    expected_output = r'\frac{s^2 + s \frac{  2}{  3} +   2^2}{s^2 + s \frac{  5}{  4} +   5^2}'
    assert test_module.pretty_print_bicuad_omegayq(num, den, displaystr=False) == expected_output

def test_pretty_print_bicuad_omegayq_biquad_passband():
    num = [ 1./3.*5./4., 0 ]
    den = [1, 5/4, 25]
    expected_output = r'\frac{s\,0.3333\,\frac{  5}{  4}}{s^2 + s \frac{  5}{  4} +   5^2}'
    assert test_module.pretty_print_bicuad_omegayq(num, den, displaystr=False) == expected_output

def test_pretty_print_bicuad_omegayq_notch():
    num = [3.*1., 0, 3.*25.]
    den = [1, 5/4, 25]
    expected_output = r'\frac{  3(s^2 +   5^2)}{s^2 + s \frac{  5}{  4} +   5^2}'
    assert test_module.pretty_print_bicuad_omegayq(num, den, displaystr=False) == expected_output

def test_pretty_print_bicuad_omegayq_lowpass():
    num = [3.*25.]
    den = [1, 5/4, 25]
    expected_output = r'\frac{ 75 }{s^2 + s \frac{  5}{  4} +   5^2}'
    assert test_module.pretty_print_bicuad_omegayq(num, den, displaystr=False) == expected_output

def test_pretty_print_bicuad_omegayq_highpass():
    num = [3., 0, 0]
    den = [1, 5/4, 25]
    expected_output = r'\frac{s^2 \,\,   3 }{s^2 + s \frac{  5}{  4} +   5^2}'
    assert test_module.pretty_print_bicuad_omegayq(num, den, displaystr=False) == expected_output

def test_pretty_print_SOS_default_mode():
    mySOS = np.array([[1., 7., 1., 1., 9., 1.],
                      [1., 2., 1., 1., 3., 1.]])
    expected_output = r' \frac{s^2 + s \,\,   7 +   1 }{s^2 + s \,\,   9 +   1 } . \frac{s^2 + s \,\,   2 +   1 }{s^2 + s \,\,   3 +   1 }'
    
    assert test_module.pretty_print_SOS(mySOS, displaystr=False) == expected_output

def test_pretty_print_SOS_omegayq_mode():
    mySOS = np.array([[1., 7., 1., 1., 9., 1.],
                      [1., 2., 1., 1., 3., 1.]])
    expected_output = r' \frac{s^2 + s \frac{  1}{0.1429} +   1^2}{s^2 + s \frac{  1}{0.1111} +   1^2} . \frac{s^2 + s \frac{  1}{0.5} +   1^2}{s^2 + s \frac{  1}{0.3333} +   1^2}'
    
    assert test_module.pretty_print_SOS(mySOS, mode='omegayq', displaystr=False) == expected_output

def test_pretty_print_SOS_invalid_mode():
    with pytest.raises(ValueError):
        test_module.pretty_print_SOS(np.array([[1, 2, 1, 1, 2, 1], [1, 1, 1, 1, 1, 1]]), mode='invalid_mode')

def test_pretty_print_SOS_invalid_matrix():
    with pytest.raises(ValueError):
        test_module.pretty_print_SOS([[1, 2, 1, 1, 2, 1], [1, 1, 1, 1, 1, 1]])

def test_pretty_print_SOS_invalid_displaystr():
    with pytest.raises(ValueError):
        test_module.pretty_print_SOS(np.array([[1, 2, 1, 1, 2, 1], [1, 1, 1, 1, 1, 1]]), displaystr='not_bool')

def test_pretty_print_SOS_invalid_shape():
    with pytest.raises(ValueError):
        test_module.pretty_print_SOS(np.array([[1, 2, 1, 1, 2]]))


def test_analyze_sys_single_transfer_function():
    # Crear una única función de transferencia para analizar
    num = [ 1./3.*5./4., 0 ]
    den = [1, 5/4, 25]
    H = TransferFunction(num, den)
    # Llamar a la función analyze_sys con la función de transferencia única
    result = test_module.analyze_sys([H], sys_name='Single Transfer Function')
    # Verificar el tipo y la longitud del resultado
    assert isinstance(result, list)
    assert len(result) == 4  # debe haber cuatro pares de handles de figuras y ejes
    
    freq_resp, pzmap, retardo_digital, retardo_analog = result

    assert len(freq_resp[1]) == 2  # debe haber dos handles de ejes
    assert isinstance(freq_resp[1][0], plt.Axes)
    assert isinstance(freq_resp[1][1], plt.Axes)

    assert isinstance(pzmap[1], plt.Axes)
    assert isinstance(retardo_analog[1], plt.Axes)
    
    
@pytest.mark.parametrize(
    "all_sys, all_lbls",
    [
       ([tf_arb, tf_bili], ['l1', 'l1']),
       ([tf_arb, tf_arb_dig], ['l1', 'l1_dig']),
       ([tf_lp1, tf_hp1], ['l1', 'l2']),
       ([tf_lp1_dig, tf_hp1_dig], ['l1', 'l2']),
       ([tf_lp1, tf_hp1, tf_lp1_dig, tf_hp1_dig], ['l1', 'l2', 'l1_d', 'l2_d']),
       ([tf_notch], ['l1']),
       ([tf_notch_dig], ['l1_d']),
       ([tf_arb, tf_lpnotch, tf_hpnotch], ['l1', 'l2', 'l3'])
       
   ])  
def test_analyze_sys_multiple_transfer_functions(all_sys, all_lbls):
    # Crear varias funciones de transferencia para analizar

    # Llamar a la función analyze_sys con las funciones de transferencia múltiples
    result = test_module.analyze_sys(all_sys, sys_name=all_lbls)
    # Verificar el tipo y la longitud del resultado
    assert isinstance(result, list)
    assert len(result) == 4  # debe haber cuatro pares de handles de figuras y ejes
    
    freq_resp, pzmap_s, pzmap_z, retardo_analog = result

    assert len(freq_resp[1]) == 2  # debe haber dos handles de ejes
    assert isinstance(freq_resp[1][0], plt.Axes)
    assert isinstance(freq_resp[1][1], plt.Axes)

    if not isinstance(pzmap_s[1], plt.Axes):
        assert isinstance(pzmap_z[1], plt.Axes)
        
    assert isinstance(retardo_analog[1], plt.Axes)

@pytest.mark.parametrize(
    "all_sys, all_lbls",
    [
       (test_module.tfcascade(tf_arb, tf_bili), ['l1']),
       (test_module.tfcascade(tf_lp1, tf_hpnotch), ['l2']),
       (test_module.tfcascade(test_module.tfcascade(tf_lp, tf_lpnotch), tf_lp1), ['l1']),
       (test_module.tfcascade(test_module.tfcascade(tf_arb, tf_lpnotch), tf_hpnotch), ['l3'])
   ])  
def test_analyze_sys_sos_matrix(all_sys, all_lbls):
    # Crear una matriz SOS para analizar

    # Llamar a la función analyze_sys con la matriz SOS
    result = test_module.analyze_sys( test_module.tf2sos_analog(all_sys), sys_name=all_lbls)
    # Verificar el tipo y la longitud del resultado
    assert isinstance(result, list)
    assert len(result) == 4  # debe haber cuatro pares de handles de figuras y ejes
    
    freq_resp, pzmap, retardo_digital, retardo_analog = result

    assert len(freq_resp[1]) == 2  # debe haber dos handles de ejes
    assert isinstance(freq_resp[1][0], plt.Axes)
    assert isinstance(freq_resp[1][1], plt.Axes)

    assert isinstance(pzmap[1], plt.Axes)
    assert isinstance(retardo_analog[1], plt.Axes)

def test_analyze_sys_invalid_xaxis():
    # Llamar a la función analyze_sys con un valor de xaxis no válido
    with pytest.raises(ValueError):
        test_module.analyze_sys([], xaxis='invalid_xaxis_value')

# Prueba para asegurar que la función arroja ValueError cuando se le pasan arreglos de diferentes longitudes
def test_group_delay_different_lengths():
    freq = np.linspace(0, 10, 100)
    phase = np.sin(freq)
    with pytest.raises(ValueError):
        test_module.group_delay(freq[:-1], phase)

# Prueba para asegurar que la función arroja ValueError cuando se le pasan argumentos que no son arreglos NumPy
def test_group_delay_non_numpy_arrays():
    freq = [0, 1, 2, 3]
    phase = [0, 1, 2, 3]
    with pytest.raises(ValueError):
        test_module.group_delay(freq, phase)

# Prueba para asegurar que la función devuelve un arreglo con la misma longitud que los arreglos de entrada
def test_group_delay_output_length():
    freq = np.linspace(0, 10, 100)
    phase = np.sin(freq)
    group_delay_result = test_module.group_delay(freq, phase)
    assert len(group_delay_result) == len(freq)

# Prueba para asegurar que la función calcula correctamente el retardo de grupo
def test_group_delay_calculation():
    # tolerancia numérica
    tol = 1e-10
    freq = np.arange(0, -10, step = -1/10)
    group_delay_result = test_module.group_delay(-freq, freq)

    assert np.max(np.abs( group_delay_result - 1)) < tol


# Prueba para verificar si se genera un gráfico sin errores cuando se pasan argumentos válidos
def test_plot_plantilla_valid_arguments():
    try:
        test_module.plot_plantilla()
        test_module.plot_plantilla(filter_type='lowpass', fpass=0.25, ripple=0.5, fstop=0.6, attenuation=40, fs=2)
        test_module.plot_plantilla(filter_type='highpass', fpass=0.25, ripple=0.5, fstop=0.6, attenuation=40, fs=2)
        test_module.plot_plantilla(filter_type='bandpass', fpass=(0.2, 0.4), ripple=0.3, fstop=(0.1, 0.5), attenuation=50, fs=2)
        test_module.plot_plantilla(filter_type='bandstop', fpass=(0.2, 0.4), ripple=0.3, fstop=(0.1, 0.5), attenuation=50, fs=2)
    except Exception as e:
        pytest.fail(f"Se generó una excepción: {e}")

# Prueba para verificar si se genera un ValueError cuando se pasan argumentos inválidos
def test_plot_plantilla_invalid_arguments():
    with pytest.raises(ValueError):
        test_module.plot_plantilla(filter_type='lowpass', fpass='invalid', ripple=0.5, fstop=0.6, attenuation=40, fs=2)
    with pytest.raises(ValueError):
        test_module.plot_plantilla(filter_type='bandpass', fpass=(0.4,), ripple=0.3, fstop=(0.1, 0.5), attenuation=50, fs=2)
    with pytest.raises(ValueError):
        test_module.plot_plantilla(filter_type='bandstop', fpass=(0.4,), ripple=0.3, fstop=(0.1, 0.5), attenuation=50, fs=2)

# Prueba para verificar si se genera un gráfico sin errores cuando se pasa un tipo de filtro desconocido
def test_plot_plantilla_unknown_filter_type():
    try:
        test_module.plot_plantilla(fpass=0.25, ripple=0.5, fstop=0.6, attenuation=40, fs=2)
    except Exception as e:
        pytest.fail(f"Se generó una excepción: {e}")

# Prueba para verificar si se genera un ValueError cuando se pasan argumentos no numéricos
def test_plot_plantilla_non_numeric_arguments():
    with pytest.raises(ValueError):
        test_module.plot_plantilla(filter_type='lowpass', fpass='invalid', ripple=0.5, fstop=0.6, attenuation=40, fs=2)
    with pytest.raises(ValueError):
        test_module.plot_plantilla(filter_type='lowpass', fpass=0.25, ripple='invalid', fstop=0.6, attenuation=40, fs=2)
    with pytest.raises(ValueError):
        test_module.plot_plantilla(filter_type='lowpass', fpass=0.25, ripple=0.5, fstop='invalid', attenuation=40, fs=2)
    with pytest.raises(ValueError):
        test_module.plot_plantilla(filter_type='lowpass', fpass=0.25, ripple=0.5, fstop=0.6, attenuation='invalid', fs=2)
    with pytest.raises(ValueError):
        test_module.plot_plantilla(filter_type='lowpass', fpass=0.25, ripple=0.5, fstop=0.6, attenuation=40, fs='invalid')


# cerramos todas las figuras que se hayan abierto, para no perjudicar otras pruebas
def test_last_call():

    plt.close('all')
