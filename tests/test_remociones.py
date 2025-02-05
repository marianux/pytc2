#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:58:49 2024

@author: mariano
"""

import pytest
import sympy as sp
import numpy as np
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



@pytest.mark.parametrize(
    "S21_fac, bNumeric",
    [
        (   ( s**2 + sp.Rational(4)) * ( s**2 + sp.Rational(1/15)*s + sp.Rational(1)) / ( s**2 + sp.Rational(1/2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(5)*s + sp.Rational(1)) / (s+sp.Rational(1)),
            True),
        (   ( s**2 + sp.Rational(4)) * ( s**2 + sp.Rational(1/15)*s + sp.Rational(1)) / ( s**2 + sp.Rational(1/2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(5)*s + sp.Rational(1)) / (s+sp.Rational(1)), 
            False),
        (   ( s**4 ) / ( s**2 + sp.sqrt(2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(1/2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(5)*s + sp.Rational(1)) / (s+sp.Rational(1)), 
            False),
        (   ( s**3 * ( s**2 + sp.Rational(9))) / ( s**2 + sp.sqrt(2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(1/2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(5)*s + sp.Rational(1)) / (s+sp.Rational(1)), 
            False),
        (   ( s**2 * ( s**2 + sp.Rational(9))) / ( s**2 + sp.sqrt(2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(1/2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(5)*s + sp.Rational(1)) / (s+sp.Rational(1)), 
            False),
        (   (s**2 * ( s**2 + 9)) / ( s**2 + np.sqrt(2)*s + 1) / ( s**2 + 1/2*s + 1) / ( s**2 + 5*s + 1) / (s+1), 
            True),
        (   (s * ( s**2 + sp.Rational(9))) / ( s**2 + sp.sqrt(2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(1/2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(5)*s + sp.Rational(1)) / (s+sp.Rational(1)), 
            False),
        (   (s * ( s**2 + 9)) / ( s**2 + np.sqrt(2)*s + 1) / ( s**2 + 1/2*s + 1) / ( s**2 + 5*s + 1) / (s+1), 
            True),
        (   ( s**2 + sp.Rational(1/4)) * ( s**2 + sp.Rational(1/15)*s + sp.Rational(1)) / ( s**2 + sp.Rational(1/2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(5)*s + sp.Rational(1)) / s, 
            False),
        (   (s**2 * ( s**2 + sp.Rational(1/15)*s + sp.Rational(1))) / ( s**2 + sp.Rational(1/2)*s + sp.Rational(1)) / ( s**2 + sp.Rational(5)*s + sp.Rational(1)) / (s+sp.Rational(1)), 
            False),
    ]
)   
def test_modsq2mod_s_valid(S21_fac, bNumeric):
    

    # print(sp.simplify(sp.expand(S21_fac)))
        
    S21sq = sp.simplify(sp.expand(S21_fac * S21_fac.subs(s, -s)))

    print(sp.simplify(sp.expand(S21sq)))

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        factor_func = test_module.modsq2mod_s( S21sq, bTryNumeric = bNumeric )
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")

    # Verificar que el tipo de resultado sea sp.Matrix
    assert isinstance(factor_func, sp.Expr)

    # Verificar el correcto resultado
    if bNumeric:
        s21_num, s21_den = sp.fraction(S21_fac)

        s21_num_coeffs = np.array(s21_num.as_poly(s).all_coeffs(), dtype=float)
        s21_den_coeffs = np.array(s21_den.as_poly(s).all_coeffs(), dtype=float)
        
        f_num, f_den = sp.fraction(factor_func)

        f_num_coeffs = np.array(f_num.as_poly(s).all_coeffs(), dtype=float)
        f_den_coeffs = np.array(f_den.as_poly(s).all_coeffs(), dtype=float)
        
        assert np.allclose(s21_num_coeffs, f_num_coeffs)
    
        assert np.allclose(s21_den_coeffs, f_den_coeffs)
    
    else:
                
        assert sp.simplify(S21_fac - factor_func) == sp.Rational(0)
    
    
def test_modsq2mod_s_invalid_input():

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        factor_func = test_module.modsq2mod_s( 1.0 )


def test_isFRP_valid():

    immit = (s**2 + 4*s + 3)/(s**2 + 2*s)    

    # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
    try:
        test_module.isFRP(immit)
    except ValueError:
        pytest.fail("Se levantó un ValueError incorrectamente.")

    # immit si es FRP
    assert test_module.isFRP(immit)

    immit = (s**2 - 4*s + 3)/(s**2 - 2*s)
    
    # immit no es FRP
    assert not test_module.isFRP(immit)

    
def test_isFRP_invalid_input():
    
    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        test_module.isFRP(1.)



@pytest.mark.parametrize(
    "funciones_parametros, true_val",
    [
        ({
            test_module.remover_polo_sigma: {'immit' : (s**2 + 13*s + 32)/(2*(s+1)*(s+6)), 
                                             'sigma' : -1, 
                                             'isImpedance' : True,
                                             'isRC' : True,
                                             'sigma_zero' : None,
                                             },
        }, ((s + 8)/(2*(s + 6)), 2/(s + 1)) ),
        ({
            test_module.remover_polo_sigma: {'immit' : (s**2 + 13*s + 32)/(2*(s+1)*(s+6)), 
                                             'sigma' : -6, 
                                             'isImpedance' : True,
                                             'isRC' : True,
                                             'sigma_zero' : -5.,
                                             },
        }, (2.5*(0.2*s + 1.0)/(s + 1), 1.0/(s + 6.0), 1/6, 1.) ),
        ({
            test_module.remover_polo_sigma: {'immit' : (s**2 + 13*s + 32)/(2*(s+1)*(s+6)), 
                                             'sigma' : -1, 
                                             'isImpedance' : False,
                                             'isRC' : False,
                                             'sigma_zero' : None,
                                             },
        }, (2/3*(1/8*s + 1.0)/(1/6*s + 1.), 2.0/(s + 1.0), 1/2, 1/2) ),
        ({
            test_module.remover_polo_sigma: {'immit' : (s**2 + 13*s + 32)/(2*(s+1)*(s+6)), 
                                             'sigma' : -6., 
                                             'isImpedance' : False,
                                             'isRC' : False,
                                             'sigma_zero' : -5.,
                                             },
        }, (2.5*(0.2*s + 1.0)/(s + 1), 1.0/(s + 6.0), 6., 1.) ),
        ({
            test_module.remover_polo_sigma: {'immit' : 2*(s + 1)*(s + 3)/(s + 2)/(s + 6), 
                                             'sigma' : -6., 
                                             'isImpedance' : False,
                                             'isRC' : True,
                                             'sigma_zero' : None,
                                             },
        }, (0.5*(0.75*s + 1.0)/(0.5*s + 1.0), 1.25*s/(s + 6.0), 0.8, 0.208333333333333) ),
        ({
            test_module.remover_polo_sigma: {'immit' : 2*(s + 1)*(s + 3)/(s + 2)/(s + 6), 
                                             'sigma' : -6., 
                                             'isImpedance' : True,
                                             'isRC' : False,
                                             'sigma_zero' : None,
                                             },
        }, (0.5*(0.75*s + 1.0)/(0.5*s + 1.0), 1.25*s/(s + 6.0), 1.25, 0.208333333333333) ),
        ({
            test_module.remover_polo_sigma: {'immit' : 2*(s + 1)*(s + 3)/(s + 2)/(s + 6), 
                                             'sigma' : -6., 
                                             'isImpedance' : False,
                                             'isRC' : True,
                                             'sigma_zero' : None,
                                             },
        }, (0.5*(0.75*s + 1.0)/(0.5*s + 1.0), 1.25*s/(s + 6.0), 0.8, 0.208333333333333) ),
        ({
            test_module.remover_polo_sigma: {'immit' : (s + 8)/(s+4), 
                                             'sigma' : -4., 
                                             'isImpedance' : True,
                                             'isRC' : True,
                                             'sigma_zero' : None,
                                             },
         }, (1., 4./(s + 4.), 1., 0.25) ),
        ({
            test_module.remover_valor_en_dc: {'immit' : (3*s**2 + 27*s+ 44)/(s**2 + 13*s + 32), 
                                             'sigma_zero' : None,
                                             },
          }, (s*(13*s + 73)/(8*(s**2 + 13*s + 32)), 11/8) ),
        ({
            test_module.remover_valor_en_infinito: {'immit' : (s**2 + 13*s + 32)/(3*s**2 + 27*s+ 44), 
                                             'sigma_zero' : None
                                             },
         }, (4*(3*s + 13)/(3*(3*s**2 + 27*s + 44)), 1/3) ),
        ({
            test_module.remover_polo_infinito: {'immit' : ((s**2+2)*(s**2+5))/(3*s*(s**2+sp.Rational(7,3))), 
                                             'omega_zero' : None,
                                             'isSigma' : False
                                             },
         }, (2*(7*s**2 + 15)/(3*s*(3*s**2 + 7)), s/3) ),
        ({
            test_module.remover_polo_dc: {'immit' : ((s**2+2)*(s**2+5))/(3*s*(s**2+sp.Rational(7,3))), 
                                             'omega_zero' : 1., 
                                             'isSigma' : False,
                                             },
         }, ((s**2 + 1)*(s**2 + 3)/(s*(3*s**2 + 7)), 1/s) ),
        ({
            test_module.remover_polo_jw: {'immit' : (s * (3*s**2+7) )/((s**2+1)*(s**2+3)), 
                                             'omega' : 1., 
                                             'isImpedance' : True,
                                             'omega_zero' : None,
                                             },
         }, (s/(s**2 + 3), 2*s/(s**2 + 1)) ),
        ({
            test_module.remover_polo_jw: {'immit' : (s * (3*s**2+7) )/((s**2+1)*(s**2+3)), 
                                             'omega' : 1., 
                                             'isImpedance' : False,
                                             'omega_zero' : None,
                                             },
         }, (1/3*s/(1/3*s**2 + 1.0), 2.0*s/(s**2 + 1.0), 0.5, 2.) ),

   ]
)   
def test_remociones_valid(funciones_parametros, true_val):

    for funcion, parametros in funciones_parametros.items():

        # Verificar que no se levante un ValueError al pasar funciones de transferencia válidas
        try:
            out = funcion(**parametros)
        except ValueError:
            pytest.fail("Se levantó un ValueError incorrectamente.")
        
        # Verificar que el tipo de resultado sea sp.Matrix
        assert isinstance(out[0], sp.Expr)
        assert isinstance(out[1], sp.Expr)
    
        # Verificar el correcto resultado
        assert sp.simplify(out[0] - true_val[0]) == sp.Rational(0)
        assert sp.simplify(out[1] - true_val[1]) == sp.Rational(0)


# @pytest.mark.parametrize(
#     "func_ptr, args, true_val",
#     [
#        (   
#         ), 
     
        
       
#     ]
# )   
# def test_remociones_valid(func_ptr, args, true_val):
    
    
#     immit, sigma_omega = args

#     # verificar que no se levante un valueerror al pasar funciones de transferencia válidas
#     try:
#         out = func_ptr(immit, sigma_omega)
#     except valueerror:
#         pytest.fail("se levantó un valueerror incorrectamente.")
        
#     # verificar que el tipo de resultado sea sp.matrix
#     assert isinstance(out[0], sp.expr)
#     assert isinstance(out[1], sp.expr)

#     # verificar el correcto resultado
#     assert sp.simplify(out[0] - true_val[0]) == sp.rational(0)
#     assert sp.simplify(out[1] - true_val[1]) == sp.rational(0)
    
    
def test_remover_polo_sigma_invalid_input():
    
    func_ptr = test_module.remover_polo_sigma
    
    immit = (s**2 + 13*s + 32)/(2*(s+1)*(s+6))

    sigma_R1C1 = -1

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a', sigma_R1C1)

    with pytest.raises(ValueError):
        out = func_ptr(immit, 'a')

    with pytest.raises(ValueError):
        out = func_ptr(immit, sigma_R1C1, isImpedance = 2)

    with pytest.raises(ValueError):
        out = func_ptr(immit, sigma_R1C1, isRC = 2)

    with pytest.raises(ValueError):
        out = func_ptr(immit, sigma_R1C1, sigma_zero = 'a')

def test_remover_polo_jw_invalid_input():
    
    func_ptr = test_module.remover_polo_jw
    
    immit = (s * (3*s**2+7) )/((s**2+1)*(s**2+3))

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a')

    with pytest.raises(ValueError):
        out = func_ptr(immit, omega = 'a')

    with pytest.raises(ValueError):
        out = func_ptr(immit, isImpedance = 2)

    with pytest.raises(ValueError):
        out = func_ptr(immit, omega_zero = 'a')


def test_remover_polo_jw_invalid_input():
    
    func_ptr = test_module.remover_polo_jw
    
    immit = (s * (3*s**2+7) )/((s**2+1)*(s**2+3))

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a')

    with pytest.raises(ValueError):
        out = func_ptr(immit, omega = 'a')

    with pytest.raises(ValueError):
        out = func_ptr(immit, isImpedance = 2)


    with pytest.raises(ValueError):
        out = func_ptr(immit, omega_zero = 'a')

@pytest.mark.parametrize(
    "func_ptr, immit",
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
def test_remover_polo_dc_invalid_input(func_ptr, immit):

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a')

    with pytest.raises(ValueError):
        out = func_ptr(immit, omega_zero = 'a')

    with pytest.raises(ValueError):
        out = func_ptr(immit, isSigma = 2)

@pytest.mark.parametrize(
    "func_ptr, immit",
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
def test_remover_valor_invalid_input(func_ptr, immit):

    # Verificar que se levante un ValueError 
    with pytest.raises(ValueError):
        out = func_ptr('a')

    with pytest.raises(ValueError):
        out = func_ptr(immit, sigma_zero = 'a')


# def test__valid():
    
# def test__invalid_input():
