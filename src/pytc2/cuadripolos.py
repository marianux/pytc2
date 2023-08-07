#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:14:13 2023

@author: mariano
"""

import numpy as np

import sympy as sp

from .general import print_latex


'''
    Funciones de conversión de matrices de cuadripolos lineales
'''

def S2Ts_s(Spar):
    '''
    Convierte una matriz de parámetros scattering (S) simbólica 
    al modelo de parámetros transferencia de scattering (Ts).

    Parameters
    ----------
    Spar : Symbolic Matrix
        Matriz de parámetros S.

    Returns
    -------
    Ts : Symbolic Matrix
        Matriz de parámetros de transferencia scattering.

    '''
    
    Ts = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Ts[0,0] = sp.Rational('1')
    # B = DZ/Z21
    Ts[0,1] = -Spar[1,1]
    # C = 1/Z21
    Ts[1,0] = Spar[0,0]
    # D = Z22/Z21
    Ts[1,1] = -sp.simplify(sp.expand(sp.Determinant(Spar)))
    
    return( sp.simplify(sp.expand(1 / Spar[1,0] * Ts) ) ) 

def Ts2S_s(Ts):
    '''
    Convierte una matriz de transferencia de scattering (Ts) simbólica 
    al modelo de parámetros scattering (S).

    Parameters
    ----------
    Ts : Symbolic Matrix
        Matriz de parámetros S.

    Returns
    -------
    Spar : Symbolic Matrix
        Matriz de parámetros de scattering.

    '''
    
    Spar = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Spar[0,0] = Ts[1,0]
    # B = DZ/Z21
    Spar[0,1] = sp.simplify(sp.expand(sp.Determinant(Ts)))
    # C = 1/Z21
    Spar[1,0] = sp.Rational('1') 
    # D = Z22/Z21
    Spar[1,1] = -Ts[0,1] 
    
    return( sp.simplify(sp.expand( 1 / Ts[0,0] * Spar ) ) ) 

def Ts2Tabcd_s(Ts, Z0 = sp.Rational('1') ):
    '''
    Convierte una matriz de transferencia de scattering (Ts) simbólica 
    al modelo ABCD o Tabcd simbólico (Tabcd).

    Parameters
    ----------
    Ts : Symbolic Matrix
        Matriz de parámetros Ts.

    Returns
    -------
    Tabcd : Symbolic Matrix
        Matriz de parámetros de ABCD.

    '''
    
    return(sp.simplify(sp.expand(S2Tabcd_s(Ts2S_s(Ts), Z0))))


def Tabcd2S_s(Tabcd, Z0 = sp.Rational('1') ):
    '''
    Convierte una matriz de parámetros ABCD (Tabcd) simbólica 
    al modelo de parámetros scattering (S).

    Parameters
    ----------
    Ts : Symbolic Matrix
        Matriz de parámetros S.

    Returns
    -------
    Spar : Symbolic Matrix
        Matriz de parámetros de scattering.

    '''
    
    Spar = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Spar[0,0] = Tabcd[0,0] + Tabcd[0,1]/Z0 - Tabcd[1,0]*Z0 - Tabcd[1,1]
    # B = DZ/Z21
    Spar[0,1] = 2 * sp.simplify(sp.expand(sp.Determinant(Tabcd)))
    # C = 1/Z21
    Spar[1,0] = sp.Rational('2') 
    # D = Z22/Z21
    Spar[1,1] = -Tabcd[0,0] + Tabcd[0,1]/Z0 - Tabcd[1,0]*Z0 + Tabcd[1,1]
    
    common = Tabcd[0,0] + Tabcd[0,1]/Z0 + Tabcd[1,0]*Z0 + Tabcd[1,1]
    
    return( sp.simplify(sp.expand( 1 / common * Spar ) ) ) 

def S2Tabcd_s(Spar, Z0 = sp.Rational('1') ):
    '''
    Convierte una matriz de parámetros ABCD (Tabcd) simbólica 
    al modelo de parámetros scattering (S).

    Parameters
    ----------
    Ts : Symbolic Matrix
        Matriz de parámetros S.

    Returns
    -------
    Spar : Symbolic Matrix
        Matriz de parámetros de scattering.

    '''
    
    Tabcd = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Tabcd[0,0] = (1 - Spar[0,0]) * (1 + Spar[1,1]) + Spar[1,0] * Spar[0,1]
    # B = DZ/Z21
    Tabcd[0,1] = Z0*((1 + Spar[0,0]) * (1 + Spar[1,1]) - Spar[1,0] * Spar[0,1])
    # C = 1/Z21
    Tabcd[1,0] = 1/Z0*((1 - Spar[0,0]) * (1 - Spar[1,1]) - Spar[1,0] * Spar[0,1])
    # D = Z22/Z21
    Tabcd[1,1] = (1 - Spar[0,0]) * (1 + Spar[1,1]) + Spar[1,0] * Spar[0,1]
    
    return( sp.simplify(sp.expand( 1 / 2 / Spar[1,0] * Tabcd ) ) ) 

def SparZ_s(Zexc, Z01=sp.Rational('1'), Z02=sp.Rational('1')):
    '''
    Convierte una matriz de transferencia de scattering (Ts) simbólica 
    al modelo de parámetros scattering (S).

    Parameters
    ----------
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.

    Z01 : Symbolic impedance
          Impedancia de referencia en el plano 1

    Z02 : Symbolic impedance
          Impedancia de referencia en el plano 2

    Returns
    -------
    Spar : Symbolic Matrix
           Matriz de parámetros de scattering de Z.

    '''
    
    Spar = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Spar[0,0] = Zexc
    # B = DZ/Z21
    Spar[0,1] = 2*Z01
    # C = 1/Z21
    Spar[1,0] = 2*Z01
    # D = Z22/Z21
    Spar[1,1] = Zexc 
    
    return( sp.simplify(sp.expand( 1 / (Zexc + 2*Z01) * Spar) ) ) 

def SparY_s(Yexc, Y01=sp.Rational('1'), Y02=sp.Rational('1')):
    '''
    Convierte una matriz de transferencia de scattering (Ts) simbólica 
    al modelo de parámetros scattering (S).

    Parameters
    ----------
    Yexc : Symbolic impedance
           Función de excitación de la admitancia a representar.

    Y01 : Symbolic impedance
          Admitancia de referencia en el plano 1

    Y02 : Symbolic impedance
          Admitancia de referencia en el plano 2

    Returns
    -------
    Spar : Symbolic Matrix
           Matriz de parámetros de scattering de Y.

    '''
    
    Spar = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Spar[0,0] = -Yexc
    # B = DZ/Z21
    Spar[0,1] = 2*Y01
    # C = 1/Z21
    Spar[1,0] = 2*Y01
    # D = Z22/Z21
    Spar[1,1] = -Yexc 
    
    return( sp.simplify(sp.expand( 1 / (Yexc + 2*Y01) * Spar) ) ) 

def TabcdLYZ_s(Yexc, Zexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido por una Y en derivación seguida  por 
    una Z en serie.

    Parameters
    ----------
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.
    
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.

    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    '''
    
    Tpar = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Tpar[0,0] = sp.Rational('1')  
    # B = DZ/Z21
    Tpar[0,1] = Zexc
    # C = 1/Z21
    Tpar[1,0] = Yexc
    # D = Z22/Z21
    Tpar[1,1] = sp.Rational('1') + sp.simplify(sp.expand(Zexc * Yexc))
    
    return( Tpar ) 

def TabcdLZY_s(Zexc, Yexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido por una Z en serie seguida una Y en 
    derivación.

    Parameters
    ----------
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.
    
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.

    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    '''
    
    Tpar = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Tpar[0,0] = sp.Rational('1') + sp.simplify(sp.expand(Zexc * Yexc))
    # B = DZ/Z21
    Tpar[0,1] = Zexc
    # C = 1/Z21
    Tpar[1,0] = Yexc
    # D = Z22/Z21
    Tpar[1,1] = sp.Rational('1')  
    
    return( Tpar ) 

def TabcdZ_s(Zexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido únicamente por una Z en serie.

    Parameters
    ----------
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.


    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    '''
    
    Tpar = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Tpar[0,0] = sp.Rational('1') 
    # B = DZ/Z21
    Tpar[0,1] = Zexc
    # C = 1/Z21
    Tpar[1,0] = sp.Rational('0') 
    # D = Z22/Z21
    Tpar[1,1] = sp.Rational('1')  
    
    return( Tpar ) 

def TabcdY_s(Yexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido únicamente por una Y en derivación.

    Parameters
    ----------
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.


    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    '''
    
    Tpar = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Tpar[0,0] = sp.Rational('1') 
    # B = DZ/Z21
    Tpar[0,1] = sp.Rational('0')
    # C = 1/Z21
    Tpar[1,0] = Yexc
    # D = Z22/Z21
    Tpar[1,1] = sp.Rational('1')  
    
    return( Tpar ) 

def Y2Tabcd_s(YY):
    """
    
    Parameters
    ----------
    tfa : TYPE
        DESCRIPTION.
    tfb : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    Example
    -------

    """
    
    TT = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Y22/Y21
    TT[0,0] = sp.simplify(sp.expand(-YY[1,1]/YY[1,0]))
    # B = -1/Y21
    TT[0,1] = sp.simplify(sp.expand(-1/YY[1,0]))
    # C = -DY/Y21
    TT[1,0] = sp.simplify(sp.expand(-sp.Determinant(YY)/YY[1,0]))
    # D = Y11/Y21
    TT[1,1] = sp.simplify(sp.expand(-YY[1,1]/YY[1,0]))
    
    return(TT)

def Z2Tabcd_s(ZZ):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    
    TT = sp.Matrix([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    TT[0,0] = sp.simplify(sp.expand(ZZ[0,0]/ZZ[1,0]))
    # B = DZ/Z21
    TT[0,1] = sp.simplify(sp.expand(sp.Determinant(ZZ)/ZZ[1,0]))
    # C = 1/Z21
    TT[1,0] = sp.simplify(sp.expand(1/ZZ[1,0]))
    # D = Z22/Z21
    TT[1,1] = sp.simplify(sp.expand(ZZ[1,1]/ZZ[1,0]))
    
    return(TT)

def Tabcd2Z_s(TT):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    
    ZZ = sp.Matrix([[0, 0], [0, 0]])
    
    # Z11 = A/C
    ZZ[0,0] = sp.simplify(sp.expand(TT[0,0]/TT[1,0]))
    # Z11 = DT/C
    ZZ[0,1] = sp.simplify(sp.expand(sp.Determinant(TT)/TT[1,0]))
    # Z21 = 1/C
    ZZ[1,0] = sp.simplify(sp.expand(1/TT[1,0]))
    # Z22 = D/C
    ZZ[1,1] = sp.simplify(sp.expand(TT[1,1]/TT[1,0]))
    
    return(ZZ)

def Tabcd2Y_s(TT):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    
    YY = sp.Matrix([[0, 0], [0, 0]])
    
    # Y11 = D/B
    YY[0,0] = sp.simplify(sp.expand(TT[1,1]/TT[0,1]))
    # Y12 = -DT/B
    YY[0,1] = sp.simplify(sp.expand(-sp.Determinant(TT)/TT[0,1]))
    # Y21 = -1/B
    YY[1,0] = sp.simplify(sp.expand(-1/TT[0,1]))
    # Y22 = A/B
    YY[1,1] = sp.simplify(sp.expand(TT[0,0]/TT[0,1]))
    
    return(YY)

def I2Tabcd_s(gamma, z01, z02 = None):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    gamma : Symbol
        Ganancia compleja expresada en neppers (Re{gamma}) y rad (Im{gamma}).
        
    z01 : Symbol
        Nivel de impedancia del puerto 1.

    z02 : Symbol
        (Opcional) Nivel de impedancia del puerto 2. Default: z02 = z01 = Zo (imp. característica). 

    Returns
    -------
    TT : Symbolic Matrix
        Matriz ABCD en función de los parámetros imagen 

    '''
    
    if z02 is None:
        z02 = z01
    
    TT = sp.Matrix([[sp.cosh(gamma)*sp.sqrt(z01/z02),
                     sp.sinh(gamma)*sp.sqrt(z01*z02)], 
                    [sp.sinh(gamma)/sp.sqrt(z01*z02),
                     sp.cosh(gamma)*sp.sqrt(z02/z01)]])
    
    
    return(TT)


def Model_conversion( src_model, dst_model  ):
    '''
    Convierte modelos de cuadripolos lineales.

    Parameters
    ----------
    src_model : Dict. 
        Diccionario que describe al modelo de origen.
        
    dst_model : Dict. 
        Diccionario que describe al modelo de salida.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    Example
    -------

    # Parámetros Z (impedancia - circ. abierto)
    ZZ = sp.Matrix([[z11, z12], [z21, z22]])
    # vars. dependientes
    vv = sp.Matrix([[v1], [v2]])
    # vars. INdependientes
    ii = sp.Matrix([[i1], [i2]])
    
    
    # Parámetros Tdcba (Transmisión inversos, DCBA)
    TTi = sp.Matrix([[Ai, Bi], [-Ci, -Di]])
    # vars. dependientes
    ti_dep = sp.Matrix([[v2], [i2]])
    # vars. INdependientes. (Signo negativo de corriente)
    ti_ind = sp.Matrix([[v1], [i1]])
    
    # Diccionario con la definición de cada modelo
    src_model = { 'model_name': 'Z', 'matrix': ZZ, 'dep_var': vv, 'indep_var':ii }
    dst_model = { 'model_name': 'T', 'matrix': TT, 'dep_var': t_dep, 'indep_var':t_ind, 'neg_i2_current': True }
    
    T_z = convert_params( src_model, dst_model )
    
    print_latex(a_equal_b_latex_s('T_Z', T_z ))


    '''

    
    if src_model['model_name'] == dst_model['model_name']:
        
        # mismo modelo -> salimos
        QQ = src_model['matrix']
        
        return({'matrix': QQ, 'name': dst_model['model_name'] + '_{' + src_model['model_name'] + '}' })

    if 'proxy_matrix' in src_model:
        # Uso una representación Tabcd en el caso de parámetros S o Ts
        src_matrix = src_model['proxy_matrix']
        
    else:
        src_matrix = src_model['matrix']
    
    aa = sp.solve([ src_matrix * src_model['indep_var'] - src_model['dep_var']
                ], 
                dst_model['dep_var'])

    
    # reemplazaremos el determinante por Delta.
    det_src_matrix = sp.det(src_matrix)
    
    # i2 se define al revés en este modelo
    if 'neg_i2_current' in src_model:
        det_src_matrix = -det_src_matrix
    
    
    dd = sp.Symbol('\Delta')
    jj = 0
    
    QQ = sp.Matrix([[0,0],[0,0]])




    for dep_var in dst_model['dep_var']:

        jj += 1
        
        yyy = sp.collect(sp.expand(aa[dep_var]), dst_model['indep_var'][0])
        yyy = sp.collect(yyy, dst_model['indep_var'][1])

        # i2 se define al revés en este modelo
        if dep_var.name == 'i2' and 'neg_i2_current' in dst_model:
            yyy = -yyy

        kk = 0
        for indep_var in dst_model['indep_var']:
            kk += 1
            
            bb = sp.cancel(yyy.coeff(indep_var,1))
            # i2 se define al revés en este modelo
            if indep_var.name == 'i2' and 'neg_i2_current' in dst_model:
                bb = -bb
                
            QQ[jj-1,kk-1] = bb.subs(det_src_matrix, dd)
            #print_latex(a_equal_b_latex_s( param_lab + '{:d}{:d}'.format(jj,kk), QQ[jj-1,kk-1] ))

    return({'matrix': QQ, 'name': dst_model['model_name'] + '_{' + src_model['model_name'] + '}' })


def Y2Tabcd(YY):
    """
    
    Parameters
    ----------
    tfa : TYPE
        DESCRIPTION.
    tfb : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    Example
    -------

    """
    
    TT = np.zeros_like(YY)
    
    # A = Y22/Y21
    TT[0,0] = -YY[1,1]/YY[1,0]
    # B = -1/Y21
    TT[0,1] = -1/YY[1,0]
    # C = -DY/Y21
    TT[1,0] = -np.linalg.det(YY)/YY[1,0]
    # D = Y11/Y21
    TT[1,1] = -YY[1,1]/YY[1,0]
    
    return(TT)

def Z2Tabcd(ZZ):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    
    TT = np.zeros_like(ZZ)
    
    # A = Z11/Z21
    TT[0,0] = ZZ[0,0]/ZZ[1,0]
    # B = DZ/Z21
    TT[0,1] = np.linalg.det(ZZ)/ZZ[1,0]
    # C = 1/Z21
    TT[1,0] = 1/ZZ[1,0]
    # D = Z22/Z21
    TT[1,1] = ZZ[1,1]/ZZ[1,0]
    
    return(TT)

def Tabcd2Z(TT):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    
    ZZ = np.zeros_like(TT)
    
    # Z11 = A/C
    ZZ[0,0] = TT[0,0]/TT[1,0]
    # Z11 = DT/C
    ZZ[0,1] = np.linalg.det(TT)/TT[1,0]
    # Z21 = 1/C
    ZZ[1,0] = 1/TT[1,0]
    # Z22 = D/C
    ZZ[1,1] = TT[1,1]/TT[1,0]
    
    return(ZZ)

def Tabcd2Y(TT):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    YY = np.zeros_like(TT)
    
    # Y11 = D/B
    YY[0,0] = TT[1,1]/TT[0,1]
    # Y12 = -DT/B
    YY[0,1] = -np.linalg.det(TT)/TT[0,1]
    # Y21 = -1/B
    YY[1,0] = -1/TT[0,1]
    # Y22 = A/B
    YY[1,1] = TT[0,0]/TT[0,1]
    
    return(YY)

def y2mai(YY):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    
    Ymai = YY.row_insert(YY.shape[0], sp.Matrix([-sum(YY[:,ii] ) for ii in range(YY.shape[1])]).transpose() )
    Ymai = Ymai.col_insert(Ymai.shape[1], sp.Matrix([-sum(Ymai[ii,:] ) for ii in range(Ymai.shape[0])]) )
    Ymai[-1] = sum(YY)
    
    return(Ymai)

def may2y(Ymai, nodes2del):
    '''
    Convierte la MAI en MAD luego de remover filas y columnas indicadas en nodes2del

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    
    YY = Ymai
    
    for ii in nodes2del:
        YY.row_del(ii)
    
    for ii in nodes2del:
        YY.col_del(ii)
    
    return(YY)

def TabcdLYZ(Yexc, Zexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido por una Y en derivación seguida  por 
    una Z en serie.

    Parameters
    ----------
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.
    
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.

    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    '''
    
    Tpar = np.array([[0, 0], [0, 0]])
    
    # A = Z11/Z21
    Tpar[0,0] = 1 
    # B = DZ/Z21
    Tpar[0,1] = Zexc
    # C = 1/Z21
    Tpar[1,0] = Yexc
    # D = Z22/Z21
    Tpar[1,1] = 1 + Zexc * Yexc
    
    return( Tpar ) 

def TabcdLZY(Zexc, Yexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido por una Z en serie seguida una Y en 
    derivación.

    Parameters
    ----------
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.
    
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.

    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    '''
    
    Tpar = np.array([[0.0, 0], [0, 0]])
    
    # A = Z11/Z21
    Tpar[0,0] = 1 + Zexc * Yexc 
    # B = DZ/Z21
    Tpar[0,1] = Zexc
    # C = 1/Z21
    Tpar[1,0] = Yexc
    # D = Z22/Z21
    Tpar[1,1] = 1
    
    return( Tpar ) 

def TabcdZ(Zexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido únicamente por una Z en serie.

    Parameters
    ----------
    Zexc : Symbolic impedance
           Función de excitación de la impedancia a representar.


    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    '''
    
    Tpar = np.array([[0.0, 0.0], [0, 0]])
    
    # A = Z11/Z21
    Tpar[0,0] = 1 
    # B = DZ/Z21
    Tpar[0,1] = Zexc
    # C = 1/Z21
    Tpar[1,0] = 0
    # D = Z22/Z21
    Tpar[1,1] = 1

    return( Tpar ) 

def TabcdY(Yexc):
    '''
    Implementa una matriz de transferencia ABCD (Tabcd) a partir de 
    un cuadripolo constituido únicamente por una Y en derivación.

    Parameters
    ----------
    Yexc : Symbolic admitance
           Función de excitación de la admitancia a representar.


    Returns
    -------
    Tabcd : Symbolic Matrix
           Matriz de parámetros ABCD.

    '''
    
    Tpar = np.array([[0.0, 0], [0, 0]])
    
    # A = Z11/Z21
    Tpar[0,0] = 1 
    # B = DZ/Z21
    Tpar[0,1] = 0
    # C = 1/Z21
    Tpar[1,0] = Yexc
    # D = Z22/Z21
    Tpar[1,1] = 1
    
    return( Tpar ) 

def I2Tabcd(gamma, z01, z02 = None):
    '''
    Convierte la MAD en MAI luego de levantar de referencia.

    Parameters
    ----------
    Ymai : Symbolic Matrix
        Matriz admitancia indefinida.
    nodes2del : list or integer
        Nodos que se van a eliminar.

    Returns
    -------
    YY : Symbolic Matrix
        Matriz admitancia 

    '''
    if z02 is None:
        z02 = z01

    # if np.sqrt(z02/z01)
    
    TT = np.matrix([[np.cosh(gamma)*np.sqrt(z01/z02),
                     np.sinh(gamma)*np.sqrt(z01*z02)], 
                    [np.sinh(gamma)/np.sqrt(z01*z02),
                     np.cosh(gamma)*np.sqrt(z02/z01)]])
    
    return(TT)


def calc_MAI_ztransf_ij_mn(Ymai, ii=2, jj=3, mm=0, nn=1, verbose=False):
    """Calcula la transferencia de impedancia V_ij / I_mn

    Parameters
    ----------
    tfa : TYPE
        DESCRIPTION.
    tfb : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    Example
    -------
    
    """

    """
    """
    
    if ii > jj:
        max_ouput_idx = ii
        min_ouput_idx = jj
    else:
        max_ouput_idx = jj
        min_ouput_idx = ii
    
    if mm > nn:
        max_input_idx = mm
        min_input_idx = nn
    else:
        max_input_idx = nn
        min_input_idx = mm
    
    # cofactor de 2do orden
    num = Ymai.minor_submatrix(max_ouput_idx, max_input_idx).minor_submatrix(min_ouput_idx, min_input_idx)
    # cualquier cofactor de primer orden
    den = Ymai.minor_submatrix(min_input_idx, min_input_idx)

    num_det = sp.simplify(num.det())
    den_det = sp.simplify(den.det())
    
    sign_correction = mm+nn+ii+jj
    Tz = sp.simplify(-1**(sign_correction) * num_det/den_det)
    
    if( verbose ):
    
        print_latex(r' [Y_{MAI}] = ' + sp.latex(Ymai) )
        
        print_latex(r' [Y^{{ {:d}{:d} }}_{{ {:d}{:d} }} ] = '.format(mm,nn,ii,jj) + sp.latex(num) )
    
        print_latex(r'[Y^{{ {:d} }}_{{ {:d} }}] = '.format(mm,mm) + sp.latex(den) )
    
        print_latex(r'\mathrm{{Tz}}^{{ {:d}{:d} }}_{{ {:d}{:d} }} = \frac{{ \underline{{Y}}^{{ {:d}{:d} }}_{{ {:d}{:d} }} }}{{ \underline{{Y}}^{{ {:d} }}_{{ {:d} }} }} = '.format(ii,jj,mm,nn,mm,nn,ii,jj,mm,mm) + r' -1^{{ {:d} }} '.format(sign_correction)  + r'\frac{{ ' + sp.latex(num_det) + r'}}{{' + sp.latex(den_det) + r'}} = ' + sp.latex(Tz))
    
    return(Tz)

def calc_MAI_vtransf_ij_mn(Ymai, ii=2, jj=3, mm=0, nn=1, verbose=False):
    """Calcula la transferencia de tensión V_ij / V_mn
    
    Parameters
    ----------
    tfa : TYPE
        DESCRIPTION.
    tfb : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    Example
    -------
    
    """
    
    if ii > jj:
        max_ouput_idx = ii
        min_ouput_idx = jj
    else:
        max_ouput_idx = jj
        min_ouput_idx = ii
    
    if mm > nn:
        max_input_idx = mm
        min_input_idx = nn
    else:
        max_input_idx = nn
        min_input_idx = mm
    
    # cofactores de 2do orden
    num = Ymai.minor_submatrix(max_ouput_idx, max_input_idx).minor_submatrix(min_ouput_idx, min_input_idx)

    den = Ymai.minor_submatrix(max_input_idx, max_input_idx).minor_submatrix(min_input_idx, min_input_idx)
    
    num_det = sp.simplify(num.det())
    den_det = sp.simplify(den.det())
    
    sign_correction = mm+nn+ii+jj
    Av = sp.simplify(-1**(sign_correction) * num_det/den_det)
    
    if( verbose ):
    
        print_latex(r' [Y_{MAI}] = ' + sp.latex(Ymai) )
        
        print_latex(r' [Y^{{ {:d}{:d} }}_{{ {:d}{:d} }} ] = '.format(mm,nn,ii,jj) + sp.latex(num) )
    
        print_latex(r'[Y^{{ {:d}{:d} }}_{{ {:d}{:d} }} ] = '.format(mm,nn,mm,nn) + sp.latex(den) )
    
        print_latex(r'T^{{ {:d}{:d} }}_{{ {:d}{:d} }} = \frac{{ \underline{{Y}}^{{ {:d}{:d} }}_{{ {:d}{:d} }} }}{{ \underline{{Y}}^{{ {:d}{:d} }}_{{ {:d}{:d} }} }} = '.format(ii,jj,mm,nn,mm,nn,ii,jj,mm,nn,mm,nn) + r' -1^{{ {:d} }} '.format(sign_correction)  + r'\frac{{ ' + sp.latex(num_det) + r'}}{{' + sp.latex(den_det) + r'}} = ' + sp.latex(Av) )
    
    return(Av)

def calc_MAI_impedance_ij(Ymai, ii=0, jj=1, verbose=False):
    """Calcula la transferencia de tensión V_ij / V_mn
    
    Parameters
    ----------
    tfa : TYPE
        DESCRIPTION.
    tfb : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    Example
    -------
    
    """

    if ii > jj:
        max_idx = ii
        min_idx = jj
    else:
        max_idx = jj
        min_idx = ii
 
    # cofactor de 2do orden
    num = Ymai.minor_submatrix(max_idx, max_idx).minor_submatrix(min_idx, min_idx)
    # cualquier cofactor de primer orden
    den = Ymai.minor_submatrix(min_idx, min_idx)
    
    ZZ = sp.simplify(num.det()/den.det())
    
    if( verbose ):

        print_latex(r' [Y_{MAI}] = ' + sp.latex(Ymai) )
        
        print_latex(r' [Y^{{ {:d}{:d} }}_{{ {:d}{:d} }} ] = '.format(ii,ii,jj,jj) + sp.latex(num) )

        print_latex(r'[Y^{{ {:d} }}_{{ {:d} }}] = '.format(ii,ii) + sp.latex(den) )

        print_latex(r'Z_{{ {:d}{:d} }} = \frac{{ \underline{{Y}}^{{ {:d}{:d} }}_{{ {:d}{:d} }} }}{{ \underline{{Y}}^{{ {:d} }}_{{ {:d} }} }} = '.format(ii,jj,ii,ii,jj,jj,ii,ii) + sp.latex(ZZ))

    return(ZZ)

