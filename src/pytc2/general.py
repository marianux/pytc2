#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:22:31 2023

Originally based on the work of Combination of 2011 Christopher Felton
Further modifications were added for didactic purposes
by Mariano Llamedo llamedom _at_ frba_utn_edu_ar

@author: marianux
"""

import sympy as sp
import numpy as np

from IPython.display import display, Math, Markdown

  ##########################################
 ## Variables para el análisis simbólico ##
##########################################

# Laplace complex variable. s = σ + j.ω
s = sp.symbols('s', complex=True) # Laplace complex variable. s = σ + j.ω
# Fourier real variable ω 
w = sp.symbols('w', complex=False) # Fourier real variable ω 

#%%
  #########################
 ## Funciones generales ##
#########################
#%%
  
def pp(z1, z2):
    """
    Asocia en paralelo dos impedancias o en serie dos admitancias.
    
    Parameters
    ----------
    z1 : Symbolic o float
        Inmitancia 1.
    z2 : Symbolic o float
        Inmitancia 2.
    

    Returns
    -------
    zp : Symbolic o float
        Inmitancia resultante.


    Raises
    ------
      TypeError: Si alguno de los argumentos no es de tipo `Symbolic`.
    
    
    See Also
    -----------
    :func:`print_latex`
    :func:`to_latex`
    :func:`a_equal_b_latex_s`

    
    Examples
    --------
    >>> # Asociación en paralelo de dos impedancias
    >>> z1 = sp.symbols('z1')
    >>> z2 = sp.symbols('z2')
    >>> zp = pp(z1, z2)
    >>> print(zp)
    z1*z2/(z1 + z2)
    >>> # Asociación en serie de dos admitancias
    >>> y1 = 1/z1
    >>> y2 = 1/z2
    >>> yp = pp(y1, y2)
    >>> print(yp)
    (y1 + y2)**(-1)
    
    """
    if not ((isinstance(z1, sp.Symbol) and isinstance(z2, sp.Symbol)) or
            (isinstance(z1, (np.floating, np.complexfloating)) and isinstance(z2, (np.floating, np.complexfloating)))):
        raise ValueError('z1 y z2 deben ser AMBOS de tipo Symbolic o np.floating')
        
    return(z1*z2/(z1+z2))
    

#%%
  ##################################
 ## Funciones para uso simbólico ##
##################################

#%%


def simplify_n_monic(tt):
    '''
    Simplifica un polinomio de fracciones en forma monica.

    Parameters
    ----------
    tt : Expr
        Polinomio de fracciones a simplificar.

    Returns
    -------
    Expr
        Polinomio simplificado en forma monica.

    Raises
    ------
    TypeError
        Si la entrada no es una expresión simbólica.

    See Also
    --------
    :func:`print_latex`
    :func:`to_latex`
    :func:`a_equal_b_latex_s`

    Examples
    --------
    >>> from sympy import symbols
    >>> s = symbols('s')
    >>> tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
    >>> simplified_tt = simplify_n_monic(tt)
    >>> print(simplified_tt)
    (s + 1)/(2*s + 1)

    '''
    if not isinstance(tt, sp.Expr):
        raise TypeError("La entrada debe ser una expresión simbólica.")

    # Obtener el numerador y el denominador de la expresión y convertirlos en polinomios
    num, den = sp.fraction(sp.simplify(sp.expand(tt)))
    num = sp.poly(num, s)
    den = sp.poly(den, s)
    
    # Calcular el coeficiente principal del numerador y el denominador
    k = num.LC() / den.LC()
    
    # Convertir el numerador y el denominador a forma monica
    num = num.monic()
    den = den.monic()

    # Devolver el polinomio simplificado en forma monica
    return sp.Mul(k, num/den, evaluate=False)

def Chebyshev_polynomials(nn):
    '''
    Calcula el polinomio de Chebyshev de grado nn.

    Parameters
    ----------
    nn : int
        Grado del polinomio de Chebyshev.

    Returns
    -------
    Ts : Symbolic Matrix
        Matriz de parámetros de transferencia scattering.

    Raises
    ------
    ValueError
        Si nn no es un entero positivo.

    See Also
    --------
    :func:`print_latex`
    :func:`to_latex`
    :func:`a_equal_b_latex_s`
    

    Examples
    --------
    >>> Ts = Chebyshev_polynomials(3)
    >>> print(Ts)
    4*w**3 - 3*w

    '''
    
    if not isinstance(nn, int) or nn < 0:
        raise ValueError("nn debe ser un entero positivo.")

    if nn == 0:
        return sp.Rational(1)
    elif nn == 1:
        return w
    else:
        Cn_pp = sp.Rational(1)
        Cn_p = w
        
        for ii in range(nn-1):
            Cn = sp.Rational(2) * w * Cn_p - Cn_pp
            Cn_pp = Cn_p
            Cn_p = Cn
            
        return sp.simplify(sp.expand(Cn))
    

def a_equal_b_latex_s(a, b):
    '''
    A partir de un string o expresión de SymPy (a), y otra expresión de SymPy (b):
    
    a = b
    
    en un nuevo string formateado para visualizarse en LaTeX.

    Parameters
    ----------
    a : Symbolic or str
        Símbolo o cadena para el lado izquierdo de la igualdad.
    b : Symbolic, str o lista de ambas
        Símbolo o cadena para el lado derecho de la igualdad.

    Returns
    -------
    str
        String formateado en LaTeX representando la igualdad.

    Raises
    ------
    TypeError
        Si a no es un símbolo ni una cadena.
        Si b no es un símbolo.

    See Also
    --------
    :func:`expr_simb_expr`
    :func:`print_latex`
    :func:`to_latex`

    Examples
    --------
    >>> print(a_equal_b_latex_s(sp.symbols('x'), sp.symbols('y')))
    $x=y$
    >>> print(a_equal_b_latex_s('a', sp.symbols('b')))
    $a=b$

    '''

    if not (isinstance(a, (sp.Expr, str)) and isinstance(b, (list, sp.Expr))):
        raise TypeError("a debe ser un símbolo o una cadena y b debe ser un símbolo.")
    
    a_str = sp.latex(a) if isinstance(a, sp.Basic) else a
    
    return '$' + a_str + '=' + sp.latex(b) + '$'

def expr_simb_expr(a, b, symbol='='):
    '''
    A partir de un string o expresión de SymPy (a), y otra expresión de SymPy (b):
    
    a symbol b
    
    en un nuevo string formateado para visualizarse en LaTeX.

    Parameters
    ----------
    a : Symbolic or str
        Símbolo o cadena para el lado izquierdo de la expresión.
    b : Symbolic or str
        Símbolo o cadena para el lado derecho de la expresión.
    symbol : str, optional
        Símbolo de operación entre a y b (por defecto es '=').

    Returns
    -------
    str
        String formateado en LaTeX representando la expresión.

    Raises
    ------
    TypeError
        Si a no es un símbolo ni una cadena.
        Si b no es un símbolo.


    See Also
    --------
    :func:`a_equal_b_latex_s`
    :func:`print_latex`
    :func:`to_latex`


    Examples
    --------
    >>> print(expr_simb_expr(sp.symbols('x'), sp.symbols('y')))
    $x=y$
    >>> print(expr_simb_expr('a', sp.symbols('b'), '!='))  # Usando otro símbolo
    $a\neq b$

    '''

    if not (isinstance(a, (sp.Expr, str)) and isinstance(b, sp.Expr)):
        raise TypeError("a debe ser un símbolo o una cadena y b debe ser un símbolo.")
    
    a_str = sp.latex(a) if isinstance(a, sp.Basic) else a
    
    return '$' + a_str + symbol + sp.latex(b) + '$'

def to_latex(unsimbolo):
    '''
    Convierte un símbolo en un string formateado para visualizarse en LaTeX.

    Parameters
    ----------
    unsimbolo : Symbolic or str
        Símbolo o cadena a convertir a formato LaTeX.

    Returns
    -------
    str
        String formateado en LaTeX.

    Raises
    ------
    TypeError
        Si unsimbolo no es un símbolo ni una cadena.


    See Also
    --------
    :func:`print_latex`
    :func:`str_to_latex`
    :func:`to_latex`


    Examples
    --------
    >>> print(to_latex(sp.symbols('x')))
    $x$
    >>> print(to_latex('a'))  
    $a$

    '''

    if not isinstance(unsimbolo, (sp.Expr, str)):
        raise TypeError("unsimbolo debe ser un símbolo o una cadena.")
    
    return '$' + sp.latex(unsimbolo) + '$'

def str_to_latex(unstr):
    '''
    Formatea un string para visualizarse en LaTeX.

    Parameters
    ----------
    unstr : str
        Cadena a formatear para visualización en LaTeX.

    Returns
    -------
    str
        String formateado en LaTeX.

    Raises
    ------
    TypeError
        Si unstr no es una cadena.


    See Also
    --------
    :func:`print_latex`
    :func:`print_subtitle`
    :func:`to_latex`


    Examples
    --------
    >>> print(str_to_latex('x'))
    $x$
    >>> print(str_to_latex('a'))  
    $a$

    '''

    if not isinstance(unstr, str):
        raise TypeError("unstr debe ser una cadena.")
    
    return '$' + unstr + '$'

def print_latex(unstr):
    '''
    Muestra una expresión LaTeX en formato matemático.

    Parameters
    ----------
    unstr : str
        Cadena que representa la expresión LaTeX.

    Returns
    -------
    None
        Esta función no devuelve nada, simplemente muestra la expresión en formato LaTeX.

    Raises
    ------
    TypeError
        Si unstr no es una cadena.


    See Also
    --------
    :func:`print_subtitle`
    :func:`str_to_latex`
    :func:`to_latex`


    Examples
    --------
    >>> print_latex('\\frac{1}{2}')
    $\frac{1}{2}$

    '''
    if not isinstance(unstr, str):
        raise TypeError("unstr debe ser una cadena.")
    

    display(Math(unstr))
    
#%%
  ###############################################
 ## funciones para presentación de resultados ##
###############################################

#%%
 
def print_console_alert(unstr):
    '''
    Imprime una cadena rodeada por símbolos de alerta en la consola.

    Parameters
    ----------
    unstr : str
        Cadena a imprimir.

    Returns
    -------
    None
        Esta función no devuelve nada, simplemente imprime la cadena en la consola.

    Raises
    ------
    TypeError
        Si unstr no es una cadena.


    See Also
    --------
    :func:`print_subtitle`
    :func:`print_latex`
    :func:`print_alert`


    Examples
    --------
    >>> print_console_alert('Advertencia: Datos incompletos')
    ##################################
    # Advertencia: Datos incompletos #
    ##################################

    '''

    if not isinstance(unstr, str):
        raise TypeError("unstr debe ser una cadena.")

    unstr = '# ' + unstr + ' #\n'
    unstr1 =  '#' * (len(unstr)-1) + '\n' 
    
    print( '\n\n' + unstr1 + unstr + unstr1 )

    
def print_console_subtitle(unstr):
    '''
    Imprime un subtítulo en la consola.

    Parameters
    ----------
    unstr : str
        Cadena que representa el subtítulo.

    Returns
    -------
    None
        Esta función no devuelve nada, simplemente imprime el subtítulo en la consola.

    Raises
    ------
    TypeError
        Si unstr no es una cadena.


    See Also
    --------
    :func:`print_subtitle`
    :func:`print_latex`
    :func:`print_console_alert`


    Examples
    --------
    >>> print_console_subtitle('Subtítulo')
    Subtítulo
    ---------

    '''

    if not isinstance(unstr, str):
        raise TypeError("unstr debe ser una cadena.")

    unstr = unstr + '\n'
    unstr1 =  '-' * (len(unstr)-1) + '\n' 
    
    print( '\n\n' + unstr + unstr1 )

def print_subtitle(unstr):
    '''
    Imprime un subtítulo.

    Parameters
    ----------
    unstr : str
        Cadena que representa el subtítulo.

    Returns
    -------
    None
        Esta función no devuelve nada, simplemente imprime el subtítulo.

    Raises
    ------
    TypeError
        Si unstr no es una cadena.


    See Also
    --------
    :func:`print_latex`
    :func:`print_console_alert`
    :func:`print_console_subtitle`


    Examples
    --------
    >>> print_subtitle('Subtítulo')
    #### Subtítulo

    '''

    if not isinstance(unstr, str):
        raise TypeError("unstr debe ser una cadena.")
    
    display(Markdown('#### ' + unstr))
    
#%%

  ###########################################
 ## funciones para conversión de unidades ##
###########################################

#%%

def db2nepper(at_en_db):
    '''
    Convierte una magnitud en decibels a su equivalente en nepers.

    Parameters
    ----------
    at_en_db : float or numpy.ndarray
        Magnitud en decibelios a convertir.

    Returns
    -------
    float or numpy.ndarray
        Equivalente en nepers.

    Raises
    ------
      TypeError: Si at_en_db no es de tipo `np.floating`.


    See Also
    --------
    :func:`nepper2db`


    Examples
    --------
    >>> db2nepper(20)
    8.685889638065036
    >>> db2nepper(np.array([10, 30, 50]))
    array([4.34294482, 13.02883445, 21.71472408])

    '''

    if not isinstance(at_en_db, np.floating):
        raise TypeError('at_en_db debe ser np.floating')

    return at_en_db / (20 * np.log10(np.exp(1)))

def nepper2db(at_en_np):
    '''
    Convierte una magnitud en neperios a su equivalente en decibelios.

    Parameters
    ----------
    at_en_np : float or numpy.ndarray
        Magnitud en neperios a convertir.

    Returns
    -------
    float or numpy.ndarray
        Equivalente en decibelios.

    Raises
    ------
      TypeError: Si at_en_db no es de tipo `np.floating`.
    

    See Also
    --------
    :func:`db2nepper`


    Examples
    --------
    >>> nepper2db(10)
    22.046285227092834
    >>> nepper2db(np.array([4, 8, 12]))
    array([8.81851409, 17.63702818, 26.45554226])

    '''

    if not isinstance(at_en_np, np.floating):
        raise TypeError('at_en_np debe ser np.floating')

    return at_en_np * (20 * np.log10(np.exp(1)))
    