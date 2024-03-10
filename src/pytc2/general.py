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

s = sp.symbols('s', complex=True) 
"""
Variable compleja de Laplace s = σ + j.ω
En caso de necesitar usarla, importar el símbolo desde este módulo.
"""

w = sp.symbols('w', complex=False) 
"""
Fourier real variable ω 
En caso de necesitar usarla, importar el símbolo desde este módulo.
"""

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
    >>> import sympy as sp
    >>> from pytc2.general import pp
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
    1/(z1*z2*(1/z2 + 1/z1))
    
    """
    if not ( (isinstance(z1, sp.Expr) and isinstance(z2, sp.Expr)) or
              (isinstance(z1, (float, complex)) and isinstance(z2, (float, complex))) 
            ):
        raise ValueError('z1 y z2 deben ser AMBOS de tipo Symbolic o float')
        
    return(z1*z2/(z1+z2))
    

#%%
  ##################################
 ## Funciones para uso simbólico ##
##################################

#%%


def simplify_n_monic(tt):
    '''
    Simplifica un polinomio de fracciones en forma mónica.


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
    >>> import sympy as sp
    >>> from pytc2.general import simplify_n_monic
    >>> s = sp.symbols('s')
    >>> tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
    >>> simplified_tt = simplify_n_monic(tt)
    >>> print(simplified_tt)
    (s + 2)/(2*s + 3)
    
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
    >>> from pytc2.general import Chebyshev_polynomials
    >>> Ts = Chebyshev_polynomials(3)
    >>> print(Ts)
    w*(4*w**2 - 3)

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
    
    .. math:: a = b
    
    en un nuevo string formateado para visualizarse en LaTeX.


    Parameters
    ----------
    a : Symbolic or str
        Símbolo o cadena para el lado izquierdo de la igualdad.
    b : Symbolic, str o lista de ambas
        Símbolo o cadena para el lado derecho de la igualdad.


    Returns
    -------
    str: string
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
    >>> import sympy as sp
    >>> from pytc2.general import a_equal_b_latex_s, print_latex
    >>> s = sp.symbols('s')
    >>> tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
    >>> print(a_equal_b_latex_s(sp.symbols('tt'), tt))
    'tt=\\frac{s^{2} + 3 s + 2}{2 s^{2} + 5 s + 3}$'
    >>> print_latex(a_equal_b_latex_s(sp.symbols('tt'), tt))
    [LaTex formated equation]

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
    >>> import sympy as sp
    >>> from pytc2.general import expr_simb_expr, print_latex
    >>> s = sp.symbols('s')
    >>> tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
    >>> tt1 = (s**2 + 4*s + 7) / (2*s**2 + 5*s + 3)
    >>> print_latex(expr_simb_expr('tt', tt1, r'\\neq'))  
    [LaTex formated equation]
    >>> print_latex(expr_simb_expr('tt', tt))
    [LaTex formated equation]


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
    >>> import sympy as sp
    >>> from pytc2.general import to_latex, print_latex
    >>> print(to_latex(sp.symbols('x')))
    $x$
    >>> print_latex(to_latex(sp.symbols('x')))
    [LaTex formated equation]

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
    
    >>> import sympy as sp
    >>> from pytc2.general import str_to_latex, print_latex
    >>> print(str_to_latex('x'))
    $x$
    >>> print_latex(str_to_latex('x'))  
    [LaTex formated equation]

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
    >>> import sympy as sp
    >>> from pytc2.general import str_to_latex, print_latex
    >>> print(str_to_latex('x'))
    $x$
    >>> print_latex(str_to_latex('x'))  
    [LaTex formated equation]

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
    >>> from pytc2.general import print_console_alert
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
    >>> from pytc2.general import print_console_subtitle
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
    >>> from pytc2.general import print_subtitle
    >>> print_subtitle('Subtítulo')
    <IPython.core.display.Markdown object>

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
      TypeError: Si at_en_db no es de tipo `float`.


    See Also
    --------
    :func:`nepper2db`


    Examples
    --------
    >>> from pytc2.general import db2nepper
    >>> db2nepper(20.)
    2.3025850929940455
    >>> db2nepper(1.)
    0.11512925464970228

    '''

    if not isinstance(at_en_db, float):
        raise TypeError('at_en_db debe ser float')

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
      TypeError: Si at_en_db no es de tipo `float`.
    

    See Also
    --------
    :func:`db2nepper`


    Examples
    --------
    >>> from pytc2.general import nepper2db
    >>> nepper2db(1.)
    8.685889638065037
    >>> nepper2db(2.3025850929940455)
    20.

    '''

    if not isinstance(at_en_np, float):
        raise TypeError('at_en_np debe ser float')

    return at_en_np * (20 * np.log10(np.exp(1)))
    