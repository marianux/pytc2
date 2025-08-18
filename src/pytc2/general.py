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
from sympy import Pow
import numpy as np
from scipy.signal import TransferFunction
from numbers import Integral, Real, Complex

from IPython.display import display, Math, Markdown
import os
import warnings

#%%
  ##############################################
 ## Variables para el funcionamiento general ##
##############################################
#%%

# Obtener el directorio base o raíz del código fuente
pytc2_full_path = os.path.dirname(os.path.abspath(__file__))
"""
Path a donde se encuentra pyTC2 localmente.
"""

small_val = np.finfo(float).eps
"""
Es un valor muy pequeño para que las funciones que tienen restringido su evaluación 
en 0 no arrojen warnings ni errores. e.g. los cálculos de los logaritmos
"""

#%%
  ##########################################
 ## Variables para el análisis simbólico ##
##########################################
#%%

s = sp.symbols('s', complex=True) 
"""
Variable compleja de Laplace s = σ + j.ω
En caso de necesitar usarla, importar el símbolo desde este módulo.
"""

w = sp.symbols('w', real=True) 
"""
Fourier real variable ω 
En caso de necesitar usarla, importar el símbolo desde este módulo.
"""

#%%
  #########################
 ## Funciones generales ##
#########################
#%%

def get_home_directory():
    
    if os.name == 'posix':  # Linux/MacOS
        return os.environ['HOME']
    elif os.name == 'nt':  # Windows
        return os.path.expanduser('~')
    else:
        raise RuntimeError("Unsupported operating system")


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
      ValueError: Si alguno de los argumentos no es de tipo `Symbolic`.
    
    
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
              (isinstance(z1, (Real, Complex)) and isinstance(z2, (Real, Complex))) 
            ):
        raise ValueError('z1 y z2 deben ser AMBOS de tipo Symbolic o float')
        
    return(z1*z2/(z1+z2))



#%%
  ##################################
 ## Funciones para uso simbólico ##
##################################

#%%

def factorSOS(ratfunc, decimals = 4):
    '''
    Factoriza una función racional simbólica, en polinomios de segundo y primer
    orden. 


    Parameters
    ----------
    ratfunc : Expr. simbólica
        Función racional simbólica.
    decimals : entero
        Cantidad de decimales para la evaluación simbólica.


    Returns
    -------
    Expr. simbólica
        Función racional simbólica factorizada.


    Raises
    ------
    ValueError
        Si la entrada no es una expresión simbólica.


    See Also
    --------
    :func:`symbfunc2tf`
    :func:`simplify_n_monic`
    :func:`a_equal_b_latex_s`


    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.general import s, factorSOS
    >>> tt = (s**4 + 8*s**3 + 18*s**2 + 11*s + 2)/(s**3 + 16*s**2 + 65*s + 14)
    >>> factorized_tt, _, _ = factorSOS(tt)
    >>> print(factorized_tt)
    (s + 0.382)*(s + 0.438)*(s + 2.62)*(s + 4.56)/((s + 0.228)*(s + 7.0)*(s + 8.77))
    
    '''
    if not isinstance(ratfunc, sp.Expr):
        raise ValueError("La entrada debe ser una expresión simbólica.")
    
    if not isinstance(decimals, Integral):
        raise ValueError("La cantidad de decimales debe ser un número entero.")

    num, den = sp.fraction(ratfunc)
    
    num = sp.Poly(num,s)
    den = sp.Poly(den,s)

    polySOS = num.LC() / den.LC()
    
    raices = sp.roots(num, s)
    
    # Separa las raíces complejas conjugadas y las raíces reales
    raices_complejas_conjugadas_num = []
    raices_reales_num = []
    
    for raiz, multiplicidad in raices.items():
        if raiz.is_real:
            raices_reales_num.extend([raiz]*multiplicidad)
            polySOS = polySOS * (s - raiz.evalf(decimals))**(multiplicidad)
        else:
           # Busca si ya hay un grupo para la parte real
            grupo_existente = False
            for grupo in raices_complejas_conjugadas_num:
                # pregunto por la parte Real.
                if sp.ask(sp.Q.real((grupo + raiz))):
                    grupo_existente = True
                    break
            if not grupo_existente:
                raices_complejas_conjugadas_num.extend([raiz]*multiplicidad)
                raices_complejas_conjugadas_num.extend([sp.conjugate(raiz)]*multiplicidad)
                this_sos = sp.simplify(sp.expand((s - raiz) * (s - sp.conjugate(raiz))) )
                polySOS = polySOS * this_sos.evalf(decimals)**(multiplicidad)
                

    raices = sp.roots(den, s)
    
    # Separa las raíces complejas conjugadas y las raíces reales
    raices_complejas_conjugadas_den = []
    raices_reales_den = []
    
    for raiz, multiplicidad in raices.items():
        if raiz.is_real:
            raices_reales_den.extend([raiz]*multiplicidad)
            polySOS = polySOS / (s - raiz.evalf(decimals))**(multiplicidad)
        else:
           # Busca si ya hay un grupo para la parte real
            grupo_existente = False
            for grupo in raices_complejas_conjugadas_den:
                # pregunto por la parte Real.
                if sp.ask(sp.Q.real((grupo + raiz))):
                    grupo_existente = True
                    break
            if not grupo_existente:
                raices_complejas_conjugadas_den.extend([raiz]*multiplicidad)
                raices_complejas_conjugadas_den.extend([sp.conjugate(raiz)]*multiplicidad)
                this_sos = sp.simplify(sp.expand((s - raiz) * (s - sp.conjugate(raiz))) )
                polySOS = polySOS / this_sos.evalf(decimals)**(multiplicidad)


    return(polySOS, [ [raices_reales_num],[raices_reales_den] ], [[raices_complejas_conjugadas_num], [raices_complejas_conjugadas_den]])

def symbfunc2tf(tt):
    '''
    Convierte una función racional simbólica, con coeficientes numéricos 
    (convertibles a flotante), en un objeto transfer function.


    Parameters
    ----------
    tt : Expr. simbólica
        Función racional simbólica.


    Returns
    -------
    TransferFunction
        TransferFunction que representa numéricamente la función.


    Raises
    ------
    ValueError
        Si la entrada no es una expresión simbólica.


    See Also
    --------
    :func:`simplify_n_monic`
    :func:`to_latex`
    :func:`factorSOS`


    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.general import s, symbfunc2tf
    >>> tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
    >>> simplified_tt = symbfunc2tf(tt)
    >>> print(simplified_tt)
    TransferFunctionContinuous(
    array([0.5, 1.5, 1. ]),
    array([1. , 2.5, 1.5]),
    dt: None
    )    
    '''
    if not isinstance(tt, sp.Expr):
        raise ValueError("La entrada debe ser una expresión simbólica.")

    num, den = sp.fraction(tt)
        
    aa = np.array( [ float(ii) for ii in num.as_poly(s).all_coeffs()])
    bb = np.array( [ float(ii) for ii in den.as_poly(s).all_coeffs()])
    
    cc = TransferFunction(aa, bb)

    return cc


def flatten_pow(expr):
    """
    Convierte potencias anidadas (a**b)**c  ->  a**(b*c)
    Repite la transformación hasta que ya no cambie la expresión.
    """
    prev = None
    cur = expr
    # Reemplazo repetido hasta estabilidad
    while prev != cur:
        prev = cur
        cur = cur.replace(
            lambda e: isinstance(e, Pow) and isinstance(e.base, Pow),
            lambda e: Pow(e.base.base, sp.simplify(e.base.exp * e.exp))
        )
    return cur

def _symbolic_gt(a, b):
    """Intenta decidir si a > b (para exponentes simbólicos) con fallback numérico y lexicográfico."""
    diff = sp.simplify(a - b)
    if diff.is_positive:
        return True
    if diff.is_negative:
        return False
    # fallback: probar con una sustitución numérica (asumiendo símbolos positivos)
    syms = list(diff.free_symbols)
    if syms:
        subs = {s: 5 for s in syms}   # valor arbitrario >0
        try:
            val = float(diff.subs(subs))
            return val > 0
        except Exception:
            pass
    # último recurso: comparación por cadena (determinista, no matemática)
    return str(a) > str(b)

def leading_coeff(terms, poly_val):
    """
    Dado un iterable de términos (sumandos), devuelve (coeficiente_líder, exponente_líder)
    relativas a poly_val, aun cuando los exponentes sean simbólicos.
    """
    max_coeff = None
    max_exp = None

    for t in terms:
        # expandir y "aplanar" potencias anidadas
        t_flat = flatten_pow(sp.expand(t))
        coeff, exp = t_flat.as_coeff_exponent(poly_val)

        # si el término no contiene poly_val, as_coeff_exponent devuelve (t_flat, 0)
        # en ese caso tratamos exp == 0 y coeff == t_flat (no es relevante si no hay poly_val)
        if max_coeff is None:
            max_coeff, max_exp = coeff, exp
            continue

        # comparar exponentes simbólicos con fallback
        if _symbolic_gt(exp, max_exp):
            max_coeff, max_exp = coeff, exp

    return sp.simplify(max_coeff), sp.simplify(max_exp)


def numeric_equiv(expr1, expr2, tol=1e-9, trials=5):
    """
    Verifica si dos expresiones simbólicas son equivalentes numéricamente
    evaluándolas en valores aleatorios en el rango (0,1).

    Parameters
    ----------
    expr1, expr2 : sympy.Expr
        Expresiones simbólicas a comparar.
    tol : float
        Tolerancia relativa para la comparación.
    trials : int
        Número de evaluaciones aleatorias.

    Returns
    -------
    bool
        True si las expresiones son equivalentes numéricamente dentro de la tolerancia.
    """
    syms = list((expr1 - expr2).free_symbols)
    if not syms:
        # No hay símbolos → comparar directamente
        return abs(float(expr1 - expr2)) <= tol

    for _ in range(trials):
        # Generar sustitución aleatoria en (0,1)
        subs = {s: random.random() for s in syms}
        try:
            v1 = complex(expr1.subs(subs).evalf())
            v2 = complex(expr2.subs(subs).evalf())
        except Exception:
            # Si falla (ej: división por cero), intentar otro set
            continue

        if v1 == v2 == 0:
            continue  # ambos cero → ok
        if abs(v1 - v2) > tol * max(1, abs(v1), abs(v2)):
            return False
    return True


def flatten_exp_expr(expr, poly_val):

    expr_terms = sp.Add.make_args(sp.expand(expr))

    flat_expr = sp.Rational(0)

    for t in expr_terms:
        # expandir y "aplanar" potencias anidadas
        this_term = flatten_pow(sp.expand(t))

        flat_expr += sp.powsimp(this_term)



    if sp.simplify( sp.expand()) == sp.Rational('0') :

        return(flat_expr)
    
    else:
        
        warnings.warn("No se pudo aplanar la expresión. Revisar", UserWarning)
        return(expr)
        

def simplify_n_monic2(tt, poly_val=s):

    
    num, den = sp.fraction(sp.together(tt))

    num = flatten_exp_expr(num, poly_val)
    den = flatten_exp_expr(den, poly_val)

    # detectar coeficientes principales sin forzar polinomio entero
    num_terms = sp.Add.make_args(sp.expand(num))
    den_terms = sp.Add.make_args(sp.expand(den))
    
    
    knum, num_exp = leading_coeff(num_terms, poly_val)
    kden, den_exp = leading_coeff(den_terms, poly_val)
    
    k = knum / kden
    num = sp.simplify(num / knum)
    den = sp.simplify(den / kden)

    return k, num, den


def simplify_n_monic(tt, poly_val = s):
    '''
    Factoriza una función racional tt, en polinmios numerador y denominador 
    mónicos multiplicados por un escalar k.


    Parameters
    ----------
    tt : Expr
        Polinomio de fracciones a simplificar.


    Returns
    -------
    k
        escala o factor de la función tt.
    num
        Polinomio numerador simplificado en forma monica.
    den
        Polinomio denominador simplificado en forma monica.


    Raises
    ------
    ValueError
        Si la entrada no es una expresión simbólica.


    See Also
    --------
    :func:`print_latex`
    :func:`to_latex`
    :func:`a_equal_b_latex_s`


    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.general import s, simplify_n_monic
    >>> tt = (s**2 + 3*s + 2) / (2*s**2 + 5*s + 3)
    >>> k, num, den = simplify_n_monic(tt)
    >>> simplified_tt = k * num / den
    >>> print(simplified_tt)
    (s + 2)/(2*s + 3)
    
    '''
    if not isinstance(tt, sp.Expr):
        raise ValueError("La entrada debe ser una expresión simbólica.")

    if not isinstance(poly_val, sp.Symbol):
        raise ValueError("La variable del polinomio debe ser un símbolo SymPy.")

    # Obtener el numerador y el denominador de la expresión y convertirlos en polinomios
    num, den = sp.fraction(sp.simplify(sp.expand(tt)))
    
    if num.has(poly_val):

        num = sp.poly(num, poly_val)
        
        knum = num.LC() 

        # Convertir el numerador y el denominador a forma monica
        num = num.monic()

        num = num.expr
        
    else:
        
        knum = sp.Rational(1)
        
    if den.has(poly_val):

        den = sp.poly(den, poly_val)
        
        kden = den.LC() 

        # Convertir el numerador y el denominador a forma monica
        den = den.monic()
        
        den = den.expr
        
    else:
        
        kden = sp.Rational(1)
        
    
    # Calcular el coeficiente principal del numerador y el denominador
    k = knum / kden
    
    # Devolver el polinomio simplificado en forma monica
    return(k, num, den)

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
    
    if not isinstance(nn, Integral) or nn < 0:
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
    ValueError
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

    if not isinstance(a, (sp.Expr, str)):
        raise ValueError("a debe ser un símbolo o una cadena.")
    
    if not isinstance(b, (str, list, sp.Expr, sp.MatrixBase)):
        raise ValueError("b debe ser un símbolo o una lista de símbolos.")
    
    a_str = sp.latex(a) if isinstance(a, sp.Expr) else a
    b_str = sp.latex(b) if isinstance(b, (sp.Expr, sp.MatrixBase, list)) else b
    
    return '$' + a_str + '=' + b_str + '$'

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
    ValueError
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
    >>> print_latex(expr_simb_expr('tt', tt1, '\\neq'))  
    [LaTex formated equation]
    >>> print_latex(expr_simb_expr('tt', tt))
    [LaTex formated equation]


    '''

    if not isinstance(a, (sp.Expr, str)):
        raise ValueError("a debe ser un símbolo o una cadena.")
    
    if not isinstance(b, (list, sp.Expr)):
        raise ValueError("b debe ser un símbolo o una lista de símbolos.")
    
    if not isinstance(symbol, str):
        raise ValueError("symbol debe ser un string que represente un comando interpretable en LaTex.")
    
    a_str = sp.latex(a) if isinstance(a, sp.Expr) else a
    
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
    ValueError
        Si unsimbolo no es un símbolo ni una cadena.


    See Also
    --------
    :func:`print_latex`
    :func:`to_latex`
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
        raise ValueError("unsimbolo debe ser un símbolo o una cadena.")
    
    a_str = sp.latex(unsimbolo) if isinstance(unsimbolo, sp.Expr) else unsimbolo

    
    return '$' + a_str + '$'

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
    ValueError
        Si unstr no es una cadena.


    See Also
    --------
    :func:`print_subtitle`
    :func:`to_latex`
    :func:`to_latex`


    Examples
    --------
    >>> import sympy as sp
    >>> from pytc2.general import to_latex, print_latex
    >>> print(to_latex('x'))
    $x$
    >>> print_latex(to_latex('x'))  
    [LaTex formated equation]

    '''
    if not isinstance(unstr, str):
        raise ValueError("unstr debe ser una cadena.")
    

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
    ValueError
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
        raise ValueError("unstr debe ser una cadena.")

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
    ValueError
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
        raise ValueError("unstr debe ser una cadena.")

    # unstr = unstr + '\n'
    # unstr1 =  '-' * (len(unstr)-1) + '\n' 
    # print( '\n\n' + unstr + unstr1 )

    print( '\n\n\033[4m' + unstr + '\033[0m' )


def running_in_spyder():
    import sys
    return "spyder_kernels.console" in sys.modules


def running_in_jupyter():
    try:
        from IPython import get_ipython
        shell = get_ipython().__class__.__name__
        return shell == "ZMQInteractiveShell"  # típico de Jupyter Notebook/Lab
    except:
        return False

     
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
    ValueError
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
        raise ValueError("unstr debe ser una cadena.")

    if running_in_jupyter():
    
        display(Markdown('#### ' + unstr))
        
    else:
    
        print( '\n\n\033[4m' + unstr + '\033[0m' )
        
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
      ValueError: Si at_en_db no es de tipo `float`.


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

    if not isinstance(at_en_db, Real):
        raise ValueError('at_en_db debe ser float')

    return at_en_db / (20. * np.log10(small_val+np.exp(1.)))

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
      ValueError: Si at_en_db no es de tipo `float`.
    

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

    if not isinstance(at_en_np, Real):
        raise ValueError('at_en_np debe ser float')

    return at_en_np * (20. * np.log10(small_val+np.exp(1.)))
    