#--------------------------------------------------------------------------------------------------------------#

import sympy as sp
#from .general import s
s = sp.symbols( 's', complex = True )

#--------------------------------------------------------------------------------------------------------------#

def remover_polo_dc( imit, omega_zero = None, isSigma = False ):
    '''
    Se removerá el residuo en continua (s=0) de la imitancia ($I$) de forma 
    completa, o parcial en el caso que se especifique una omega_zero. 
    Como resultado de la remoción, quedará otra función racional definida
    como:
        
    $$ I_R = I - k_0/s  $$
    
    siendo 

    $$ k_0=\lim\limits _{s\to0}I.s $$
    
    En cuanto se especifique omega_zero, la remoción parcial estará definida 
    como, siempre que isSigma = False

    $$ I_{R}\biggr\rfloor_{s^{2}=-\omega_z^{2}}=0=I-s.k_{0}^{'}\biggr\rfloor_{s^{2}=-\omega_z^{2}} $$
    
    siendo 
    
    $$ k_{0}^{'}=I.s\biggr\rfloor_{s^{2}=-\omega_z^{2}} $$

    De lo contrario, con isSigma = True

    $$ I_{R}\biggr\rfloor_{s=-\omega_z}=0=I-s.k_{0}^{'}\biggr\rfloor_{s=-\omega_z} $$
    
    siendo 
    
    $$ k_{0}^{'}=I.s\biggr\rfloor_{s=-\omega_z} $$
    

    Parameters
    ----------
    imit : Symbolic
        Imitancia que se utilizará para la remoción. Es una función racional 
        simbólica que tendrá un polo de orden 1 en 0, es decir la 
        diferencia de grados entre num y den será exactamente -1.
    omega_zero : Symbolic
        Frecuencia a la que la imitancia será cero luego de la remoción.

    Returns
    -------
    imit_r : Symbolic
        Imitancia luego de la remoción
    k_inf : Symbolic
        Valor del residuo en infinito
    '''

    if omega_zero is None:
        # remoción total
        k_cero = sp.limit(imit*s, s, 0)
        
    else:
        # remoción parcial en el eje j\omega
    	if isSigma is False:
	        k_cero = sp.simplify(sp.expand(imit*s)).subs(s**2, -(omega_zero**2) )

    	# remoción parcial en el eje \sigma
    	else:
	        k_cero = sp.simplify(sp.expand(imit*s)).subs(s, -omega_zero )

    k_cero = k_cero/s
    
    # extraigo C3
    imit_r = sp.factor(sp.simplify(sp.expand(imit - k_cero)))

    return( [imit_r, k_cero] )


#--------------------------------------------------------------------------------------------------------------#

def remover_polo_infinito( imit, omega_zero = None, isSigma = False ):
    '''
    Se removerá el residuo en infinito de la imitancia ($I$) de forma 
    completa, o parcial en el caso que se especifique una omega_zero. 
    Como resultado de la remoción, quedará otra función racional definida
    como:
        
    $$ I_R = I - s.k_\infty  $$
    
    siendo 

    $$ k_{\infty}=\lim\limits _{s\to\infty}I.\nicefrac{1}{s} $$
    
    En cuanto se especifique omega_zero, la remoción parcial estará definida 
    como, siempre que isSigma = False

    $$ I_{R}\biggr\rfloor_{s^{2}=-\omega_z^{2}}=0=I-s.k_{\infty}^{'}\biggr\rfloor_{s^{2}=-\omega_z^{2}} $$
    
    siendo 
    
    $$ k_{\infty}^{'}=I.\nicefrac{1}{s}\biggr\rfloor_{s^{2}=-\omega_z^{2}} $$

    De lo contrario, con isSigma = True

    $$ I_{R}\biggr\rfloor_{s=-\omega_z}=0=I-s.k_{\infty}^{'}\biggr\rfloor_{s=-\omega_z} $$
    
    siendo 
    
    $$ k_{\infty}^{'}=I.\nicefrac{1}{s}\biggr\rfloor_{s=-\omega_z} $$
    

    Parameters
    ----------
    imit : Symbolic
        Imitancia que se utilizará para la remoción. Es una función racional 
        simbólica que tendrá un polo de orden 1 en infinito, es decir la 
        diferencia de grados entre num y den será exactamente 1.
    omega_zero : Symbolic
        Frecuencia a la que la imitancia será cero luego de la remoción.
    isSigma : Boolean
	Indica si la remoción parcial se realiza en el eje j\omega ( isSigma = False )
	o bien sobre el eje \sigma ( isSigma = True )

    Returns
    -------
    imit_r : Symbolic
        Imitancia luego de la remoción
    k_inf : Symbolic
        Valor del residuo en infinito
    '''

    if omega_zero is None:
        # remoción total
        k_inf = sp.limit(imit/s, s, sp.oo)
        
    else:
        # remoción parcial en el eje j\omega
        if isSigma is False:
        	k_inf = sp.simplify(sp.expand(imit/s)).subs(s**2, -(omega_zero**2) )
	
    	# remoción parcial en el eje \sigma
        else:
        	k_inf = sp.simplify(sp.expand(imit/s)).subs(s, -omega_zero )		

    k_inf = k_inf * s

    # extraigo C3
    imit_r = sp.factor(sp.simplify(sp.expand(imit - k_inf)))

    return( [imit_r, k_inf] )

#--------------------------------------------------------------------------------------------------------------#