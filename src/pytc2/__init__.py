# read version from installed package
from importlib.metadata import version
__version__ = version("pytc2")

# Aplanamos toda la estructura para facilitar su uso. 
# (Hay que pensarlo un poco más ...)
# from .sistemas_lineales import analyze_sys, pzmap, bodePlot, plot_plantilla, tf2sos_analog, pretty_print_SOS, pretty_print_lti
# from .cuadripolos import S2Tabcd_s, S2Ts_s, SparY_s, SparZ_s, T2Y, T2Y_s, T2Z, T2Z_s, Tabcd2S_s, TabcdLYZ, TabcdLYZ_s, TabcdLZY, TabcdLZY_s, TabcdY, TabcdY_s, TabcdZ, TabcdZ_s, Ts2S_s, Y2T, Y2T_s, Z2T, Z2T_s, calc_MAI_impedance_ij, calc_MAI_vtransf_ij_mn, calc_MAI_ztransf_ij_mn, display, may2y, print_console_alert, print_console_subtitle, print_latex, print_subtitle, y2mai 
# from .dibujar import dibujar_Pi, dibujar_Tee, dibujar_cauer_LC, dibujar_cauer_RC_RL, dibujar_elemento_derivacion, dibujar_elemento_serie, dibujar_espaciador, dibujar_espacio_derivacion, dibujar_foster_derivacion, dibujar_foster_serie, dibujar_funcion_exc_abajo, dibujar_funcion_exc_arriba, dibujar_puerto_entrada, dibujar_puerto_salida, dibujar_tanque_RC_derivacion, dibujar_tanque_RC_serie, dibujar_tanque_RL_derivacion, dibujar_tanque_RL_serie, dibujar_tanque_derivacion, dibujar_tanque_serie, display, str_to_latex, to_latex
# from .general import Chebyshev_polynomials, pp
# from .imagen import I2T, I2T_s, db2nepper, nepper2db
# from .ltspice import ltsp_capa_derivacion, ltsp_etiquetar_nodo, ltsp_ind_serie, ltsp_nuevo_circuito
# from .remociones import modsq2mod, modsq2mod_s, remover_polo_dc, remover_polo_infinito, remover_polo_jw, remover_polo_sigma, remover_valor, remover_valor_en_dc, remover_valor_en_infinito, simplify_n_monic, tanque_y, tanque_z, trim_func_s, trim_poly_s
# from .sintesis_dipolo import cauer_LC, cauer_RC, foster, remover_polo_dc, remover_polo_infinito, remover_valor_en_dc, remover_valor_en_infinito, trim_func_s
