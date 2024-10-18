import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqz
from pytc2.filtros_digitales import fir_design_pm, fir_design_ls

# Parámetros comunes del transformador Hilbert FIR
Ftype = 'h'  # Transformador de Hilbert
order = 70  # Orden del filtro (debe ser par para tipo III)
fs = 2.0  # Frecuencia de muestreo
lgrid = 16  # Densidad del grid

# Bandas y especificaciones del transformador Hilbert
band_edges = [0.05, 0.95]  # Bandas de diseño (en términos de frecuencia normalizada)
desired = [1., 1.]  # Respuesta ideal en magnitud
W = [1., 1.]  # Peso en las bandas

# Diseño del filtro con PM
b_pm, _, _ = fir_design_pm(order, band_edges, desired, grid_density=lgrid, fs=fs, filter_type=Ftype)

# Diseño del filtro con LS
b_ls = fir_design_ls(order, band_edges, desired, grid_density=lgrid, fs=fs, filter_type=Ftype)

# Comparar la respuesta en frecuencia de ambos filtros
w, h_pm = freqz(b_pm, worN=1024)
w, h_ls = freqz(b_ls, worN=1024)

# Graficar ambas respuestas en magnitud
plt.figure(figsize=(10, 6))
plt.plot(w/np.pi, np.abs(h_pm), label='PM - Parks-McClellan', linestyle='--')
plt.plot(w/np.pi, np.abs(h_ls), label='LS - Mínimos Cuadrados', linestyle='-')
plt.title('Comparación de Transformadores de Hilbert FIR - PM vs LS (Magnitud)')
plt.xlabel('Frecuencia Normalizada (×π rad/sample)')
plt.ylabel('Amplitud')
plt.legend()
plt.grid()
plt.show()

# Graficar la respuesta en fase
plt.figure(figsize=(10, 6))
plt.plot(w/np.pi, np.angle(h_pm) * 180 / np.pi, label='PM - Parks-McClellan', linestyle='--')
plt.plot(w/np.pi, np.angle(h_ls) * 180 / np.pi, label='LS - Mínimos Cuadrados', linestyle='-')
plt.title('Comparación de Transformadores de Hilbert FIR - PM vs LS (Fase)')
plt.xlabel('Frecuencia Normalizada (×π rad/sample)')
plt.ylabel('Fase (grados)')
plt.legend()
plt.grid()
plt.show()

# Análisis de la planicidad de la banda de paso
passband_pm = np.abs(h_pm[(w >= band_edges[0] * np.pi) & (w <= band_edges[1] * np.pi)])
passband_ls = np.abs(h_ls[(w >= band_edges[0] * np.pi) & (w <= band_edges[1] * np.pi)])

ripple_pm = np.max(passband_pm) - np.min(passband_pm)
ripple_ls = np.max(passband_ls) - np.min(passband_ls)

print(f'Ripple en la banda de paso - PM: {ripple_pm:.4e}')
print(f'Ripple en la banda de paso - LS: {ripple_ls:.4e}')
