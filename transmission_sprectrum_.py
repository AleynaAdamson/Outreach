
import ipywidgets as widgets
from ipywidgets import HBox, VBox
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import FloatLogSlider, HBox, VBox, interactive_output, FloatSlider
from IPython.display import display, clear_output
from ipywidgets import Layout
import random

#simulated spectrum
def part_1():
	random.seed(1)
	noise_level = 0.003
	
	# Convert log-abundance (VMR) → feature height (VMR = (Volume of a given gas) / (Total volume of the gas mixture))
	def abundance_to_depth(log_abundance, scale=0.08):
	    return scale * (log_abundance - np.log10(1e-10)) / (np.log10(1) - np.log10(1e-10))
	
	# Molecules with center wavelengths (µm), sigma, and scaling factor
	MOLECULES = {
	    'H₂O': [(1.4, 0.05, 1.4), (1.97, 0.05, 1.0), (2.57, 0.07, 1.5), (1.85, 0.05, 0.97), (0.95, 0.045, 0.45), (1.17, 0.045, 0.45), (3.35, 0.06, 0.9), (3.2, 0.06, 0.9), (3, 0.1, 1)],
	    'CH₄': [(2.3, 0.1, 0.2), (2.5, 0.01, 0.3), (3.5, 0.2, 1.05), (3.8, 0.03, 0.3)],
	
	    'CO₂': [(2.7, 0.04, 2.1), (2.83, 0.07, 2.25), (4.38, 0.1, 3.7)],
	    'O2': [(0.1, 0.01, 0.3), (0.7, 0.01, 0.8)],
	    'Na': [(0.61, 0.03, 0.9)],
	    # 'Na/K': [(0.59, 0.01, 0.3), (0.77, 0.01, 0.3)],
	    'K': [(0.76, 0.01, 0.8)],
	    'SO2': [(4.03, 0.05, 0.9)],
	    "CO": [(4.79, 0.3, 1.5), (2.36, 0.1, 0.8)],
	    'O3': [(4.79, 0.1, 2)]
	    # SO₂ band at 7.3 µm (outside current plot)
	}
	
	# Generate spectrum using Gaussian absorption for each molecule
	def generate_spectrum(log_vmr):
	    wavelengths = np.linspace(0.3, 5.0, 2500)
	    spectrum = np.ones_like(wavelengths) + noise_level * np.random.normal(size=len(wavelengths))
	    spectrum_no_noise = np.ones_like(wavelengths)
	
	    for mol, log_ab in log_vmr.items():
	        depth = abundance_to_depth(log_ab, scale=0.08)
	        for center, sigma, rel in MOLECULES[mol]:
	            gauss = rel * depth * np.exp(-((wavelengths - center)**2) / (2 * sigma**2))
	            spectrum    -= gauss
	            spectrum_no_noise -= gauss
	
	    return wavelengths, spectrum, spectrum_no_noise
	
	# Log abundance sliders
	sliders = {
	    mol: FloatSlider(value=0, min=0, max=1, step=0.1, description=f'{mol}') #VMR:
	    for mol in ['H₂O', 'CH₄', 'CO₂', 'O2', 'Na', 'K', 'SO2', 'CO', 'O3'] #'O2']
	}
	
	def update_plot(**vals):
	    log_vmr = {}
	    for mol, val in vals.items():
	        if val == 0:
	            log_vmr[mol] = -10  # minimum log abundance for zero slider value
	        else:
	            log_vmr[mol] = np.log10(val)
	
	    wl, spec, spec0 = generate_spectrum(log_vmr)
	
	    clear_output(wait=True)
	
	    plt.figure(figsize=(8, 5))
	    grad = np.linspace(0, 1, 500).reshape(1, -1)
	    plt.imshow(grad, extent=[wl.min(), wl.max(), (1-spec).min(), 0.35],
	               aspect='auto', cmap='nipy_spectral', alpha=0.2)
	
	    '''
	    for mol, peaks in MOLECULES.items():
	        for c, sigma, _ in peaks:
	            if 0.3 < c < 5.0:
	                color_map = {
	                    'H₂O': 'blue',
	                    'CH₄': 'forestgreen',
	                    'CO₂': 'red',
	                    'O2': 'darkviolet',
	                    'Na': 'yellow',
	                    'K': 'orange',
	                    'SO2': 'hotpink',
	                    'CO': 'gray',
	                    'O3': 'brown'
	                }
	                plt.axvspan(c - sigma, c + sigma, alpha=0.2, label=mol, color=color_map[mol])
	    '''
	    plt.scatter(wl, 1-spec, s=3, c='k', label='Noisy')
	    plt.plot(wl, 1-spec0, c='red', lw=1.5, label='Noiseless')
	    plt.xlabel('Wavelength (µm)')
	    plt.ylabel('Transit Depth')
	    #plt.title('Simulated Transit Spectrum — Scientific Molecules')
	    plt.tight_layout()
	    plt.show()
	
	slider_box = VBox([sliders[mol] for mol in sliders],layout=Layout(margin='100px 0px 0 0'))
	out = interactive_output(update_plot, sliders)
	layout = HBox([out, slider_box],layout=Layout(gap='0px', padding='0px', margin='0px') )
	return display(layout)
#different molecules
def part_2(csv):
	random.seed(1)
	noise_level = 0.003
	
	# Convert log-abundance (VMR) → feature height (VMR = (Volume of a given gas) / (Total volume of the gas mixture))
	def abundance_to_depth(log_abundance, scale=0.08):
	    return scale * (log_abundance - np.log10(1e-10)) / (np.log10(1) - np.log10(1e-10))
	
	# Molecules with center wavelengths (µm), sigma, and scaling factor
	MOLECULES = {
	    'H₂O': [(1.4, 0.05, 1.4), (1.97, 0.05, 1.0), (2.57, 0.07, 1.5), (1.85, 0.05, 0.97), (0.95, 0.045, 0.45), (1.17, 0.045, 0.45), (3.35, 0.06, 0.9), (3.2, 0.06, 0.9), (3, 0.1, 1)],
	    'CH₄': [(2.3, 0.1, 0.2), (2.5, 0.01, 0.3), (3.5, 0.2, 1.05), (3.8, 0.03, 0.3)],
	
	    'CO₂': [(2.7, 0.04, 2.1), (2.83, 0.07, 2.25), (4.38, 0.1, 3.7)],
	    'O2': [(0.1, 0.01, 0.3), (0.7, 0.01, 0.8)],
	    'Na': [(0.61, 0.03, 0.9)],
	    # 'Na/K': [(0.59, 0.01, 0.3), (0.77, 0.01, 0.3)],
	    'K': [(0.76, 0.01, 0.8)],
	    'SO2': [(4.03, 0.05, 0.9)],
	    "CO": [(4.79, 0.3, 1.5), (2.36, 0.1, 0.8)],
	    'O3': [(4.79, 0.1, 2)]
	    # SO₂ band at 7.3 µm (outside current plot)
	}
	
	# Generate spectrum using Gaussian absorption for each molecule
	def generate_spectrum(log_vmr):
	    wavelengths = np.linspace(0.3, 5.0, 2500)
	    spectrum = np.ones_like(wavelengths) + noise_level * np.random.normal(size=len(wavelengths))
	    spectrum_no_noise = np.ones_like(wavelengths)
	
	    for mol, log_ab in log_vmr.items():
	        depth = abundance_to_depth(log_ab, scale=0.08)
	        for center, sigma, rel in MOLECULES[mol]:
	            gauss = rel * depth * np.exp(-((wavelengths - center)**2) / (2 * sigma**2))
	            spectrum    -= gauss
	            spectrum_no_noise -= gauss
	
	    return wavelengths, spectrum, spectrum_no_noise
	
	# Log abundance sliders
	sliders = {
	    mol: FloatSlider(value=0, min=0, max=1, step=0.1, description=f'{mol}')
	    for mol in ['H₂O', 'CH₄', 'CO₂', 'O2', 'Na', 'K', 'SO2', 'CO', 'O3'] #'O2']
	}
	
	def update_plot(**vals):
	    log_vmr = {}
	    for mol, val in vals.items():
	        if val == 0:
	            log_vmr[mol] = -10  # minimum log abundance for zero slider value
	        else:
	            log_vmr[mol] = np.log10(val)
	
	    wl, spec, spec0 = generate_spectrum(log_vmr)
	
	    clear_output(wait=True)
	
	    plt.figure(figsize=(8, 5))
	    norm = (csv[:, 1] - np.min(csv[:, 1])) / (np.max(csv[:, 1])-np.min(csv[:, 1]))
	    plt.scatter(csv[:, 0], norm/3, color='purple')
	    grad = np.linspace(0, 1, 500).reshape(1, -1)
	    plt.imshow(grad, extent=[wl.min(), wl.max(), (1-spec).min(), 0.35],
	               aspect='auto', cmap='nipy_spectral', alpha=0.2)
	    '''
	    for mol, peaks in MOLECULES.items():
	        for c, sigma, _ in peaks:
	            if 0.3 < c < 5.0:
	                color_map = {
	                    'H₂O': 'blue',
	                    'CH₄': 'forestgreen',
	                    'CO₂': 'red',
	                    'O2': 'darkviolet',
	                    'Na': 'yellow',
	                    'K': 'orange',
	                    'SO2': 'hotpink',
	                    'CO': 'gray',
	                    'O3': 'brown'
	                }
	                plt.axvspan(c - sigma, c + sigma, alpha=0.2, label=mol, color=color_map[mol])
	      '''
	    plt.scatter(wl, 1-spec, s=3, c='k', label='Noisy')
	    plt.plot(wl, 1-spec0, c='red', lw=1.5, label='Noiseless')
	    plt.xlabel('Wavelength (µm)')
	    plt.ylabel('Transit Depth')
	    #plt.title('Simulated Transit Spectrum — Scientific Molecules')
	    plt.xlim(0.5, 5)
	    #plt.legend(loc='upper right', ncol=2)
	    #plt.grid(True)
	    plt.tight_layout()
	    plt.show()
	
	slider_box = VBox([sliders[mol] for mol in sliders],layout=Layout(margin='100px 0px 0 0'))
	out = interactive_output(update_plot, sliders)
	layout = HBox([out, slider_box],layout=Layout(gap='0px', padding='0px', margin='0px') )
	return display(layout)
  
