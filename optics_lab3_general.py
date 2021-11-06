import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# TODO: add peak detection to infer wavelength
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
from scipy.signal import savgol_filter

class SlitData:
    """
    Gathers and cleans position, intensity data. Plots and generates fits of data.
    Assumes data stored in a .txt file, one column of position, one column of intensity data.
    Use position_unit_conversion to convert from your units (e.g. mm) to m.
    Ammend the get_data method at will in order to properly scrape data.
    All units are in SI (m) 
    Infinitesimal addnl_offset used for non-divergence at zero position
    
    """
    def __init__(self, file, num_slits=0, slit_spacing=0, slit_width=0, screen_dist=1.0, wavelength = 632.8e-9, addnl_offset=1e-9, position_unit_conversion=1e-3):
        self.file = str(file) + str('.txt')
        self.num_slits = num_slits
        self.slit_spacing = slit_spacing
        self.slit_width = slit_width
        self.screen_dist = screen_dist
        self.wavelength = wavelength
        self.addnl_offset = addnl_offset
        self.position_unit_conversion = position_unit_conversion

        self.slit_spacing_fit = None
        self.slit_width_fit = None
        self.I0_fit = None
        self.fit_position = None
        self.fit_intensity = None

        self.get_data()
        self.clean_data()
        

    def get_data(self):
        with open(self.file, 'r') as f:
            lines = f.readlines()
            self.position = []
            self.intensity = []
            for line in lines[1:]:
                pos, inten = line.split('\t')
                self.position.append(float(pos)*self.position_unit_conversion)
                self.intensity.append(float(inten.split('\n')[0]))

    def clean_data(self):
        # self.center = np.mean(self.position)
        self.position = np.array(self.position)
        self.intensity = np.array(self.intensity)
        self.center = np.argmax(self.intensity)
        self.position -= self.position[self.center]
        self.position += self.addnl_offset
        self.intensity -= min(self.intensity)
        self.intensity = self.intensity/100 # % of intensity data
        # self.intensity /= self.intensity.max() # can optionally normalize intensity
        self.fit_position = np.linspace(min(self.position), max(self.position), 10_000) # 10_000 points on x-axis of fit
    

    def generate_fit(self):
        if self.num_slits == 1:
            self.fit_intensity = self.fit_single_slit()
            self.I0_fit, self.slit_width_fit = self.popt
        else:
            self.fit_intensity = self.fit_N_slit()
            self.I0_fit, self.slit_spacing_fit, self.slit_width_fit = self.popt
        
        
    def N_slit(self, pos, I0, A, D):
        theta = np.arctan((pos)/self.screen_dist)
        beta = np.pi / self.wavelength * D * np.sin(theta)
        gamma = np.pi / self.wavelength * A * np.sin(theta)
        return I0 * (np.sin(beta)/beta * np.sin(self.num_slits * gamma)/np.sin(gamma))**2

    def single_slit(self, pos, I0, D):
        theta = np.arctan((pos)/self.screen_dist)
        beta = np.pi / self.wavelength * D * np.sin(theta)
        return I0*(np.sin(beta)/beta)**2

    def fit_single_slit(self):
        self.popt, self.pcov = curve_fit(self.single_slit, self.position, self.intensity, [0.05, self.slit_width], bounds=([0,0],[1,0.5])) # guesses and bounds provided to help the fit converge.
        return self.single_slit(self.fit_position, *self.popt)

    def fit_N_slit(self):
        self.popt, self.pcov = curve_fit(self.N_slit, self.position, self.intensity, [0.05, self.slit_spacing, self.slit_width], bounds=([0,0,0],[0.5,1,0.1]))
        return self.N_slit(self.fit_position, *self.popt)
   
    def generate_plot(self):
        self.generate_fit()
        plt.figure()
        plt.xlabel('Position (m)', fontsize=16)
        plt.ylabel('Intensity (a.u.)', fontsize=16)
        if self.num_slits == 1:
            plt.title(f'{int(self.num_slits)} slit, Width: {self.slit_width}m', fontsize=20)
        else:
            plt.title(f'{int(self.num_slits)} slits, Width: {self.slit_width}m, Spacing: {self.slit_spacing}m', fontsize=20)
        plt.plot(self.position, self.intensity, label = 'Data')
        plt.plot(self.fit_position, self.fit_intensity, label = 'Fit')
        plt.legend(fontsize=16)
        plt.show()

    def report_params(self):
        self.generate_fit()
        self.perr = np.sqrt(np.diag(self.pcov))
        print('=========')
        print('Observed peak intensity:', self.intensity.max())
        print('Fitted peak intensity:', self.I0_fit*self.num_slits**2, '+/-', self.perr[0])
        # limit as pos -> 0  [sin(N gamma) / sin(gamma)]^2 = N^2
        if self.num_slits > 1:
            # slit spacing only valid for multiple slits
            print('=========')
            print('True slit spacing:', self.slit_spacing)
            print('Fitted slit spacing:', self.slit_spacing_fit, '+/-', self.perr[1])
            print('=========')
            print('True slit width:', self.slit_width)
            print('Fitted slit width:', self.slit_width_fit, '+/-', self.perr[2])
        else:
            print('=========')
            print('True slit width: ', self.slit_width)
            print('Fitted slit width:', self.slit_width_fit, '+/-', self.perr[1])
        print('=========')

x = SlitData('2 slits a0_04d_25', num_slits=2, slit_spacing=0.25e-3, slit_width=0.04e-3, addnl_offset=0.0001)

x.generate_plot()
x.report_params()
