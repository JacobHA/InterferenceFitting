import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
from scipy.signal import savgol_filter

class SlitData:
    def __init__(self, file, num_slits, a, d, screen_dist=90*10, wavelength = 632.8e-9 * 1e3, addnl_offset=0):
        self.file = str(file) + str('.txt')
        self.num_slits = num_slits
        self.a = a
        self.d = d
        self.screen_dist = screen_dist
        self.wavelength = wavelength
        self.addnl_offset = addnl_offset

        self.a_fit = None
        self.d_fit = None
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
                self.position.append(float(pos))
                self.intensity.append(float(inten.split('\n')[0]))

    def clean_data(self):
        # self.center = np.mean(self.position)
        self.position = np.array(self.position)
        self.intensity = np.array(self.intensity)
        self.center = np.argmax(self.intensity)
        self.position -= self.position[self.center]
        self.position += self.addnl_offset
        self.intensity -= min(self.intensity)
        self.intensity = self.intensity/100
        self.fit_position = np.linspace(min(self.position), max(self.position), 10_000)
    

    def generate_fit(self):
        if self.num_slits == 1:
            self.fit_intensity = self.fit_single_slit()
            self.I0_fit, self.d_fit = self.popt
        else:
            self.fit_intensity = self.fit_N_slit()
            self.I0_fit, self.a_fit, self.d_fit = self.popt
        
        
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
        self.popt, self.pcov = curve_fit(self.single_slit, self.position, self.intensity, [0.05, self.d])
        return self.single_slit(self.fit_position, *self.popt)

    def fit_N_slit(self):
        self.popt, self.pcov = curve_fit(self.N_slit, self.position, self.intensity, [0.05, self.a, self.d], bounds=([0,0,0],[0.5,1,0.1]))
        return self.N_slit(self.fit_position, *self.popt)
   
    def generate_plot(self):
        self.generate_fit()
        plt.figure()
        plt.xlabel('Position (mm)', fontsize=16)
        plt.ylabel('Intensity (a.u.)', fontsize=16)
        if self.num_slits == 1:
            plt.title(f'{int(self.num_slits)} slit, Width: {self.d}mm', fontsize=20)
        else:
            plt.title(f'{int(self.num_slits)} slits, Width: {self.d}mm, Spacing: {self.a}mm', fontsize=20)
        plt.plot(self.position, self.intensity, label = 'Data')
        plt.plot(self.fit_position, self.fit_intensity, label = 'Fit')
        plt.legend(fontsize=16)
        plt.show()

    def report_params(self):
        self.generate_fit()
        self.perr = np.sqrt(np.diag(self.pcov))
        print('=========')
        print('I0:', round(self.I0_fit,4))
        print('I0 err', round(self.perr[0],4))
        if self.num_slits > 1:
            # slit spacing only valid for multiple slits
            print('=========')
            print('a true:', round(self.a,4))
            print('A:', round(self.a_fit,4))
            print('A err:', round(self.perr[1],4))
            print('=========')
            print('d true:', round(self.d,4))
            print('D:', round(self.d_fit,4))
            print('D err:', round(self.perr[2],4))
        else:
            print('=========')
            print('d true: ', round(self.d,4))
            print('D:', round(self.d_fit,4))
            print('D err:', round(self.perr[1],4))


def linear_function(x, a):
    return a * x 

    
# x = SlitData('2 slits a0_08d_50', 2, 0.5, 0.08, addnl_offset=0.01)
# doubleslit = SlitData('2 slits a0_04d_50', 2, 0.5, 0.08, addnl_offset=0.01)
# peaks = [-6.495, -5.1326, -3.5253, -2.529, -1.224, 0.003, 1.183, 2.446, 3.824, 5.165, 6.511]
# modes = np.arange(-5,6,1)
# slope,pcov = curve_fit(linear_function, modes, peaks)
# slope=slope[0]
# plt.figure()
# plt.plot(peaks, 'bo', label='data')
# plt.plot(linear_function(modes, slope), 'r-', label=f'fit\nInferred wavelength: {round(slope*1e-3 * x.a * 1e-3 / (.9)*1e9,5)} nm')
# plt.legend()

# x = SlitData('2 slits a0_08d_25', 2, 0.25, 0.08, addnl_offset=0.01)
# peaks = [-10.4, -7.44, -5.02, -2.488, 0.01, 2.371, 4.92, 7.32, 10.41]
# modes = np.arange(-4,5,1)
# print(modes)
# slope,pcov = curve_fit(linear_function, modes, peaks)
# print(slope)

# plt.figure()
# plt.plot(peaks, 'bo', label='data')
# plt.plot(linear_function(modes, slope), 'r-', label=f'fit\nInferred wavelength: {round(np.mean(np.diff(peaks))*1e-3 * x.a * 1e-3 / (1.)*1e9,5)} nm')
# plt.legend()


# x = SlitData('2 slits a0_04d_50', 2, 0.5, 0.04, addnl_offset=1.2)
# peaks = [-7.95, -6.61, -3.94, -2.642, -1.289, -0.001, 1.198, 2.405, 3.704, 6.42, 7.71]
# peaks=peaks[2:-2]
# modes = np.arange(-3,4,1)
# print(modes)
# slope,pcov = curve_fit(linear_function, modes, peaks)
# print(slope)

# plt.figure()
# plt.plot(peaks, 'bo', label='data')
# plt.plot(linear_function(modes, slope), 'r-', label=f'fit\nInferred wavelength: {round(np.mean(np.diff(peaks))*1e-3 * x.a * 1e-3 / (1.)*1e9,5)} nm')
# plt.legend()

# x = SlitData('2 slits a0_04d_25', 2, 0.25, 0.04, addnl_offset=0.01)
# peaks = [-7.823, -5.15, -2.543, -0.014, 2.340,  5.0, 7.616]
# modes = np.arange(-3,4,1)
# print(modes)
# slope,pcov = curve_fit(linear_function, modes, peaks)
# print(slope)

# plt.figure()
# plt.plot(peaks, 'bo', label='data')
# plt.plot(linear_function(modes, slope), 'r-', label=f'fit\nInferred wavelength: {round(np.mean(np.diff(peaks))*1e-3 * x.a * 1e-3 / ()*1e9,5)} nm')
# plt.legend()


x = SlitData('Single W = 0.16mm', 1, 0.25, 0.16, addnl_offset=-0.1, screen_dist=110*10)
peaks = [-13.37, -9.54, -5.51, -0.095, 5.64, 9.68, 13.4]
modes = np.arange(-3,4,1)
slope,pcov = curve_fit(linear_function, modes, peaks)
slope=slope[0]
plt.figure()
plt.title(r'Mode Spacing for W=0.16 mm', fontsize=20)
plt.plot(modes, peaks, 'bo', label='Mode gap from intensity data')
plt.plot(modes, linear_function(modes, slope), 'r-', label=f'Linear fit\nInferred wavelength: {round(slope*1e-3 * x.d * 1e-3 / (1.1)*1e9,1)} nm')
plt.legend(fontsize=12)
plt.xlabel('Mode Number', fontsize=16)
plt.ylabel('Mode Spacing (mm)', fontsize=16)
plt.show()
factor = 1e-3 * x.d * 1e-3 / (1.1)*1e9
print(f'{slope*factor} +/- {factor*np.diag(np.sqrt(pcov))[0]}') 


# x = SlitData('Single W = .08mm', 1, 0.25, 0.08, addnl_offset=-0.7)
# # singleslit = SlitData('Single W = .04mm', 1, 0.25, 0.08, addnl_offset=0.01)
# peaks = [-18.7, -10.6, 0.8, 11.92, 19.87]
# modes = np.arange(-2,3,1)
# slope,pcov = curve_fit(linear_function, modes, peaks)
# slope=slope[0]
# plt.figure()
# plt.title('Mode Spacing for W=0.08mm', fontsize=20)
# plt.plot(modes, peaks, 'bo', label='Mode gap from intensity data')
# plt.plot(modes, linear_function(modes, slope), 'r-', label=f'Linear fit\nInferred wavelength: {round(slope*1e-3 * x.d * 1e-3 / (1.1)*1e9,5)} nm')
# plt.legend(fontsize=12)
# plt.xlabel('Mode Number', fontsize=16)
# plt.ylabel('Mode Spacing (mm)', fontsize=16)
# plt.show()
# factor = 1e-3 * x.d * 1e-3 / (1.1)*1e9
# print(f'{slope*factor} +/- {factor*np.diag(np.sqrt(pcov))[0]}') 



# x = SlitData('Single W = .04mm', 1, 0.25, 0.04, addnl_offset=0.01)
# peaks = [-23, 0.7, 24.12]
# modes = np.arange(-1,2,1)
# slope,pcov = curve_fit(linear_function, modes, peaks)
# slope=slope[0]
# plt.figure()
# plt.plot(modes, peaks, 'bo', label='Mode gap from intensity data')
# plt.plot(modes, linear_function(modes, slope), 'r-', label=f'Linear fit\nInferred wavelength: {round(slope*1e-3 * x.d * 1e-3 / (1.1)*1e9,1)} nm')
# plt.title(r'Mode Spacing for W=0.04 mm', fontsize=20)
# plt.legend(fontsize=12)
# plt.xlabel('Mode Number', fontsize=16)
# plt.ylabel('Mode Spacing (mm)', fontsize=16)
# plt.show()
# factor = 1e-3 * x.d * 1e-3 / (1.1)*1e9
# print(f'{slope*factor} +/- {factor*np.diag(np.sqrt(pcov))[0]}') 



# x = SlitData('4 slits 0.125mm separation', 4, 0.125, 0.04, addnl_offset=0.035)
# quadslit = SlitData('4 slits 0.125mm separation', 4, 0.125, 0.04, addnl_offset=0.035)


x.generate_plot()
x.report_params()


# plt.figure()
# plt.title('Different Slits - Normalized Interference Patterns', fontsize=20)
# for exp in [singleslit, doubleslit, quadslit]:
#     exp.generate_fit()
#     # plt.plot(exp.position, exp.intensity, label=f'{exp.num_slits} Slit(s) Data')
#     plt.plot(exp.fit_position, exp.fit_intensity/max(exp.fit_intensity), label=f'{exp.num_slits} Slit(s) Fit')

# plt.legend()
# plt.show()