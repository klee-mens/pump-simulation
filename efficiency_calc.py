# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 17:06:11 2021

@author: mens
"""

import numpy as np
from inversion_calc import inversion
import matplotlib.pyplot as plt
from init_param import param_struct

# =============================================================================
# Berechung
# =============================================================================
def efficency(param):
	I_p_Max = param.I_p
	I_p_vec = np.linspace(0, I_p_Max)
	all_effi = []
	for I_p in I_p_vec:
		param.I_p = I_p
		egal = inversion(param)
		all_effi.append(param.effi)

	all_effi = np.array(all_effi)
	param.effi_array = all_effi
	return all_effi


# =============================================================================
# Darstellung
# =============================================================================
def show_effi(param, einheiten=None):
	all_effi = np.array(param.effi_array)
	all_effi[all_effi<0] = 0
	effi_img = np.flipud(np.transpose(all_effi))
	plt.figure()
	if einheiten == None:
		ext = [0, param.I_p, param.tau_p, 0]
		plt.imshow(effi_img*100, aspect='auto', extent=ext)
		plt.colorbar(label="storage efficiency in %")
		plt.ylabel("Pumpduration in s")
		plt.xlabel("Pumpintensity in W/m^2")
	else:
		xkey = "pump density"
		ykey = "pump duration"
		ext = [0, param.woerterbuch[xkey] / einheiten[xkey][1],
		 0, param.woerterbuch[ykey] / einheiten[ykey][1]]
		plt.imshow(all_effi*100, aspect='auto', extent=ext)
		plt.colorbar(label="storage efficiency in %")
		plt.xlabel("Pumpintensity in " + einheiten[xkey][0])
		plt.ylabel("Pumpduration in " + einheiten[ykey][0])
	plt.show()


def show_all_results(param, einheiten=None):
	all_effi = np.array(param.effi_array)
# 	effi_img = np.flipud(np.transpose(all_effi))
	effi_img = np.transpose(all_effi)

	if einheiten == None:
		plt.figure()
		ext = [0, param.I_p, param.tau_p, 0]
		plt.imshow(effi_img*100, aspect='auto', extent=ext)
		plt.clim(0)
		plt.colorbar(label="storage efficiency in %")
		plt.ylabel("Pumpduration in s")
		plt.xlabel("Pumpintensity in W/m^2")
		plt.title(r'efficiency vs. pump duration and $I_p$')

		plt.figure()
		s = all_effi.shape
		ip = np.linspace(0, param.I_p, s[0])
		eoi = all_effi[:,-1]
		plt.plot(ip, eoi*100)
		plt.ylim((0, np.max(eoi)*105))
		plt.ylabel("storage efficiency in %")
		plt.xlabel("Pumpintensity in W/m^2")
		plt.grid()
		plt.title(r'efficiency vs. pump intensity at $\tau_p$ =' + str(param.tau_p*1000) + 'ms')

		plt.figure()
		s = all_effi.shape
		t = np.linspace(0, param.I_p, s[1])
		eot = all_effi[-1,:]
		plt.plot(t, eot*100)
		plt.ylim((0, np.max(eot)*105))
		plt.ylabel("storage efficiency in %")
		plt.xlabel("Pumpduration in s")
		plt.grid()
		plt.title(r'efficiency vs. pump duration at $I_p$ =' + str(param.I_p*1E-7) + 'kW/cm^2')
	else:
		plt.figure()
		xkey = "pump density"
		ykey = "pump duration"
		ext = [0, param.woerterbuch[xkey] / einheiten[xkey][1],
		 param.woerterbuch[ykey] / einheiten[ykey][1], 0]
		plt.imshow(effi_img*100, aspect='auto', extent=ext)
		plt.clim(0)
		plt.colorbar(label="storage efficiency in %")
		plt.xlabel("Pumpintensity in " + einheiten[xkey][0])
		plt.ylabel("Pumpduration in " + einheiten[ykey][0])
		plt.title(r'efficiency vs. pump duration and $I_p$')

		plt.figure()
		s = all_effi.shape
		ip = np.linspace(0, param.woerterbuch[xkey] / einheiten[xkey][1], s[0])
		eoi = all_effi[:,-1]
		plt.plot(ip, eoi*100)
		plt.ylim((0, np.max(eoi)*105))
		plt.ylabel("storage efficiency in %")
		plt.xlabel("Pumpintensity in " + einheiten[xkey][0])
		plt.grid()
		plt.title(r'efficiency vs. pump intensity at $\tau_p$ =' +
			str(param.tau_p / einheiten[ykey][1]) + einheiten[ykey][0])

		plt.figure()
		s = all_effi.shape
		ip = np.linspace(0, param.woerterbuch[ykey] / einheiten[ykey][1], s[1])
		eot = all_effi[-1,:]
		plt.plot(ip, eot*100)
		plt.ylim((0, np.max(eot)*105))
		plt.ylabel("storage efficiency in %")
		plt.xlabel("Pumpduration in " + einheiten[ykey][0])
		plt.grid()
		plt.title(r'efficiency vs. pump duration at $I_p$ =' +
			str(param.I_p/einheiten[xkey][1]) + einheiten[xkey][0])

	plt.show()

# =============================================================================
# Tests
# =============================================================================
if __name__ == "__main__":
	plt.close("all")
	param = param_struct()
	all_effi = efficency(param)
	#show_effi(param)
	show_all_results(param)