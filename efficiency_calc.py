# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 17:06:11 2021

@author: mens
"""

import numpy as np
from inversion_calc import inversion
import matplotlib.pyplot as plt
from init_param import param_struct

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


def show_effi(param, einheiten=None):
	all_effi = np.array(param.effi_array)
	plt.figure()
	all_effi[all_effi<0] = 0
	all_effi = all_effi.transpose()
	if einheiten == None:
		ext = [0, param.I_p, 0, param.tau_p]
		plt.imshow(all_effi*100, aspect='auto', extent=ext)
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


if __name__ == "__main__":
	param = param_struct()
	all_effi = efficency(param)
	show_effi(param)