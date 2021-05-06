# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 11:09:18 2021

@author: mens
"""

import tkinter as tk
import init_param
import inversion_calc
import extract
import efficiency_calc
import tkinter.filedialog as fd
import numpy as np

# =============================================================================
# Eingabefelder und Einheiten festlegen festlegen
# =============================================================================
felder = ["pumpwavelength",
		    "pump absorbtion cross section",
 			"pump emission cross section",
 			"laser wavelength",
 			"laser absorbtion cross section",
 			"laser emission cross section",
 			"upper state lifetime",
 			"doping concentration",
 			"pump duration",
 			"material thickness",
 			"pump density",
 			"seed fluence",
 			"number of material passes",
 			"seed duration",
 			"order of temporal Gaussian pulse shape (2n)",
 			"losses per round trip",
 			"numerical resolution"
 			]

einheiten = [("nm", 1e-9),
			("cm^2", 1e-4),
			("cm^2", 1e-4),
			("nm", 1e-9),
			("cm^2", 1e-4),
			("cm^2", 1e-4),
			("ms", 1e-3),
			("cm^-3", 1e+6),
			("ms", 1e-3),
			("mm", 1e-3),
			("kW/cm^2", 1e7),
			("mJ/cm^2", 10),
			("", 1),
			("s", 1),
			("", 1),
			("%", 1e-2),
			("", 1)]
eh = []
for i in range(len(felder)):
	eh.append([felder[i], einheiten[i]])
einheiten = dict(eh)


# =============================================================================
# Gui-Aufbau und Funktionen der Buttons
# =============================================================================
def eintrage_auf_param_set(param, eintraege):
	wb = param.woerterbuch
	for feld in felder:
		v = wb[feld] / einheiten[feld][1]
# 		p = np.floor(np.log(v) / np.log(10))
# 		v = np.round(v * 10**(-p), 3) * 10**p
		if v >= 1 and v < 1000:
			eintraege[feld].set(int(np.round(v)))
		else:
			eintraege[feld].set(format(v, ".2e"))

def eintraege_erstellen(root, felder, param=None):
	eintraege = {}
	for feld in felder:
		reihe = tk.Frame(root)
		labtext = feld + " [" + einheiten[feld][0] + "]"
		lab = tk.Label(reihe, width=40, text=labtext, anchor="w")
		text_var = tk.StringVar()
		ent = tk.Entry(reihe, textvariable=text_var)
		reihe.pack(side=tk.TOP, padx=5, pady=5)
		lab.pack(side=tk.LEFT)
		ent.pack(side=tk.RIGHT)
		eintraege[feld] = text_var
	return eintraege

def param_aktualisieren(eintraege):
	nutzereintrag = {}
	for feld in felder:
		nutzereintrag[feld] = float(eintraege[feld].get()) * einheiten[feld][1]
	param2 = init_param.param_struct(nutzereintrag)
	return param2

def inverion_berechnen():
	param = param_aktualisieren(eintraege)
	beta_0 = inversion_calc.inversion(param)
	#inversion_calc.zeige_beta0(param, einheiten)
	inversion_calc.show_all_results(param, einheiten)

def effizienz_berechnen():
    param = param_aktualisieren(eintraege)
    all_effi = efficiency_calc.efficency(param)
    efficiency_calc.show_all_results(param, einheiten)

def extraction_berechnen():
	param = param_aktualisieren(eintraege)
	F_out, pulse_out = extract.extracted_pulse(param)
	extract.zeige_ergebnisse(param, einheiten)

def speichern():
	param = param_aktualisieren(eintraege)
	fh = fd.asksaveasfile(initialfile = 'Untitled.txt', defaultextension=".txt",filetypes=[("All Files","*.*"),("Text Documents","*.txt")])
	fh.write(init_param.param2txt(param))
	fh.close()

def laden():
	name = fd.askopenfilename()
	param = init_param.load_txt(name)
	eintrage_auf_param_set(param, eintraege)


# =============================================================================
# GUI starten
# =============================================================================
if __name__ == "__main__":
	param = init_param.param_struct()

	root = tk.Tk()
	eintraege = eintraege_erstellen(root, felder)
	eintrage_auf_param_set(param, eintraege)
	b1 = tk.Button(root, text="Inversion", command=inverion_berechnen)
	b1.pack(side=tk.LEFT, padx=5, pady=5)
	b2 = tk.Button(root, text="Efficinecy", command=effizienz_berechnen)
	b2.pack(side=tk.LEFT, padx=5, pady=5)
	b3 = tk.Button(root, text="Amplification", command=extraction_berechnen)
	b3.pack(side=tk.LEFT, padx=5, pady=5)
	b4 = tk.Button(root, text="Save", command=speichern)
	b4.pack(side=tk.RIGHT, padx=5, pady=5)
	b5 = tk.Button(root, text="Load", command=laden)
	b5.pack(side=tk.RIGHT)
	root.mainloop()