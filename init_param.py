import numpy as np
import pickle
import json


class param_struct:
	# große Parameter Struktur, die alle Parameter des Materials, der Pumpe
	# und vom Seed enthält und an die die Simualtionsergebnisse angehängt werden
	# stellt das Interface für so ziemlich alle Funktionen dar
	def __init__(self, nutzereingabe=None):

		if nutzereingabe == None:
			# das große Wörterbuch
			self.woerterbuch = {"pumpwavelength": 940E-9, # [m]
				"pump absorbtion cross section": 0.8E-24, # [m^2]
				"pump emission cross section": 0.16E-24, # [m^2]
				"laser wavelength": 1030E-9, # [m]
				"laser absorbtion cross section": 1.1e-25, # [m^2]
				"laser emission cross section": 2.3E-24, # [m^2]
				"upper state lifetime": 0.95E-3, # [s]
				"doping concentration": 6E26, # [1/m^3]
				"pump duration": 2E-3, # [s]
				"material thickness": 1*1E-2, # [m]
				"pump density": 30E7, # [W/m²]
				"seed fluence": 0.01E4, # [J/m²]
				"number of material passes": 6,
				"seed duration": 5E-9, # [s] 1/e
				"order of temporal Gaussian pulse shape (2n)": 1,
				"losses per round trip": 0.02, # loss coefficient per round-trip
				"numerical resolution": 200
				}
		else:
			self.woerterbuch = nutzereingabe

		# Naturkonst
		self.c = 3E8 # [m/s]
		self.h = 6.626e-34 # [Js]

		# Pumpwellenlänge
		self.lambda_p = self.woerterbuch["pumpwavelength"]
		self.nu_p = self.c/self.lambda_p # [1/s]
		self.sigma_a = self.woerterbuch["pump absorbtion cross section"]
		self.sigma_e = self.woerterbuch["pump emission cross section"]
		self.beq = self.sigma_a / (self.sigma_a + self.sigma_e)

		# Laserwellenlänge
		self.lambda_l = self.woerterbuch["laser wavelength"]
		self.nu_l = self.c/self.lambda_l # [1/s]
		self.sigma_al = self.woerterbuch["laser absorbtion cross section"]
		self.sigma_el = self.woerterbuch["laser emission cross section"]
		self.beta_eql = self.sigma_al/(self.sigma_al+self.sigma_el)

		# Materialparameter
		self.tau_f = self.woerterbuch["upper state lifetime"]
		self.N_dop = self.woerterbuch["doping concentration"]

		# Pumpparameter
		self.tau_p = self.woerterbuch["pump duration"]
		self.z_max = self.woerterbuch["material thickness"]
		self.I_p = self.woerterbuch["pump density"]

		# Seedparameter
		self.F_in = self.woerterbuch["seed fluence"]
		self.passes = int(self.woerterbuch["number of material passes"])
		self.tau_seed = self.woerterbuch["seed duration"]
		self.seedtype = 'gauss' + str(self.woerterbuch["order of temporal Gaussian pulse shape (2n)"]) # available shapes: 'rect', 'gaussn' (n is the order for the super-Gaussian pulse of the form e^((1/x))^2n)
		self.losses = self.woerterbuch["losses per round trip"]

		# numerische Paramter
		self.numres = int(self.woerterbuch["numerical resolution"])
		self.seed_sample_length = 2 * self.tau_seed # [s]
		self.dz = self.z_max/self.numres
		self.z = np.linspace(0, (1*self.z_max), self.numres) #[m]
		self.dt = self.tau_p/self.numres
		self.t = np.linspace(0, self.tau_p, self.numres) #[s]
		self.seedres = self.numres # [pts]

		# calculate further input data based on definition above
		self.seed_time, self.pulse = self.pulse_gen() # [1/m³]

		# Simulationsergebnisse
		self.beta_0 = False
		self.pumprate = False
		self.beta = False
		self.F_ex = False
		self.effi = False
		self.effi_array = False
		self.F_out = False
		self.pulse_out = False


	def pulse_gen(self):
	# function generates a temporal pulse array
		c = 3E8 # [m/s]
		h = 6.626E-34 # [Js]
		
		if len(self.seedtype) > 5 and self.seedtype[0:5] == 'gauss':
			n = int(float(self.seedtype[5::]))
		else:
			n = 1
	
		pulse = np.zeros(self.seedres)

		stepsize = self.seed_sample_length / self.seedres
		t = np.linspace( -(self.seedres-1)/2, (self.seedres-1)/2, self.seedres)*stepsize

		if self.seedtype == 'rect' or n==0:
			k1 = int((self.seed_sample_length - self.tau_seed) / (2 * stepsize))
			k2 = int((self.seed_sample_length + self.tau_seed) / (2 * stepsize))
			pulse[k1:k2] = self.F_in / h / self.nu_l / c / self.tau_seed

		elif self.seedtype[0:5] == 'gauss':
			pulse = np.exp( -(t / self.tau_seed * 2) ** (2*n))
			pulse = self.F_in / h / c * self.lambda_l / c / np.sum(pulse) / stepsize * pulse

		return t, pulse


# =============================================================================
# Parameter speichern und laden
# Textform ist gegenüber Pickle Binärdateien zu bevorzugen
# =============================================================================
def save_pickle(param, filename):
	fh = open(filename, "bw")
	pickle.dump(param.woerterbuch, fh)
	fh.close()

def load_pickle(filename):
	fh = open(filename, "br")
	param = param_struct(pickle.load(fh))
	fh.close()
	return param

def param2txt(param):
	txt = json.dumps(param.woerterbuch, indent=4)
	return txt


def save_txt(param, filename):
	fh = open(filename, "w")
	fh.write(param2txt(param))
	fh.close()

def load_txt(filename):
	fh = open(filename, "r")
	txt = fh.read()
	fh.close()
	wb = json.loads(txt)
	param = param_struct(wb)
	return param

# =============================================================================
# Tests
# =============================================================================
if __name__ == "__main__":
	import matplotlib.pyplot as plt
	param = param_struct()
	save_pickle(param, "default_param.pickle")
	param = load_pickle("default_param.pickle")

	save_txt(param, "default_param.txt")
	param = load_txt("default_param.txt")

	plt.figure()
	plt.plot(param.seed_time, param.pulse)
	plt.ylabel("Photonendichte in 1/(sm^3)")
	plt.xlabel("Zeit in s")
	plt.grid()
	plt.show()
