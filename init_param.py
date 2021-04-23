import numpy as np
import pickle
import json

class param_struct:
	def __init__(self, nutzereingabe=None):

		# Naturkonst
		self.c = 3E8 # [m/s]
		self.h = 6.626e-34 # [Js]

		if nutzereingabe == None:

			# Pumpwellenlänge
			self.lambda_p = 940E-9 # [m]
			self.nu_p = self.c/self.lambda_p # [1/s]

			self.sigma_a = 0.8E-24 # [m²]
			self.sigma_e = 0.16E-24 # [m²]
			self.beq = self.sigma_a / (self.sigma_a + self.sigma_e)
			# Laserwellenlänge
			self.lambda_l = 1030E-9 # [m]
			self.nu_l = self.c/self.lambda_l # [1/s]
			self.sigma_al = 0.11E-24 # [m²]
			self.sigma_el = 2.3E-24 # [m²]
			self.beta_eql = self.sigma_al/(self.sigma_al+self.sigma_el)
			# Materialparameter
			self.tau_f = 0.95E-3 # [s]
			self.N_dop = 6E26 # [1/m²]

			# Pumpparameter
			self.tau_p = 2E-3 # [s]
			self.z_max = 1*1E-2
			self.I_p = 30E7 # [W/m²]

			# Seedparameter
			self.F_in = 0.01E4 # [J/m²]
			self.passes = 6
			self.tau_seed = 5E-9 # [s] 1/e
			self.seedtype = 'gauss' # available shapes: 'rect', 'gaussn' (n is the order for the super-Gaussian pulse of the form e^((1/x))^2n)
			self.losses = 0.02 # loss coefficient per round-trip


			# numerische Paramter
			self.numres = 200
			self.seed_sample_length = 2 * self.tau_seed # [s]
			self.dz = self.z_max/self.numres
			self.z = np.linspace(0, (1*self.z_max), self.numres) #[m]
			self.dt = self.tau_p/self.numres
			self.t = np.linspace(0, self.tau_p, self.numres) #[s]
			self.seedres = self.numres # [pts]

		else:

			# Pumpwellenlänge
			self.lambda_p = nutzereingabe["pumpwavelength"]
			self.nu_p = self.c/self.lambda_p # [1/s]
			self.sigma_a = nutzereingabe["pump absorbtion cross section"]
			self.sigma_e = nutzereingabe["pump emission cross section"]
			self.beq = self.sigma_a / (self.sigma_a + self.sigma_e)

			# Laserwellenlänge
			self.lambda_l = nutzereingabe["laser wavelength"]
			self.nu_l = self.c/self.lambda_l # [1/s]
			self.sigma_al = nutzereingabe["laser absorbtion cross section"]
			self.sigma_el = nutzereingabe["laser emission cross section"]
			self.beta_eql = self.sigma_al/(self.sigma_al+self.sigma_el)

			# Materialparameter
			self.tau_f = nutzereingabe["upper state lifetime"]
			self.N_dop = nutzereingabe["doping concentration"]

			# Pumpparameter
			self.tau_p = nutzereingabe["pump duration"]
			self.z_max = nutzereingabe["material thickness"]
			self.I_p = nutzereingabe["pump density"]

			# Seedparameter
			self.F_in = nutzereingabe["seed fluence"]
			self.passes = int(nutzereingabe["number of material passes"])
			self.tau_seed = nutzereingabe["seed duration"]
			self.seedtype = 'gauss' # available shapes: 'rect', 'gaussn' (n is the order for the super-Gaussian pulse of the form e^((1/x))^2n)
			self.losses = nutzereingabe["losses per round trip"]

			# numerische Paramter
			self.numres = int(nutzereingabe["numerical resolution"])
			self.seed_sample_length = 2 * self.tau_seed # [s]
			self.dz = self.z_max/self.numres
			self.z = np.linspace(0, (1*self.z_max), self.numres) #[m]
			self.dt = self.tau_p/self.numres
			self.t = np.linspace(0, self.tau_p, self.numres) #[s]
			self.seedres = self.numres # [pts]


		# calculate further input data based on definition above
		self.seed_time, self.pulse = self.pulse_gen() # [1/m³]



		# das große Wörterbuch
		self.woerterbuch = {"pumpwavelength": self.lambda_p,
			"pump absorbtion cross section": self.sigma_a,
			"pump emission cross section": self.sigma_e,
			"laser wavelength": self.lambda_l,
			"laser absorbtion cross section": self.sigma_al,
			"laser emission cross section": self.sigma_el,
			"upper state lifetime": self.tau_f,
			"doping concentration": self.N_dop,
			"pump duration": self.tau_p,
			"material thickness": self.z_max,
			"pump density": self.I_p,
			"seed fluence": self.F_in,
			"number of material passes": self.passes,
			"seed duration": self.tau_seed,
			"order of temporal Gaussian pulse shape (2n)": 2,
			"losses per round trip": self.losses,
			"numerical resolution": self.numres
			}


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

		pulse = np.zeros(self.seedres)

		stepsize = self.seed_sample_length / self.seedres

		if self.seedtype == 'rect':
			k1 = int((self.seed_sample_length - self.tau_seed) / (2 * stepsize))
			k2 = int((self.seed_sample_length + self.tau_seed) / (2 * stepsize))
			pulse[k1:k2] = self.F_in / h / self.nu_l / c / self.tau_seed

		elif self.seedtype[0:5] == 'gauss':
			if len(self.seedtype) > 5:
				n = int(self.seedtype[5])
			else:
				n = 1
			t = np.linspace( -(self.seedres-1)/2, (self.seedres-1)/2, self.seedres)*stepsize
			pulse = np.exp( -(t / self.tau_seed * 2) ** (2*n))
			pulse = self.F_in / h / c * self.lambda_l / c / np.sum(pulse) / stepsize * pulse

		return t, pulse


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
	txt = json.dumps(param.woerterbuch)
	txt = txt.replace(",", "\n") # bessere Lesbarkeit
	return txt


def save_txt(param, filename):
	fh = open(filename, "w")
	fh.write(param2txt(param))
	fh.close()

def load_txt(filename):
	fh = open(filename, "r")
	txt = fh.read()
	txt = txt.replace("\n", ",")
	fh.close()
	wb = json.loads(txt)
	param = param_struct(wb)
	return param


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
