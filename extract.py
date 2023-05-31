import numpy as np
import init_param
import matplotlib.pyplot as plt
from utils import integ
from inversion_calc import inversion




def extracted_pulse(param):
	c = 3E8 # [m/s]
	h = 6.626E-34 # [Js]

	if np.any(param.beta_0):
		beta_0 = param.beta_0
	else:
		beta_0 = inversion(param)

	stepsize = param.seed_sample_length / param.seedres

	pulse_out = np.zeros( (param.passes+1, len(param.pulse)) )
	pulse_out[0,:] = param.pulse
	beta_out = np.zeros( (param.passes+1, len(beta_0)) )
	beta_out[0,:] = np.abs(beta_0)
	F_out = np.zeros(param.passes+1)
	F_out[0] =  np.sum(pulse_out[0,:] ) * stepsize * c * h *c / param.lambda_l
	Delta = param.N_dop * ((beta_out[0,:] - param.beta_eql)/(1 - param.beta_eql))

	for k in range(param.passes):
		# calculate single components
		if k > 1:
			Delta = np.flipud(Delta)

		inver_int = np.exp(-param.sigma_el * np.sum(Delta) * param.dz)
		pulse_int = np.exp(- param.sigma_el * c / (1 - param.beta_eql) * integ(pulse_out[k,:],stepsize))

		# compute pulse shape after amplification
		pulse_out[k+1,:] = pulse_out[k,:] / (1 - ((1 - inver_int) * pulse_int))*(1-param.losses)

		# compute fluence after amplification
		F_out[k+1] =  np.sum(pulse_out[k+1,:] ) * stepsize * c * h * c / param.lambda_l # J/m²

		# refit calculations
		inver_int = np.exp(-param.sigma_el * integ(Delta, param.dz))
		pulse_int = np.exp(param.sigma_el * c / (1 - param.beta_eql) * np.sum(pulse_out[k,:]) * stepsize)

		# compute beta after amplification
		Delta = Delta * inver_int / (pulse_int + inver_int - 1)*(1-param.losses)
		beta_out[k+1,:] = Delta * (1 - param.beta_eql) / param.N_dop + param.beta_eql

	#	 disp(strcat('pass: ', num2str(k), '	beta ratio: ',num2str(sum(erg.beta) / sum(beta_out(k+1,:)))))

		if k % 2:
			beta_out[k+1,:] = np.flipud(beta_out[k+1,:])

	param.F_out = F_out
	param.pulse_out = pulse_out

	max_fluence = np.max(F_out) 
	pump_fluence = param.tau_p * param.I_p
	max_gain = np.max(F_out[1::] / F_out[0:-1])
	print()
	print("Maximal laser fluence in J/cm^2:", max_fluence/1e4)
	print("Total pump fluecne in J/cm^2:", pump_fluence/1e4)
	print("Maximal gain:", max_gain)
	print()
	
	return F_out, pulse_out


def zeige_fluenz(param):
	plt.figure()
	plt.plot(param.F_out*1e-4, "-o")
	plt.ylabel("Fluenz in J/cm^2")
	plt.grid()
	plt.show()

def zeige_pulsform(param):
	plt.figure()
	counter = -1
	t = np.linspace(0, param.seed_sample_length, param.seedres)
	for zeile in param.pulse_out:
		counter += 1
		m = np.max(zeile)
		plt.plot(t*1e9, zeile/m, label="Pass Nr " +str(counter))
	plt.xlabel("Zeit in ns")
	plt.ylabel("Signal in a.u.")
	plt.grid()
	plt.legend()
	plt.show()

def zeige_ergebnisse(param, einheiten=None):
	plt.figure()
	if einheiten == None:
		plt.plot(param.F_out*1e-4, "-o")
		plt.ylabel("Fluence in J/cm^2")
	else:
		ykey = "seed fluence"
		plt.plot(param.F_out / einheiten[ykey][1], "-o")
		plt.ylabel("Fluence in " + einheiten[ykey][0])
	plt.xlabel("pass number")
	plt.grid()

	plt.figure()
	counter = -1
	t = np.linspace(0, param.seed_sample_length, param.seedres)
	for zeile in param.pulse_out:
		counter += 1
		m = np.max(zeile)
		plt.plot(t, zeile/m, label="Pass Nr " +str(counter))
	plt.xlabel("Zeit in ns")
	plt.ylabel("Signal in a.u.")
	plt.grid()
	plt.legend()
	plt.show()


if __name__ == "__main__":

	param = init_param.param_struct()
	F_out, pulse_out = extracted_pulse(param)
	zeige_fluenz(param)
	zeige_pulsform(param)
