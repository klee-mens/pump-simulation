# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 15:04:44 2021

@author: mens
"""

import numpy as np
import numpy.matlib as npm
import matplotlib.pyplot as plt #nur für Darstellung
from utils import t_integ, z_integ


# =============================================================================
# Berechnung
# =============================================================================
def inversion(param, testmode=False):
    def beta(R, beta_anf=0, beq=param.beq, tau_f=param.tau_f, dt=param.dt):
        S = np.exp(t_integ( R /beq + 1/tau_f, dt) )
        return 1/S * (t_integ( R*S, dt ) + beta_anf)

    def Rb(beta, N_dop=param.N_dop, sigma=param.sigma_a, beq=param.beq, dz=param.dz):
        return np.exp(z_integ(beta * N_dop*sigma/beq, dz))

    def abbruch(beta_L, beta_H, it):
        max_it = 20
        abweichung = np.abs(np.max(beta_H-beta_L))
        fehler = 1/param.numres
        A = (it > max_it)
        B = (abweichung < fehler)
        if A:
            print("Maximale Iteration von " + str(it) + " überschritten.")
        if B:
            print("Abweichung zwischen beta_L und beta_H =",
               abweichung, "kleiner als", fehler)
        return A or B

    # Initialisierung
    beta_L = np.zeros((param.numres, param.numres))
    beta_H = np.ones((param.numres, param.numres)) * param.beq
    R0 = param.I_p * param.sigma_a / param.h / param.nu_p
    lambert_beer = np.exp(-param.sigma_a * param.N_dop * param.z)
    Ra = R0 * npm.repmat(lambert_beer, param.numres, 1)
    it = 0
    betas = []
    pumprates = []

    # Iteration
    while not abbruch(beta_L, beta_H, it):
        it += 1
        print("Iteration Nummer: " + str(it))
        R_H = Ra * Rb(beta_H)
        R_L = Ra * Rb(beta_L)
        beta_H = beta(R_H)
        beta_L = beta(R_L)

        betas.append((beta_L, beta_H))
        pumprates.append((R_L, R_H))

    print("b_eq:", param.beq)

    if testmode:
        zeige_betas_endtime(betas)
        zeige_pumprates_endtime(pumprates)


    beta_0 = (beta_L[-1,:] + beta_H[-1,:]) / 2
    beta = (beta_L + beta_H) / 2
    pumprate = (R_H + R_L) / 2
    F_ex = param.h * param.nu_l * param.N_dop * np.sum(beta - param.beta_eql, 1) * param.dz
    effi = np.array(F_ex)
    effi[0] = -1
    effi[1::] = np.clip(F_ex[1::] / param.t[1::] / param.I_p, -1, 10.0)
    param.beta = beta
    param.pumprate = pumprate
    param.beta_0 = beta_0
    param.F_ex = F_ex
    param.effi = effi

    return beta_0



# =============================================================================
# Darstellung
# =============================================================================
def zeige_betas_endtime(betas):
    colors = [('tab:blue', '#1f77b4'),
             ('tab:orange', '#ff7f0e'),
             ('tab:green', '#2ca02c'),
             ('tab:red', '#d62728'),
             ('tab:purple', '#9467bd'),
             ('tab:brown', '#8c564b'),
             ('tab:pink', '#e377c2'),
             ('tab:gray', '#7f7f7f'),
             ('tab:olive', '#bcbd22'),
             ('tab:cyan', '#17becf')]
    colors = [p[1] for p in colors]
    plt.figure()
    it = 0
    x = np.linspace(0, param.z_max, param.numres) * 1e3 # in mm
    for betapair in betas:
         line1 = betapair[0][-1,:]
         line2 = betapair[1][-1,:]
         it += 1
         label = "Iteration " + str(it)
         plt.plot(x, line1, label=label, color=colors[it-1])
         plt.plot(x, line2, color=colors[it-1])
    plt.legend()
    plt.title("Inversion")
    plt.xlabel("$z$ in mm")
    plt.ylabel(r"$\beta$")
    plt.grid()

def zeige_pumprates_endtime(pump_rates):
    colors = [('tab:blue', '#1f77b4'),
             ('tab:orange', '#ff7f0e'),
             ('tab:green', '#2ca02c'),
             ('tab:red', '#d62728'),
             ('tab:purple', '#9467bd'),
             ('tab:brown', '#8c564b'),
             ('tab:pink', '#e377c2'),
             ('tab:gray', '#7f7f7f'),
             ('tab:olive', '#bcbd22'),
             ('tab:cyan', '#17becf')]
    colors = [p[1] for p in colors]
    plt.figure()
    it = 0
    x = np.linspace(0, param.z_max, param.numres) * 1e3 # in mm
    for pair in pump_rates:
         line1 = pair[0][-1,:]
         line2 = pair[1][-1,:]
         it += 1
         label = "Iteration " + str(it)
         plt.plot(x, line1, label=label, color=colors[it-1])
         plt.plot(x, line2, color=colors[it-1])
    plt.legend()
    plt.title("Pumprate")
    plt.xlabel("$z$ in mm")
    plt.ylabel(r"$R$ in s$^-1$")
    plt.grid()

def zeige_beta0(param, einheiten=None):
    plt.figure()

    if einheiten == None:
        plt.plot(param.z*1e3, param.beta_0, label="Inversion")
        #plt.plot(param.z*1e3, np.ones(len(param.z))*param.param.beq, color="gray", label=r"$\beta_{eq}$")
        plt.xlabel("z in mm")
    else:
        xkey = "material thickness"
        plt.plot(param.z / einheiten[xkey][1], param.beta_0,
           label="Inversion")
        plt.xlabel("z in " + einheiten[xkey][0])

    plt.ylabel("Inversion")
    plt.grid()
    #plt.legend()
    plt.show()

def show_all_results(param, einheiten=None):
    if einheiten == None:
        plt.figure()
        plt.plot(param.z*1e3, param.beta_0, label="Inversion")
        plt.xlabel("z in mm")
        plt.ylabel("Inversion")
        plt.grid()
        plt.title(r'$\beta$ at the end of the pump duration')

        plt.figure()
        ext = [0, param.z_max, param.tau_p, 0]
        plt.imshow(param.beta, aspect='auto', extent=ext)
        plt.colorbar(label="inversion")
        plt.ylabel("Pumpduration in s")
        plt.xlabel("z in m")
        plt.title(r'$\beta$ vs time and space')

        plt.figure()
        ext = [0, param.z_max, param.tau_p, 0]
        plt.imshow(param.pumprate, aspect='auto', extent=ext)
        plt.colorbar(label="pumprate in s^-1")
        plt.ylabel("Pumpduration in s")
        plt.xlabel("z in m")
        plt.title('pump rate vs time and space')

    else:
        plt.figure()
        xkey = "material thickness"
        plt.plot(param.z / einheiten[xkey][1], param.beta_0,
           label="Inversion")
        plt.xlabel("z in " + einheiten[xkey][0])
        plt.ylabel("Inversion")
        plt.grid()
        plt.title(r'$\beta$ at the end of the pump duration')

        plt.figure()
        xkey = "material thickness"
        ykey = "pump duration"
        ext = [0, param.woerterbuch[xkey] / einheiten[xkey][1],
         param.woerterbuch[ykey] / einheiten[ykey][1], 0]
        plt.imshow(param.beta, aspect='auto', extent=ext)
        plt.colorbar(label="inversion")
        plt.xlabel("z in " + einheiten[xkey][0])
        plt.ylabel("Pumpduration in " + einheiten[ykey][0])
        plt.title(r'$\beta$ vs time and space')

        plt.figure()
        xkey = "material thickness"
        ykey = "pump duration"
        ext = [0, param.woerterbuch[xkey] / einheiten[xkey][1],
         param.woerterbuch[ykey] / einheiten[ykey][1], 0]
        plt.imshow(param.pumprate, aspect='auto', extent=ext)
        plt.colorbar(label="pumprate in s^-1")
        plt.xlabel("z in " + einheiten[xkey][0])
        plt.ylabel("Pumpduration in " + einheiten[ykey][0])
        plt.title('pump rate vs time and space')

    plt.show()



if __name__ == "__main__":
    from init_param import param_struct
    plt.close("all")

    param = param_struct()
    beta0 = inversion(param, testmode=True)
#   zeige_beta0(param)

#   plt.figure()
#   plt.plot(param.t, param.effi)
#   plt.ylabel("Effizienz")
#   plt.xlabel("Pumpzeit")
#   plt.ylim([0, 1.1*np.max(param.effi)])
#   plt.show()

    show_all_results(param)