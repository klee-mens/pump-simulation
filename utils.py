import numpy as np


# Integrationsroutinen der Form integ(x) = integral(0, x, y dx')
# Konvention integ: 1D
# f(t, z) -> 2D Matrix N x M mit N in Zeit und M in Ort

def integ(y, dx):
    os = np.cumsum(y)
    us = np.array(os)
    us[1::] = os[0:-1]
    us[0] = 0
    return (os + us) / 2 * dx
    

def t_integ(mat, dt):
    os = np.cumsum(mat, axis=0)
    us = np.array(os)
    us[1::, :] = os[0:-1,:]
    us[0,:] = 0
    return (os + us) / 2 * dt

def z_integ(mat, dz):
    os = np.cumsum(mat, axis=1)
    us = np.array(os)
    us[:, 1::] = os[:,0:-1]
    us[:,0] = 0
    return (os + us) / 2 * dz


if __name__ == "__main__":
    from numpy import matlib
    v = np.arange(0, 5)
    print("v: ", v)
    print("integ(v, 1): ", integ(v, 1), " Vergleiche: 0, 0.5, 2.0, 4.5, 8\n")
    m = matlib.repmat(v, 5, 1)
    print("Funktion m in Ort und Zeit")
    print(m)
    print("z_integ(m, 1):")
    print(z_integ(m, 1))
    print("t_integ(m, 1)")
    print(t_integ(m, 1))