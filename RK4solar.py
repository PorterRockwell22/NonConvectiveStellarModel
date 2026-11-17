import matplotlib.pyplot as plt
from math import exp, sqrt 
from scipy.integrate import quad
from scipy.optimize import brentq
import numpy as np
import os

#constants(SI units)
G = 6.67e-11 #grav constant
m_p = 1.672e-27 #hydrogen mass
M_o = 1.989e30 #mass of sun
R_0 = 6.957e8 #measured radius of sun
Q = 13.35 * 1.60218e-13 #energy released in a single p-p reaction
L_s = 3.828e26 #surface luminosity
a, c = 7.5657e-16, 2.998e8  #radiation constant, speed of light
k_B_keV = 8.617333e-8 #keV per K
m_u = 1.66054e-27 #molecular mass         
k_B = 1.380649e-23 #boltzman constant

#limits for solver
T_surf = 2e4
P_frac = 1e-8
R_max  = 2.0

DATA_DIR = os.path.dirname(os.path.abspath(__file__))

def reactivity(kT, eta=3.52):
    keV, mu = 1.602176634e-16, 1.67262192e-27 / 2
    S = lambda E: 4e-22 * 1e-28 * (1 + 0.0112 * E)
    I = quad(lambda E: S(E) * np.exp(-2*np.pi*eta/np.sqrt(E) - E/kT), 0, np.inf, epsabs=0, epsrel=1e-8, limit=200)[0]
    return 4 / np.sqrt(2*np.pi*mu) * np.sqrt(keV) * kT**-1.5 * I

_Tgrid = np.logspace(6, 7.7, 120)
_logR = np.log10([reactivity(T * k_B_keV) for T in _Tgrid])
def reactivity_T(T):
    if T < _Tgrid[0]: return 0.0
    return 10**np.interp(np.log10(T), np.log10(_Tgrid), _logR)

def runge_kutta(f, x_0, y_0, h):
    k_0 = f(x_0, y_0)
    k_1 = f(x_0 + h/2, y_0 + h/2 * k_0)
    k_2 = f(x_0 + h/2, y_0 + h/2 * k_1)
    k_3 = f(x_0 + h, y_0 + h * k_2)

    k = 1/6 * (k_0 + 2.0*k_1 + 2.0*k_2 + k_3)

    x_1 = x_0 + h
    y_1 = y_0 + h * k

    return x_1, y_1

X, Y, Z = 0.71, 0.271, 0.01386 #composition of star (H, He, metals)
mu_w = 1 / (2*X + 0.75*Y + 0.5*Z)
GBF_T = 0.05

#kappa
def opacity(p, T, gbf_t=GBF_T, gff=1.0):
    kes = 0.02 * (1 + X) #electron scattering
    kff = 3.68e18 * gff * (1 - Z) * (1 + X) * p * T**-3.5 #free free
    kbf = 4.34e21 * gbf_t * Z * (1 + X) * p * T**-3.5 #bound free
    return kes + kff + kbf

def density(P, T):
    return mu_w * m_u * max(P - a*T**4/3, 0.0) / (k_B * T)

def f(r, y):
    M, P, T, L = y
    p = density(P, T)
    n_d = X * p / m_p
    dM = 4*np.pi*r**2 * p 
    dP = -G*M*p / r**2
    eps = Q * (n_d**2 / 2) * reactivity_T(T) / p if p > 0 else 0.0
    dL = 4*np.pi * r**2 * p * eps
    dT = -3 * opacity(p, T) * p * L / (16*np.pi * a * c * r**2 * T**3)
    return np.array([dM, dP, dT, dL])

def integrate(P_c, T_c, h=R_0/4000, keep=False):
    r = h
    p_c = density(P_c, T_c)
    eps_c = Q * (X*p_c/m_p)**2 / 2 * reactivity_T(T_c) / p_c
    #(typically zero but are given very low values to avoid singularity)
    M_0 = 4/3 * np.pi * p_c * r**3
    L_0 = 4/3 * np.pi * p_c * r**3 * eps_c

    y = np.array([M_0, P_c, T_c, L_0])
    out = [(r, *y)]
    while r < R_max * R_0:
        r_new, y_new = runge_kutta(f, r, y, h)
        if not np.all(np.isfinite(y_new)) or y_new[2] < T_surf:
            return "cold", r, y, out #temperture below bounds
        if y_new[1] < P_frac * P_c:
            return "hot", r_new, y_new, out #pressure below bounds
        r, y = r_new, y_new
        if keep:
            out.append((r, *y))
    return "hot", r, y, out #radii above limit

def T_csolve(P_c, lo=3e6, hi=4e7, tol=1e-12):
    while hi/lo - 1 > tol:
        mid = np.sqrt(lo*hi)
        if integrate(P_c, mid)[0] == "cold":
            lo = mid
        else:
            hi = mid
    return hi

def massmismatch(logP_c):
    P_c = np.exp(logP_c)
    T_c = T_csolve(P_c)
    _, _, y, _ = integrate(P_c, T_c)
    return (y[0] - M_o) / M_o

if __name__ == "__main__":
    #shooting for inital T_c and P_c for convectionless system (cant use known values because convection affects them)
    print("Solving for core Pressure...")
    logP_c = brentq(massmismatch, np.log(1.2e16), np.log(2.4e16), xtol=1e-6)
    P_c = np.exp(logP_c)
    print(f"Output: {P_c:.3e}")
    print("Solving for core Temperature...")
    T_c = T_csolve(P_c)
    print(f"Output: {T_c:.4g}")

    #Use solved values for a final output
    status, r, y, out = integrate(P_c, T_c, keep=True)
    r_values, M_values, P_values, T_values, L_values = np.array(out).T
    p_values = np.array([density(P, T) for P, T in zip(P_values, T_values)])
    kappa_values = np.array([opacity(p, T) for p, T in zip(p_values, T_values)])

    print("Final Results:")
    print(f"P_c = {P_c:.3e} Pa, T_c = {T_c:.4g} K, rho_c = {p_values[0]:.3g} kg/m^3")
    print(f"R = {r_values[-1]/R_0:.2f} R_sun, M = {M_values[-1]/M_o:.3f} M_sun")
    print(f"L = {L_values[-1]/L_s:.2f} L_sun, T_surf = {T_values[-1]:.3g} K")
    
    #plotting
    x = r_values / R_0
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    ax[0,0].plot(x, p_values);       ax[0,0].set_ylabel("rho (kg/m^3)")
    ax[0,1].plot(x, M_values/M_o);   ax[0,1].set_ylabel("M / M_sun")
    ax[0,2].plot(x, P_values);       ax[0,2].set_ylabel("P (Pa)"); ax[0,2].set_yscale("log")
    ax[1,0].plot(x, T_values/1e6);   ax[1,0].set_ylabel("T (MK)")
    ax[1,1].plot(x, L_values/L_s);   ax[1,1].set_ylabel("L / L_sun")
    ax[1,2].plot(x, kappa_values);   ax[1,2].set_ylabel("kappa (m^2/kg)"); ax[1,2].set_yscale("log")
    for a_ in ax.flat:
        a_.set_xlabel("r / R_sun")
    plt.tight_layout()
    plt.savefig(os.path.join(DATA_DIR, 'solar-structure.png'), dpi=200)
    plt.show()

    #temperture mapping
    r = np.linspace(0, R_0, len(T_values))
    theta = np.linspace(0, 2*np.pi, 200)

    Theta, Rgrid = np.meshgrid(theta, r, indexing="ij")
    T = np.tile(T_values, (len(theta), 1))

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    pc = ax.pcolormesh(Theta, Rgrid, T, shading="auto", cmap="inferno")
    fig.colorbar(pc, label="Temperature (K)")
    fracs = np.array([0.25, 0.5, 0.7])

    ax.set_title("Temperture Mapping of The Sun")
    ax.set_rticks(fracs * R_0)
    ax.set_yticklabels([f"{f:.2f} R₀" for f in fracs], color="cyan")
    ax.set_rlabel_position(45)
    ax.set_xticklabels([])
    ax.grid(True, axis="y", alpha=0.4, color="white")
    plt.savefig(os.path.join(DATA_DIR, 'temperturegraph.png'), dpi=200)
    plt.show()
    
    #energy released mapping
    reactivity_vals = np.array([reactivity_T(T) for T in T_values])
    eps_val = Q * (X*p_values/m_p)**2 / 2 * reactivity_vals / p_values
    r = np.linspace(0, R_0, len(eps_val))
    theta = np.linspace(0, 2*np.pi, 200)

    Theta, Rgrid = np.meshgrid(theta, r, indexing="ij")
    eps = np.tile(eps_val, (len(theta), 1))

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    pc = ax.pcolormesh(Theta, Rgrid, eps, shading="auto", cmap="inferno")
    fig.colorbar(pc, label="$\epsilon$ (J/kg*s)")
    fracs = np.array([0.25, 0.5, 0.7])
    
    ax.set_title("Energy Release/dMdt")
    ax.set_rticks(fracs * R_0)
    ax.set_yticklabels([f"{f:.2f} R₀" for f in fracs], color="cyan")
    ax.set_rlabel_position(45)
    ax.set_xticklabels([])
    ax.grid(True, axis="y", alpha=0.4, color="white")
    plt.savefig(os.path.join(DATA_DIR, 'fusiongraph.png'), dpi=200)
    plt.show()
    
