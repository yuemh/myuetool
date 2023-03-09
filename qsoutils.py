import os, sys
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

import mylib.spectrum.spec_measurement as spec

defaultcosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def Llambda(specobj, redshift, wave, dwave):

    specobj.to_luminosity(redshift)

    mask = (specobj.wavelength>wave-dwave)&(specobj.wavelength<wave+dwave)
    maskflux = specobj.value[mask]

    meanflux = np.mean(maskflux)

    lumi = meanflux * wave * specobj.units[0] * specobj.units[1]

    # erg s-1

    return lumi/(u.erg/u.s)


def Llambda_test(wave, flux, redshift, cwave, dwave):

    cwave = cwave * (1+redshift)
    dwave = dwave * (1+redshift)

    mask = (wave>cwave-dwave)&(wave<cwave+dwave)

    flux_masked = flux[mask]
    avgflux = np.mean(flux_masked)

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    dist = cosmo.luminosity_distance(redshift).to('cm').value

    lumi = avgflux * 4 * np.pi * dist**2 * cwave

    return lumi


def vLv_to_mag(vLv, wave):
    nu = 3e8 * 1e10 / wave

    Lv = vLv / nu
    dist = (10*u.pc).to('cm').value

    fv = Lv / 4 / np.pi / dist**2

    std = 3631e-23

    mag = -2.5*np.log10(fv/std)
    return mag
def MBH_MgII_O09(FWHM_MgII, L3000):
    logBH =  6.86 + \
            2 * np.log10(FWHM_MgII/1000) +\
            0.5 * np.log10(L3000/1e44)

    return logBH


def M1450_to_Lbol(M1450, factor=4.2):
    nu_1450 = 3e8*1e10/1450
    f_1450 = 1e-23 * 3631 * 10**(-0.4*M1450)

    # distance

    dist = 10 * u.pc
    dist_cm = dist.to(u.cm).value

    return nu_1450 * f_1450 * 4 * np.pi * dist_cm**2 * factor
