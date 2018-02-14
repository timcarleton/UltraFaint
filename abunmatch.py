"""
abunmatch.py

Implementations of galaxy abundance matching to dark matter halos or subhalos.
"""
import numpy as np
import numpy.random as rand

def am_stellar_logmass(halo_logmass, logm1=11.5, logms1=9.7, faintslope=1.99,
    brightnorm=3.5, brightslope=0.32):
    """
    Returns the mean of the base 10 logarithm of the stellar mass in solar
    masses of a galaxy in a halo of the given log10 mass in solar masses.
    This is the Behroozi et al. (2013) analytic form.
    
    Arguments:
      halo_logmass: base 10 log of the halo mass (Bryan & Norman defn) in
        solar masses.
    
    Optional Keyword Arguments:
      logm1: base 10 log of the characteristic halo mass in solar masses
        (default 11.5).
      logms1: base 10 log of the mean stellar mass (in solar masses) at the
        characteristic halo mass (default 9.7).
      faintslope: the spectral slope of the faint end of the abundance
        matching relation (default 1.99, the Garrison-Kimmel et al. 2017
        result for zero spread calibrated to the local group, denoted alpha).
      brightnorm: the normalization factor of the bright-end term of the
        abundance matching relation (default 3.5, denoted delta).
      brightslope: the index of the sub-power law of the bright end of the
        abundance matching relation (default 0.32, denoted gamma).
    """
    
    log2 = 0.3010299956639812 # log10(2)
    x = halo_logmass - logm1
    f0 = -log2 + brightnorm * log2**brightslope / (np.e + 1)
    # make overflow errors raise a FloatingPointError.
    np.seterr(over='raise')
    try:
        ffaint = -np.log10(10**(-faintslope*x) + 1)
    except FloatingPointError:
        # for large negative x.
        ffaint = faintslope*x
        
    try:
        fbright = brightnorm / (np.exp(10**(-x)) + 1)
    except FloatingPointError:
        # for moderately large negative x.
        fbright = 0.0
    else:
        try:
            fbright1 = np.log10(np.exp(x) + 1)**brightslope
        except FloatingPointError:
            # for large positive x.
            fbright1 = (x*np.log10(np.e))**brightslope
        finally:
            fbright *= fbright1
            
    return logms1 + ffaint + fbright - f0

def am_stellar_spread_mwsat(halo_logmass, logm1=11.5, brightspread=0.2,
    spreadgrowth=-0.2):
    """
    Returns the abundance matching logarithmic spread of stellar mass for
    satellites of Milky-Way like galaxies. This is the Garrison-Kimmel et al.
    (2017) result.
    
    Arguments:
      halo_logmass: base 10 log of the halo mass (Bryan & Norman defn) in
        solar masses.
      
    Optional Keyword Arguments:
      logm1: base 10 log of the characteristic halo mass in solar masses
        (default 11.5).
      brightspread: the logarithmic spread of the stellar mass for halos
        larger than the characteristic mass (default 0.2).
      spreadgrowth: rate of increase of the spread below the characteristic
        mass (default -0.2, denoted upsilon).
    """
    
    diff = np.array(halo_logmass - logm1)
    # Case where halo_logmass is non-iterable or is an iterable with 1 element.
    if diff.size == 1:
        if diff[0] < 0:
            return brightspread
        else:
            return brightspread + spreadgrowth * diff[0]
    # Case where halo_logmass is an iterable converted to an ndarray.
    myiter = np.nditer([diff, None])
    for di, sp in myiter:
        sp[...] = brightspread if di < 0 else brightspread+spreadgrowth*di
    return myiter.operands[1]

def am_stellar_logmass_mwsat_sample(halo_logmass, logm1=11.5, logms1=9.7,
    brightspread=0.2, spreadgrowth=-0.2, brightnorm=3.5, brightslope=0.32,
    rng=None):
    """
    Returns a random sampling of the abundance matching distribution for a
    subhalo of the given mass that is a satellite of a Milky-Way-like galaxy.
    
    Arguments:
      halo_logmass: base 10 log of the halo mass (Bryan & Norman defn) in
        solar masses.
    
    Optional Keyword Arguments:
      logm1: base 10 log of the characteristic halo mass in solar masses
        (default 11.5).
      logms1: base 10 log of the mean stellar mass (in solar masses) at the
        characteristic halo mass (default 9.7).
      brightspread: the logarithmic spread of the stellar mass for halos
        larger than the characteristic mass (default 0.2).
      spreadgrowth: rate of increase of the spread below the characteristic
        mass (default -0.2, denoted upsilon).
      brightnorm: the normalization factor of the bright-end term of the
        abundance matching relation (default 3.5, denoted delta).
      brightslope: the index of the sub-power law of the bright end of the
        abundance matching relation (default 0.32, denoted gamma).
      rng: a random number generator container, like numpy.random.RandomState
        declared with rng=RandomState(seed) (default None).
    """
    
    # The slope of matching relation for the faint end,
    # for Milky Way satellites. From Garrison-Kimmel et al. (2017).
    fslope = 1.81 - 1.48*spreadgrowth + 0.47*spreadgrowth**2
    
    # The mean log stellar mass.
    logms = am_stellar_logmass(halo_logmass, faintslope=fslope, logm1=logm1,
        logms1=logms1, brightnorm=brightnorm, brightslope=brightslope)
    # The spread of the log stellar mass.
    siglogms = am_stellar_spread_mwsat(halo_logmass, logm1=logm1, 
        brightspread=brightspread, spreadgrowth=spreadgrowth)
    
    if rng == None:
        return rand.normal(logms, siglogms)
    else:
        return rng.normal(logms, siglogms)
