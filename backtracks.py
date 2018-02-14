"""
backtracks.py

Routines to calculate the tidal stripping of a satellite galaxy in reverse
from present to infall.
"""
import scipy.optimize as opt
import profileclass
import rubiatracks
from astropy import units

def cdmgal_backtracks(mstar, rstar, halo_history):
    """
    Return the stellar mass and half-light radius of a satellite galaxy when
    it accreted onto the parent halo.
    
    Arguments:
      mstar: the stellar mass of the galaxy at present.
      rstar: the half-light radius of the galaxy at present.
      halo_history: list of profileclass objects containing the subhalo
        density profile following each pericenter passage and at infall.
    """
    
    strippedhalo = halo_history[0]
    strippedmstar = mstar
    for haloind in range(1,len(halo_history)):
        starthalo = halo_history[haloind]
        # rovera after stripping.
        strippedrovera = rstar / 1.3 / starthalo.rs.to(units.kpc).value
        # mass loss ratio of this stripping.
        delm = (strippedhalo.get_mass(strippedhalo.get_rmax()) / \
            strippedhalo.get_mass(starthalo.get_rmax())).value
        # solve for rovera before stripping with:
        #   strippedrovera=rubiatracks.getrhfinal(delm, 1, rovera,...)*rovera
        fsresult = opt.fsolve(lambda x: rubiatracks.getrhfinal(delm, 1, x,
            minrmax=True)*x/strippedrovera - 1.0, strippedrovera,full_output=True)

        if fsresult[2] == 1:
          rovera = fsresult[0][0]
        else:
          rovera = strippedrovera
          print(r"Error in cdmgal_backtracks with interpolation using fsolve: %s" % fsresult.mesg)
        # solve for the stellar mass before stripping.
        startmstar = strippedmstar / rubiatracks.getmstarfinal(delm, 1,
            rovera=rovera, minrmax=True)
        
        strippedhalo = starthalo
        strippedmstar = startmstar
    
    return startmstar, rovera * 1.3 * starthalo.rs.to(units.kpc).value
