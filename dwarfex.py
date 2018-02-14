"""
dwarfex.py

Example use of the abunmatch module and tracks to determine typical stellar
masses of a dwarf satellite with known halo mass and orbit history.
"""
import numpy as np
from astropy import units
import matplotlib.pyplot as plt
import profileclass
import getrtide
import getmlossp10
reload(getmlossp10)
import rubiatracks
import abunmatch as am
import backtracks as bt
reload(bt)
import getafromsidm
reload(getafromsidm)
import masssize

plt.rcParams["figure.figsize"]=[16,12]
plt.rcParams['font.size']=30
plt.rcParams['image.cmap'] = 'plasma'
plt.rcParams['axes.linewidth']=4
plt.rcParams['axes.labelpad']=10
plt.rcParams['xtick.major.size']=9
plt.rcParams['ytick.major.size']=9
plt.rcParams['xtick.major.width']=3
plt.rcParams['ytick.major.width']=3
plt.rcParams['xtick.minor.size']=6
plt.rcParams['ytick.minor.size']=6
plt.rcParams['xtick.minor.width']=1
plt.rcParams['ytick.minor.width']=1
plt.rcParams['ytick.major.pad']=5
plt.rcParams['xtick.major.pad']=10
plt.rcParams['axes.labelsize']=10
plt.rcParams['axes.labelweight']='bold'
plt.rcParams['savefig.bbox']='tight'
plt.rcParams['savefig.dpi']=200
plt.rcParams['lines.linewidth']=4
plt.rcParams['lines.markersize']=10
plt.rcParams['errorbar.capsize']=0
plt.rcParams['hatch.linewidth']=2

####### Parameters for our program are defined and explicitly set here. ######

# An example orbit history of the satellite.
# Pericenter passages, from first to last.
pericenters = np.array([83.3560077756, 60.5349080703, 30.7007751557])*units.kpc
# Pericenter angular velocity in Hz.
periomega = np.array([5.89522379858e-17, 1.14426400057e-16,3.50167625153e-16])/units.s
# Parent halo mass in solar masses at each pericenter.
perimass = np.array([1.42860341124e+12,1.71838324818e+12, 2.16696927838e+12])*units.M_sun
# Parent halo concentration at each pericenter.
pericon = np.array([6.5, 8.8, 11.4])

# The initial satellite halo parameters at infall.
# mass
init_msubhalo = 3.0e9*units.M_sun
# concentration
init_csubhalo = 3.5


# Sample size for sampling the abundance matching distribution of stellar
# masses at fixed subhalo mass.
nsample = 100

# Set a seed (any integer) for our random number generator.
myseed = 38729501

####### CDM subhalo evolution for the given orbit history. ###################

# Make a list of profileclass.NFW objects for the history of the subhalo.
subprofilehist = []
# Add the initial subhalo profile. Note the deltavirrhou should be calculated
# using the subhalo redshift with the Norman and Bryan virial overdensity.
# For the sake of demonstration, I am using a constant virial density.
subprofilehist.append(profileclass.NFW(mvir=init_msubhalo, c=init_csubhalo,
    deltavirrhou=27191*units.M_sun/units.kpc**3))
for i in np.arange(pericenters.size):
    # the parent halo profile during satellite pericenter passage.
    pprofile = profileclass.NFW(mvir=perimass[i], c=pericon[i], 
        deltavirrhou=27191*units.M_sun/units.kpc**3)
    # starting tidal radius of the subhalo
    startrtide = subprofilehist[i].rvir
    # tidal radius at pericenter
    perirtide = getrtide.getrtidenfw(pprofile.rho0, pprofile.rs,
        subprofilehist[i].mvir, pericenters[i], periomega[i])

    # determine subhalo density profile after stripping.
    strippedhalo = getmlossp10.getmlossp10(1, 3, 1, subprofilehist[i].mvir,
        subprofilehist[i].rs, subprofilehist[i].rvir, startrtide, perirtide,
        27191*units.M_sun/units.kpc**3)
    subprofilehist.append(strippedhalo)

####### SIDM subhalo evolution for the given orbit history. ##################

# Carry out the same procedure as with CDM but with SIDM profiles and tracks,
# starting from infall and proceeding to the present.
sidmsubprofilehist = []
# Convert the infall halo profile to a corresponding SIDM profile.


sidmsubprofilehist.append(getafromsidm.getprofilefromsidm( \
        subprofilehist[0].get_vmax(),subprofilehist[0].get_rmax(), \
            13.6*units.Gyr,subprofilehist[0].rvir,27191*units.M_sun/units.kpc**3,0))
# should have correct age, but will fix later

# Strip the halo at each pericenter passage and append the corresponding
# stripped SIDM profiles.
for i in np.arange(pericenters.size):
    # the parent halo profile during satellite pericenter passage.
    pprofile = profileclass.NFW(mvir=perimass[i], c=pericon[i], 
        deltavirrhou=27191*units.M_sun/units.kpc**3)
    # starting tidal radius of the subhalo
    startrtide = subprofilehist[i].rvir
    # tidal radius at pericenter
    perirtide = getrtide.getrtidenfw(pprofile.rho0, pprofile.rs,
        subprofilehist[i].mvir, pericenters[i], periomega[i])

    strippedhalo = getmlossp10.getmlossp10(1, 3, 0, subprofilehist[i].mvir,
        subprofilehist[i].rs, subprofilehist[i].rvir, startrtide, perirtide,
        27191*units.M_sun/units.kpc**3)
    subprofilehist.append(strippedhalo)

####### CDM abundance matching for the subhalo at z=0. #######################

# Initialize a random number generator.
rng = np.random.RandomState(myseed)

# The subhalo mass at z=0.
mhalo = subprofilehist[-2].mvir / units.M_sun
# Array of log(mhalo) to feed into abundance matching.
loghm = np.full(nsample, np.log10(mhalo))

# Generate a sample of CDM stellar masses for the dwarf satellite according
# to our abundance matching model.
logcdmsm = am.am_stellar_logmass_mwsat_sample(loghm, rng=rng)
print logcdmsm

# Use a CDM stellar radius-mass relation for satellites to determine
# half-light radii for the satellite for each stellar mass.
#cdmrstar = rmrelation(logcdmsm) ##### This needs to be implemented.
#cdmrstar = rmrelation(logcdmsm) ##### This needs to be implemented.
cdmrstar=masssize.masssizegammar(10**logcdmsm,0)


####### Generate SIDM stellar masses at z=0. #################################

# An empty array for the SIDM stellar masses.
logsidmsm = np.empty_like(logcdmsm)

for i in np.arange(len(logcdmsm[0:-1])):
    # Tidally disrupt the satellite with the given stellar mass in reverse.
    (mstarinfall, rstarinfall) = \
        bt.cdmgal_backtracks(10**logcdmsm[i]*units.M_sun,
            cdmrstar[i], subprofilehist)

    # Use mstarinfall and rstarinfall with the SIDM sidmsubprofilehist subhalo
    # profile history to determine the stellar mass and radius for SIDM at z=0.

import rubiatracks


print logcdmsm
print logsidmsm
####### Plot the results. ####################################################

plt.clf()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.hist(logcdmsm, histtype='stepfilled', alpha=.2, color='purple', 
         range=[-3,3], bins=50)
plt.hist(logsidmsm, histtype='stepfilled', alpha=.2, color='teal',
         range=[-3,3], bins=50)

plt.hist(logcdmsm, histtype='step', alpha=.2, color='purple', range=[-3,3], 
         bins=50, label='Cusp')
plt.hist(logsidmsm, histtype='step', alpha=.2, color='teal', range=[-3,3], 
         bins=50, label='Core')


plt.xlabel(r'$\log(M_*/M_\odot)$', fontsize=15)
plt.ylabel(r'$N$', fontsize=15)

plt.savefig('stellamass.png', figsize=(8,12))
