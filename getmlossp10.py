import numpy as np
import profileclass
reload(profileclass)
import rubiatracks

#get the density profile of a halo experiencing tidal mass loss, according to Penarrubia 2010

def getmlossp10(alpha0,beta0,gamma0,mvir0,rs0,rvir,rtide0,rtide1,finaldeltavirrhou):

    #these parameters specify the shape of the halo before the tidal interaction:
    #alpha0, beta0, gamma0 specify the inital profile (Penarrubia 2010, eqn 2)
    #mvir0: the initial virial mass of the subhalo
    #rs0: the initial scale radius of the subhalo
    #rvir: the initial virial radius of the subhalo

    #rtide0 and rtide1 are the tidal radii of the subhalo before and during the tidal interaction

    #finaldeltavirrhou: the virial parameter during the tidal interaction 


    #if the tidal radius is outside the virial radius, no tidal stripping occurs
    if rvir<=rtide1:
        return profileclass.Zhao(alpha0,beta0,gamma0,mvir=mvir0,rs=rs0,rvir=rvir)

    #initialize the profile before the tidal interaction
    rho0=profileclass.Zhao(alpha0,beta0,gamma0,mvir=mvir0,rs=rs0,rvir=rvir)
#    rho0=profileclass.Zhao(alpha0,beta0,gamma0,mvir=mvir0,rs=rs0,deltavirrhou=finaldeltavirrhou)

    #the mass within the tidal radius before the tidal interaction
    msa0=rho0.get_mass(rtide0)
    #the mass within the tidal radius of during the tidal interaction
    msnew=rho0.get_mass(rtide1)
    
    msms0=msnew/msa0

    #if >10% mass is lost, change beta to 5
    #if msms0<=.9:
    #    rho1=profileclass.Zhao(alpha0,5,gamma0,rho0=rho0.rho0,rs=rs0,deltavirrhou=finaldeltavirrhou)
    #else:
    #    rho1=rho0

        
    msnew1=rho0.get_mass(rtide1)
    msrt0=rho0.get_mass(rtide0)

    #the ratio between the mass within the tidal radius before the tidal interaction and during the tidal interaction
    msms01=(msnew1/msrt0).value

    #use the tracks to get the new vmax, rmax
    rmax=rubiatracks.getrmaxfinal(msms01,gamma0,minrmax=False)*rho0.get_rmax()
    vmax=rubiatracks.getvmaxfinal(msms01,gamma0,minrmax=False)*rho0.get_vmax()
    print rmax,vmax
    
    #vmax, rmvax and the virial parameter determine the final halo shapes
    if msms01<=.9:
        rhof=profileclass.Zhao(alpha0,5,gamma0,vmax=vmax,rmax=rmax,deltavirrhou=finaldeltavirrhou)
    else:
        rhof=profileclass.Zhao(alpha0,beta0,gamma0,vmax=vmax,rmax=rmax,deltavirrhou=finaldeltavirrhou)
            
    return rhof
