import numpy as np
import os
import os.path
import re
import pandas as pd
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

base_directory=os.path.dirname(os.path.realpath(__file__))

font = {'family' : 'URW Gothic',
        'weight' : 'bold',
        'size'   : 16}

plt.rc('font', **font)

##############################
##############################
######## p+p 200 GeV #########
##############################
##############################

# Plot
fig = plt.figure()
#plt.title("p-p 200 GeV")
plt.xscale('linear')
plt.xlim(0,20.)
plt.ylim(1e-11,1e-1)
plt.xlabel(r'$p_T^\gamma$ (GeV)')
plt.ylabel(r'$d\sigma_{\gamma}/dy p_T^\gamma dp_T^\gamma$ (mb GeV$^{-2}$)')
plt.yscale('log')

# Plot phenix data
# p_T, y, dy_stat, dy_syst, p_T-bin-size
#mb GeV$^-2$
filename=os.path.join("data","direct_photons_pp_200GeV_PHENIX2012.dat")
pT, dN, dy_stat, dy_syst = np.loadtxt(filename).T
plt.errorbar(pT,dN,yerr=np.sqrt(dy_stat*dy_stat+dy_syst*dy_syst),fmt="^", color='purple',capsize=4, label="PHENIX") #, barsabove=True)

#pb GeV*-2
pT, dNdpt = np.loadtxt("./pQCD_low_pT_extrapol_pp_RHIC200.txt").T
        
plt.plot(pT, dNdpt*1e-9, color='blue',label="Collinear-based pQCD")

plt.text(0.55, .55, "p+p\n"r"$\sqrt{s_{NN}}=200$ GeV",  ha='left', va='top', transform=plt.axes().transAxes, fontsize=14)

plt.legend(loc='upper right')
plt.tight_layout()

output_filename="vs_data_pp_RHIC200.pdf"
plt.savefig(output_filename)
plt.show()


##############################
##############################
######## p+p 2760 GeV #########
##############################
##############################

# Plot
plt.figure()
plt.title("p-p 2760 GeV")
plt.xscale('linear')
plt.xlim(0,60.)
#plt.ylim(1e-11,1e1)
plt.xlabel(r'$p_T^\gamma$')
plt.ylabel(r'$d\sigma_{\gamma}/dy dp_T^\gamma$')
plt.yscale('log')

# Plot ALICE data
#pT average      val     val-uncert      val+uncert
#pb GEV**-2
filename=os.path.join("data","direct_photons_pp_2760GeV_ALICE2019.dat")
pT, dN, dN_minus_err, dN_plus_err = np.loadtxt(filename).T
plt.errorbar(pT,dN*1e-9,yerr=[1e-9*(dN-dN_minus_err),1e-9*(dN_plus_err-dN)],fmt="^", color='purple',capsize=4, label="ALICE 2019") #, barsabove=True)

##ET(P=3) [GEV],ET(P=3) [GEV] LOW,ET(P=3) [GEV] HIGH,D(N)/DET [PB/GEV],stat +,stat -,sys +,sys -
filename=os.path.join("data","isolated_photons_pp_2760_CMS2013.dat")
rapidity_range=2*1.44
pTmin, pTmax, pT, dN, dN_err_stat_plus, dN_err_stat_minus, dN_err_syst_plus, dN_err_syst_minus = np.loadtxt(filename).T
norm=1e-9/(2.*np.pi*rapidity_range*pT)
plt.errorbar(pT,dN*norm,yerr=[norm*np.sqrt(dN_err_stat_minus**2+dN_err_syst_minus**2),norm*np.sqrt(dN_err_stat_plus**2+dN_err_syst_plus**2)],fmt="^", color='purple',capsize=4, label="ALICE 2019") #, barsabove=True)
#print([norm*np.sqrt(dN_err_stat_minus**2+dN_err_syst_minus**2),norm*np.sqrt(dN_err_stat_plus**2+dN_err_syst_plus**2)])
#print(norm*dN)

#pb GeV*-2
pT, dNdpt = np.loadtxt("./pQCD_low_pT_extrapol_pp_LHC2760.txt").T
        
plt.plot(pT, dNdpt*1e-9, color='blue',label="pQCD")

plt.legend(loc='upper right',fontsize=10)
plt.tight_layout()

output_filename="vs_data_pp_LHC2760.pdf"
plt.savefig(output_filename)
plt.show()
