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

# Load all calculations
# Find which systems first
system_pre_filelist=os.listdir(path=base_directory)
subdir_regex = re.compile('([A-Za-z]+)_photon_NLO_([A-Za-z]+)([0-9\.]+)_ncteq15np(.*)')

calc_dict={}
for tmp_file in system_pre_filelist:

    match=subdir_regex.match(tmp_file)
    if (match != None):
        nuclei_pair=match.group(1)
        system=match.group(2)
        sqrts_GeV=match.group(3)
        extra_label=match.group(4)
        calc_dict[nuclei_pair+"_"+system+sqrts_GeV+extra_label]={'dirname':tmp_file,'sqrts_GeV':float(sqrts_GeV)}

# Loop over systems
for system, system_tmp_dict in calc_dict.items():

    dirname=system_tmp_dict['dirname']

    system_pre_filelist=os.listdir(path=os.path.join(base_directory,dirname))
    subdir_regex = re.compile('incnloLHAPDFmodff_results_summary_var28=0.5_var30=([0-9\.]+)_var31=([0-9\.]+)_var32=([0-9\.]+)')

    calc_dict[system]['calcs']={}

    scale_list=[]

    for tmp_file in system_pre_filelist:

        match=subdir_regex.match(tmp_file)

        if (match != None):
            scale_val = match.group(1)
            tmp_data = np.genfromtxt(os.path.join(base_directory, dirname, tmp_file))
            calc_dict[system]['calcs'][scale_val] = np.column_stack((
                                tmp_data[:, 0],  # pT
                                tmp_data[:, 1],  # dNdpT
                                tmp_data[:, 3] + tmp_data[:, 6]  # dNdpT_LO
                                ))
            scale_list.append(scale_val)

    calc_dict[system]['scale_list']=scale_list



##################################################
##################################################
### Physics basis behind low p_T extrapolation ###
##################################################
##################################################
#
# 
# We look at inclusive high-pT single-particle cross-section in collinear-based pQCD
# dsigma(p_T, norm) = f(x_a, Q_fac) X f(x_b, Q_fac) X dsigma_parton(x_a,x_b,x_c,x_d, Q_ren) X D(z_c, Q_frag)
# where "norm" is defined by the way we set the factorization/renormalization/fragmentation scales equal:
# Q_fac=Q_ren=Q_frag=norm*p_T,
# We observe that the ratio of dsigma(p_T,norm1)/dsigma(p_T,norm2) is approximately a constant
# That is, changing the scales by a constant results mostly in a normalization change of dsigma(p_T,norm)
#
# See prompt_photons_at_low_pt_PhD_Paquet.pdf for more physics background
#
# Note: evidently this extrapolation is a coarse approximation that cannot replace the inclusion
#       of important physics, such as the change of number of active quarks at lower p_T,
#       high-order corrections or non-perturbative effects,
#       as well as parton interactions with the plasma in heavy ion collisions


##############################
### Find the normalization ###
##############################

# One challenge is that the normalization isn't quite constant
# One needs to make a choice as to exactly what normalization to use
# This part of the code looks at the normalization averaged over different ranges of p_T,
# and outputs some information on how much the normalization depends on the p_T cuts
#
# There's no strict way to chose the p_T range of the fit (that is, to choose the normalization)
# but the uncertainty of picking different p_T range seems acceptable

chosen_low_pT_cut_for_fit=4
chosen_high_pT_cut_for_fit=6


# In GeV
low_pT_cut_for_fit_list=[4,5]
high_pT_cut_for_fit_list=[6,7]

norm_dict={}

# Loop over systems
for system, system_tmp_dict in calc_dict.items():

    print("Preparing ", system)

    sqrts_GeV=system_tmp_dict['sqrts_GeV']
    print(sqrts_GeV)

    norm_dict[system]={}

    tmp_calc_dict=calc_dict[system]['calcs']

    #
    for low_pT_cut_for_fit in low_pT_cut_for_fit_list:

        norm_dict[system][low_pT_cut_for_fit]={}

        for high_pT_cut_for_fit in high_pT_cut_for_fit_list:
        
            norm_dict[system][low_pT_cut_for_fit][high_pT_cut_for_fit]={}

            # Normalization calculated in comparison with this calculation
            pT_ref, dNdpt_ref, dNdpt_LO_ref = np.transpose(tmp_calc_dict['0.5'])

            # Loop over scales
            for n, (scale, calc_at_scale) in enumerate(tmp_calc_dict.items()):

#                print("Doing scale",scale)

                pT, dNdpt, dNdpt_LO = np.transpose(calc_at_scale)

                ratio=dNdpt/dNdpt_ref
                spl = UnivariateSpline(pT, ratio,k=1)
                tmp_res=1./(spl.integral(low_pT_cut_for_fit, high_pT_cut_for_fit)/(high_pT_cut_for_fit-low_pT_cut_for_fit))
                norm_dict[system][low_pT_cut_for_fit][high_pT_cut_for_fit][scale]=tmp_res


plot_this=True

if (plot_this):

    # Check how much the normalization depends on the range of p_T
    for system, system_tmp_dict in calc_dict.items():

        # Plot
        plt.figure()
        plt.title(system)
        plt.xlabel("Norm of calculation at given scale over calculation at Qs=p_T/2 for different p_T cuts")

        for scale in calc_dict[system]['scale_list']:

            histogram_list=[]

            # Get all the scales for a given
            for low_pT_cut_for_fit in low_pT_cut_for_fit_list:

                for high_pT_cut_for_fit in high_pT_cut_for_fit_list:

                    histogram_list.append(norm_dict[system][low_pT_cut_for_fit][high_pT_cut_for_fit][scale])

                    #print(low_pT_cut_for_fit,high_pT_cut_for_fit,norm_dict[system][low_pT_cut_for_fit][high_pT_cut_for_fit]['8.0'])

            plt.hist(histogram_list,bins=10,label=str(scale))
            plt.arrow(norm_dict[system][chosen_low_pT_cut_for_fit][chosen_high_pT_cut_for_fit][scale], -0.1, 0, 1)

        plt.legend(loc='upper right',fontsize=10)
        plt.tight_layout()
        plt.show()



############
### Plot ###
############

plot_this=False

if (plot_this):

    font = {'family' : 'URW Gothic',
            'weight' : 'bold',
            'size'   : 16}

    plt.rc('font', **font)

    # From this point on, we assume that the use the normalization calculated with "chosen_low_pT_cut_for_fit" and "chosen_high_pT_cut_for_fit"

    for plot_ratio in [True,False]:

        normalize_iter=[True,False]
        for normalize_calc in normalize_iter:

            for system, system_tmp_dict in calc_dict.items():

                if (not system == "AuAu_RHIC200"):
                    continue

                # Get the info I need
                calcs_dict=system_tmp_dict['calcs']

                # If a ratio is taken, this is the denominator
                pT_ref, dNdpt_ref, dNdpt_LO_ref = np.transpose(calcs_dict['0.5'])
                spl_ref_pre = UnivariateSpline(pT_ref, np.log(dNdpt_ref), k=1, ext='raise')
                spl_ref = lambda pT : np.exp(spl_ref_pre(pT))

                # Plot
                plt.figure()
                plt.title(system)
                plt.xscale('log')
                x_low=0.5
                x_high=100
                plt.xlim(x_low,x_high)
                plt.xlabel(r'$p_T^\gamma$')

                if (plot_ratio):
                    plt.yscale('linear')
                    plt.ylim(0,2)
                    plt.ylabel(r'$dN_{\gamma}/dy dp_T^\gamma ratio to Q=p_T/2$')
                else:
                    plt.ylabel(r'$dN_{\gamma}/dy dp_T^\gamma$')
                    plt.yscale('log')
                    y_low=spl_ref(np.min([x_high,np.max(pT_ref)]))/2
                    y_high=spl_ref(x_low)*2 if x_low > 2 else spl_ref(4)*1e5
                    plt.ylim(y_low,y_high)


                colour_list=['blue','green','black','red','purple','orange','cyan']
#    #symbol_list=['-','--',':','-.','-','--',':']
#    symbol_list=['D','^','x','o','.','d','D']


                for n, (scale, calc_at_scale) in enumerate(calcs_dict.items()):

                    pT, dNdpt, dNdpt_LO = np.transpose(calc_at_scale)

                    if (plot_ratio):
                        y=dNdpt/dNdpt_ref
                    else:
                        y=dNdpt

                    if (normalize_calc):
                        y*=norm_dict[system][chosen_low_pT_cut_for_fit][chosen_high_pT_cut_for_fit][scale]

                    plt.plot(pT, y, color=colour_list[n],label="Q="+scale)



                plt.legend(loc='upper right',fontsize=10)
                plt.tight_layout()

                ratio_label=""
                if (plot_ratio):
                    ratio_label="ratio_"

                normalize_label=""
                if (normalize_calc):
                    normalize_label="normalized_"

                output_filename="dNdpT_prompt_photons_"+ratio_label+normalize_label+system+"_low_range_fit.pdf"

                plt.savefig(output_filename)
                plt.show()


############################################
### Output normalized pQCD extrapolation ###
############################################

scale_to_output='4.0'

for system, system_tmp_dict in calc_dict.items():

    sqrts_GeV=system_tmp_dict['sqrts_GeV']
    # Get the info I need
    calcs_dict=system_tmp_dict['calcs']

    #
    pT, dNdpt, dNdpt_LO = np.transpose(calcs_dict[scale_to_output])
    pT_Qs05, dNdpt_Qs05, dNdpt_LO_Qs05 = np.transpose(calcs_dict['0.5'])
    pT_Qs2, dNdpt_Qs2, dNdpt_LO_Qs2 = np.transpose(calcs_dict['2.0'])

    dNdpt*=norm_dict[system][chosen_low_pT_cut_for_fit][chosen_high_pT_cut_for_fit][scale_to_output]

    ###########################################
    ##### Try to estimate the uncertainty #####
    ###########################################
    # First estimate the scale uncertainty
    # Calculate relative scale uncertainty array
    relative_scale_uncertainty = (dNdpt_Qs05 / dNdpt_Qs2)
    # If p_T's are too low, that scale uncertainty isn't reliable, so I replace it by an average value taken from higher p_T's
    mask_replace = (pT_Qs05 >= 4) & (pT_Qs05 < 10.)
    replacement_value = np.mean(relative_scale_uncertainty[mask_replace])
    relative_scale_uncertainty[pT_Qs05 < 4] = replacement_value

    # The (NLO+LO)/LO ratio doesn't seem to be correlated with the pT regions where I wouldn't trust the calculations.
    # Not used at the moment
    # ratio_NLO_LO = dNdpt / dNdpt_LO
    # Normalize the ratio by its minimum value
    # normalized_ratio_NLO_LO = ratio_NLO_LO / np.min(ratio_NLO_LO)

    # Estimated relative scale uncertainty
    #estimated_rel_uncert=np.multiply(relative_scale_uncertainty,normalized_ratio_NLO_LO)
    estimated_rel_uncert=relative_scale_uncertainty
    #print(relative_scale_uncertainty)
    #print(dNdpt_Qs05)
    #print(dNdpt_Qs2)
    
    # We expect that the current calculations have corrections of order of $\Lambda_QCD/Q$ (to some power) 
    # We add an additional roughly estimated uncertainty based on that
    lambda_qcd_in_GeV=.2
    lambda_over_pT=lambda_qcd_in_GeV/pT
    # Correction factor determined by assuming that the uncertainty should be doubled for pT=1 GeV (LambdaQCD/pT=0.2).
    # Evidently, this is a rough estimate
    estimated_lambda_over_pT_coefficient=5. 
    estimated_rel_uncert*=(1+estimated_lambda_over_pT_coefficient*lambda_over_pT)
    
    dNdpt_estimated_uncert_lower=dNdpt/estimated_rel_uncert
    dNdpt_estimated_uncert_higher=dNdpt*estimated_rel_uncert

    # Save to file
    np.savetxt("pQCD_low_pT_extrapol_"+system+"_low_range_fit.txt",np.transpose([pT,dNdpt,dNdpt_estimated_uncert_lower,dNdpt_estimated_uncert_higher]))


    # Compare with raw calculations at Qs=p_T/2 as final validation

    # Also compare with old fits I made by hand
    hand_fit_dict={
        "AuAu_RHIC39":"fit_by_hand/2018_04_fits_bes/AuAu_prompt_photons_39GeV.dat",
        "AuAu_RHIC62.4":"fit_by_hand/2018_04_fits_bes/AuAu_prompt_photons_62.4GeV.dat",
        "AuAu_RHIC200":"fit_by_hand/AuAu200_ncteq15np_prompt_Qs40_times_1.84.dat",
        "PbPb_LHC2760":"fit_by_hand/PbPb_2760GeV_ncteq15np_prompt_Qs80_times_1.23.dat",
        "PbPb_LHC5020":"fit_by_hand/PbPb_5020GeV_ncteq15np_prompt_Qs80_times_1.17.dat",
        "pPb_LHC5020":"fit_by_hand/pPb_LHC_5020GeV_from_pPb_paper.dat",
        "pp_LHC5020":"fit_by_hand/pp_sqrts5020_from_old_fit_from_pPb_paper.dat",
        "pp_RHIC200":"fit_by_hand/pp_sqrts200_from_old_fit_from_pPb_paper.dat"
    }

    font = {'family' : 'URW Gothic',
            'weight' : 'bold',
            'size'   : 16}

    plt.rc('font', **font)

    # Plot
    plt.figure()
    plt.title(system)
    plt.xscale('log')

    # Set the scale
    def get_max_pT_scale(sqrts):
        if abs(sqrts - 200) < 1e-9:
            return 20
        elif sqrts > 1000:
            return 100
        elif sqrts < 10:
            return 5
        else:
            return np.max([0.2/2.*sqrts,5])
    
    x_low=0.5
    x_high=get_max_pT_scale(float(sqrts_GeV))
    plt.xlim(x_low,x_high)

    # Adjust the y scale accordingly
    mask_select = (pT < x_high)
    y_min = np.min(dNdpt[mask_select])
    y_max = np.max(dNdpt[mask_select])
    y_min_padded = y_min * 0.9
    y_max_padded = y_max * 1.1
    plt.yscale('log')
    plt.ylim(y_min_padded, y_max_padded)

    plt.xlabel(r'$p_T^\gamma$')
    plt.ylabel(r'$dN_{\gamma}/dy dp_T^\gamma$')

    # Plot original Qs=pT/2 calculations
    plt.plot(pT_Qs05, dNdpt_Qs05, 'rD', label="Qs=p_T/2")


    # Plot extrapolated
    plt.plot(pT, dNdpt, '-', color='blue', label="Extrapolated")
    plt.fill_between(pT,
                 dNdpt_estimated_uncert_lower,
                 dNdpt_estimated_uncert_higher,
                 color='blue', alpha=0.3)

    # Plot old hand-fit when available
    if system in hand_fit_dict.keys():
        pT_hand_fit, dN_hand_fit=np.loadtxt(hand_fit_dict[system]).T
        plt.plot(pT_hand_fit, dN_hand_fit, ':', color='black', label="Older hand-fitted extrapolation")

    # Plot data in the p+p case
    #if (system == "pp_photon_NLO_LHC5020_ncteq15np"):
        

    # Show the region between p_T=0 and p_T=4 GeV, which is apparently not the most trustworthy region of p_T for pQCD calculations
    plt.axvspan(0, 4, alpha=0.5, color='grey')


    plt.legend(loc='upper right',fontsize=10)
    plt.tight_layout()

    plt.savefig("validation_pQCD_low_pT_extrapol_"+system+"_low_range_fit.pdf")
    plt.show()
