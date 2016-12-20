#! /usr/bin/env python
import numpy as np
import astropy.io.fits as pf
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import math
import matplotlib.gridspec as gridspec
import FrenchMcBands
import PSFfit
import os
from scipy.special import erf
from astropy.stats import poisson_conf_interval
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import argparse

"""
For one specific config, fit the PSF for each MC simulation by a tripplegauss 
Then store the PSF tripplegauss parameters in 4D numpy table for each value of the Zenithal angle, Offset, Efficiency and Energy used for the MCs simulation
Example of commande line to run to create this 4D table with the directory of the MC simulation output and the config name as argument


./PSFtable_script.py 'ash_south_stereo' 180
./PSFtable_script.py 'ash_north_stereo' 0
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Store the PSFs from Mc simulation in a 4D numpy table')
    parser.add_argument('config', action="store", help='Config')
    parser.add_argument('az_angle', action="store", help='azimuth Angle')
    results = parser.parse_args()
    print "Store the PSF in a 4D table from the MC simulations in ", os.path.expandvars('$HESSCONFIG') , " and for the config ", results.config, "matching with a zenith angle", results.az_angle

    """
    Fonction defenition
    """

    def triplegauss(theta2,s1,s2,s3,A2,A3):
        
        s12 = s1*s1
        s22 = s2*s2
        s32 = s3*s3

        gaus1 = np.exp(-theta2/(2*s12))
        gaus2 = np.exp(-theta2/(2*s22))
        gaus3 = np.exp(-theta2/(2*s32))

        y = (gaus1 + A2*gaus2 + A3*gaus3) 
        norm =  2*math.pi*(s12+ np.abs(A2) * s22 + np.abs(A3) * s32)
        return y/norm

    def Integral_triplegauss(theta2min,theta2max,s1,s2,s3,A2,A3):
        s12 = s1*s1
        s22 = s2*s2
        s32 = s3*s3

        gaus1 = np.exp(-theta2min/(2*s12))-np.exp(-theta2max/(2*s12))
        gaus2 = np.exp(-theta2min/(2*s22))-np.exp(-theta2max/(2*s22))
        gaus3 = np.exp(-theta2min/(2*s32))-np.exp(-theta2max/(2*s32))

        y = 2*math.pi*(s12*gaus1 + A2*s22*gaus2 + A3*s32*gaus3) 
        norm =  2*math.pi*(s12+ np.abs(A2) * s22 + np.abs(A3) * s32)
        return y/norm

    def theta2_bin(data, Nev_bin, Nbinmax):
        """
        Define an adaptative theta2binning in order to have at least Nev_bin events per bin and a maximal number of bin of Nbinmax
        You give the theta2 data, the minimum number of events per bin (Nev_bin) and the maximal number of bin you want (Nbinmax).
        """
        Nev=len(data)
        while(Nev/Nev_bin >  Nbinmax):
            Nev_bin = Nev_bin * 2
        theta2sorted=np.sort(data)
        theta2bin=np.array(theta2sorted[0])
        index=np.arange(Nev_bin, Nev, Nev_bin)
        for i in index:
            theta2bin = np.append(theta2bin, theta2sorted[i]) 
        return theta2bin


    def R68(s1,s2,s3,A2,A3, th2max=0.3):
        x=np.linspace(0,th2max,3000)
        y=triplegauss(x,s1,s2,s3,A2,A3)
        res= y.cumsum()/y.sum()
        ind_R68=np.where(res>=0.68)[0][0]
        s68=np.sqrt(0.5*(x[ind_R68]+x[ind_R68-1]))
        return s68

    def R68_hist(theta2):
        Nvalue=len(theta2)
        theta2_sorted=np.sort(theta2)
        ind_R68=int(0.68*Nvalue)
        s68=np.sqrt(0.5*(theta2_sorted[ind_R68]+theta2_sorted[ind_R68-1]))
        return s68

    def king(theta2,sig, gam):
        norm = (1/(2*np.pi*sig**2))*(1-1/gam)
        king=(1+theta2/(2*gam*(sig**2)))**(-gam)
        return norm*king

    def bin_contigu(x,data,model_fun, threshold):
        """
        Regarde combien de valeur sont en-dessous ou au-dessus du fit
        """
        i_sup=np.where(data > model_fun(x))[0]
        i_inf=np.where(data < model_fun(x))[0]
        List_sup=[]
        List_inf=[]
        iband=0
        for i,i_s in enumerate(i_sup.tolist()):
            if(i==0):
                List_sup.append([i_s])
            elif(i_s == i_sup[i-1]+1):
                List_sup[iband].append(i_s)
            else:
                List_sup.append([i_s])
                iband += 1
        iband=0
        for i,i_i in enumerate(i_inf.tolist()):
            if(i==0):
                List_inf.append([i_i])
            elif(i_i == i_inf[i-1]+1):
                List_inf[iband].append(i_i)
            else:
                List_inf.append([i_i])
                iband += 1
        Npoints=len(x)
        Nsup=len(List_sup)
        Ninf=len(List_inf)
        for i in range(Nsup):
            if(len(List_sup[i])>= threshold*Npoints):
                print "WARNING: There are ",len(List_sup[i]) ," values that are superior to the fit"
        for i in range(Ninf):
            if(len(List_inf[i])>= threshold*Npoints):
                print "WARNING: There are ",len(List_inf[i]) ," values that are inferiror to the fit"

        return i_sup,List_sup,i_inf,List_inf

    # Figure definitions
    def khi2_int(x,data,err,model_fun):
        resid = (data - model_fun)**2/err**2
        khi2=np.sum(resid/len(x))
        return khi2

    
    
    """
    We keep the events that have a theta2 inferior to 0.3
    """
    theta2max=0.3

    """
    MC energy, zenithal angle, offset and efficiency
    """
    enMC = [0.02, 0.03, 0.05, 0.08, 0.125, 0.2, 0.3, 0.5, 0.8, 1.25, 2, 3, 5, 8, 12.5, 20, 30, 50, 80, 125]                        
    zenMC = [0, 18, 26, 32, 37, 41, 46, 50, 53, 57, 60, 63, 67, 70]    
    effMC = [50, 60, 70, 80, 90, 100]    
    offMC = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]    
    binEMC = len(enMC)
    binzen = len(zenMC)
    binoff = len(offMC)
    bineff = len(effMC)


    # Size of the 4D table where we will stock parameters of the tripple gaussian fit on the MC data for each true energy, zenithal angle, offset and efficiency of the MCs
    TableSigma1 = np.zeros((binEMC, binoff, binzen, bineff))
    TableSigma2 = np.zeros((binEMC, binoff, binzen, bineff))
    TableSigma3 = np.zeros((binEMC, binoff, binzen, bineff))
    TableA2 = np.zeros((binEMC, binoff, binzen, bineff))
    TableA3 = np.zeros((binEMC, binoff, binzen, bineff))


    MCband=FrenchMcBands.FrenchMcBands()
    directory=os.path.expandvars('$HESSCONFIG')
    config=results.config
    directory=directory+"/"+config
    file_nosimu = open(directory+"/file_nosimu.txt", "w")
    file_toofewevents = open(directory+"/file_toofewevents.txt", "w")
    file_khi2toohigh = open(directory+"/file_khi2sup2.txt", "w")
    file_nosimu.write("Runnumber \t E(Tev) \t Zen \t theta \t Eff \n")
    file_toofewevents.write("Runnumber \t E(Tev) \t Zen \t theta \t Eff \t Number of events \n")
    file_khi2toohigh.write("Runnumber \t E(Tev) \t Zen \t theta \t Eff \t Number of events \t khi2 \n")
    for (ieff, eff) in enumerate(effMC):
            for (ioff, off) in enumerate(offMC):
                for (izen, zen) in enumerate(zenMC):
                    #Initial parameter fot the fit
                    s1_init=0.02
                    s2_init=0.05
                    s3_init=0.08
                    A2_init=0.3
                    A3_init=0.1
                    for (ien, E) in enumerate(enMC):
                        #Calculate the runnnumber for the MC zenithal angle, offset and energy
                        run_number=MCband.run_number(int(results.az_angle),zen, off, E)
                        PSFfile=directory+"/"+str(eff)+"/run_"+run_number+"_Eff"+str(eff)+"_psf.fits"
                        try: 
                            hdu=pf.open(PSFfile)

                        except:
                            print("Cannot open file: " + PSFfile)
                            print("skipping run")
                            file_nosimu.write(run_number+"\t"+str(E)+"\t"+str(zen)+"\t"+str(off)+"\t"+str(eff)+"\n")
                            #Default value to -1000 if the fit the MCs simulation doesn t exist
                            TableSigma1[ien, ioff, izen, ieff] = -1
                            TableSigma2[ien, ioff, izen, ieff] = -1
                            TableSigma3[ien, ioff, izen, ieff] = -1
                            TableA2[ien, ioff, izen, ieff] = -1
                            TableA3[ien, ioff, izen, ieff] = -1
                            continue

                        #Select the events with a theta2 inferior to thetamax
                        theta2 = hdu[1].data["MC_ThSq"]
                        index = [theta2<theta2max]
                        theta2f = theta2[index]
                        #If there are less than 40 events, the fit is not done and we put a default value to -1000
                        if(len(theta2f)<40):
                            file_toofewevents.write(run_number+"\t"+str(E)+"\t"+str(zen)+"\t"+str(off)+"\t"+str(eff)+"\t"+str(len(theta2f))+"\n")
                            TableSigma1[ien, ioff, izen, ieff] = -1
                            TableSigma2[ien, ioff, izen, ieff] = -1
                            TableSigma3[ien, ioff, izen, ieff] = -1
                            TableA2[ien, ioff, izen, ieff] = -1
                            TableA3[ien, ioff, izen, ieff] = -1
                            continue
                        #We define for the theta2binning a minimum of 10 events per bin and a maximum of 50 bins
                        Nev_perbin=10
                        Nbinmax=50
                        theta2hist=theta2_bin(theta2f, Nev_perbin, Nbinmax)
                        hist, bin_edges = np.histogram(theta2,theta2hist)
                        #Me renvois la valeur moyenne en theta2 des evenements stockes dans les bins donc peut etre un peu mieux que de prendre thetabi=(Emin+Emax)/2
                        #theta2bintest, bin_edgestest,a = stats.binned_statistic(theta2,theta2,'mean',theta2hist)
                        #histtest, bin_edgestest,b = stats.binned_statistic(theta2,theta2,'count',theta2hist)
                        PSF=PSFfit.PSFfit(theta2f)
                        s1,s2,s3,A2,A3=PSF.minimization("triplegauss",s1_init, s2_init, s3_init, A2_init, A3_init)
                        #The initial parameters for the fit are the one fit on the previous MC energy
                        s1_init=s1
                        s2_init=s2
                        s3_init=s3
                        A2_init=A2
                        A3_init=A3
                        TableSigma1[ien, ioff, izen, ieff] = s1
                        TableSigma2[ien, ioff, izen, ieff] = s2
                        TableSigma3[ien, ioff, izen, ieff] = s3
                        TableA2[ien, ioff, izen, ieff] = A2
                        TableA3[ien, ioff, izen, ieff] = A3
                        #If the energy bin are in log, we have to take sqrt(Emin*Emax) for the center of the bin
                        #theta2bin = np.sqrt(bin_edges[:-1] * bin_edges[1:])
                        theta2bin = (bin_edges[:-1] + bin_edges[1:])/2.

                        #We have to divide by the solid angle of each bin= pi*dO^2 to normalize the histogram
                        bsize = np.diff(bin_edges)*math.pi
                        hist_norm = hist/float(np.sum(hist))/bsize
                        #bsizetest = np.diff(bin_edgestest)*math.pi
                        #hist_normtest = histtest/float(np.sum(histtest))/bsizetest
                        # use gehrels errors for low counts (http://cxc.harvard.edu/sherpa4.4/statistics/)
                        hist_err = (1+np.sqrt(hist+0.75))/float(np.sum(hist))/bsize
                        #hist_err2 = (1+np.sqrt(histtest+0.75))/float(np.sum(histtest))/bsizetest
                        #Erreur prenant en compte poisson du coup j ai des erreurs asymetrics inf et sup
                        #histerr_test=poisson_conf_interval(hist)/float(np.sum(hist))/bsize
                        Int_fitgauss = lambda x1,x2 : Integral_triplegauss(x1,x2,s1,s2,s3,A2,A3)                    
                        KHI2=khi2_int(theta2bin,hist_norm,hist_err,Int_fitgauss (bin_edges[:-1],bin_edges[1:])/(bsize))
                        if(KHI2 > 2):
                            file_khi2toohigh.write(run_number+"\t"+str(E)+"\t"+str(zen)+"\t"+str(off)+"\t"+str(eff)+"\t"+str(len(theta2f))+"\t"+str(KHI2)+"\n")
                        #R68 for the fit and for the histogram distribution
                        #R68fit=R68(s1,s2,s3,A2,A3, theta2max)
                        #R68data=R68_hist(theta2f)

    outdir=os.path.expandvars('$HESSCONFIG')+"/"+config
    np.savez(outdir+"/PSF_triplegauss_"+config+".npz", TableSigma1=TableSigma1, TableSigma2=TableSigma2, TableSigma3=TableSigma3, TableA2=TableA2, TableA3=TableA3)
    file_toofewevents.close()
    file_nosimu.close()
