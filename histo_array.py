#! /usr/bin/env python
import numpy as np
import ROOT
import root_numpy as nr
from root_numpy import hist2array
from root_numpy.testdata import get_filepath
import argparse

"""
For one specific config, store the IRF on 4D numpy table for each value of the Zenithal angle, Offset, Efficiency and Energy used for the MCs simulation
Example of commande line to run to create this 4D table with the directory of the MC simulation output and the config name as argument
./histo_array.py '/Users/jouvin/Desktop/these/WorkGAMMAPI/IRF/Brunoconfig' 'elm_south_stereo_Prod15_5'
./histo_array.py '/Users/jouvin/Desktop/these/WorkGAMMAPI/IRF/Brunoconfig' 'elm_north_stereo_Prod15_5'
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Store the IRF from Mc simulation in a 4D numpy table')
    parser.add_argument('directory', action="store", help='directory of the IRFs obtained from the MCs simulation')
    parser.add_argument('config', action="store", help='Config')
    results = parser.parse_args()
    print "Store the IRF in a 4D table from the MC simulations in ", results.directory , " and for the config ", results.config

    # Zenithal angle, offset, efficiency and energy used to simulate the MCs in HAP-fr
    enMC = [0.02, 0.03, 0.05, 0.08, 0.125, 0.2, 0.3, 0.5, 0.8, 1.25, 2, 3, 5, 8, 12.5, 20, 30, 50, 80, 125]
    lnenMC = np.log10(enMC)
    zenMC = [0, 18, 26, 32, 37, 41, 46, 50, 53, 57, 60, 63, 67, 70]
    effMC = [50, 60, 70, 80, 90, 100]
    offMC = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]

    binEMC = len(enMC)
    binzen = len(zenMC)
    binoff = len(offMC)
    bineff = len(effMC)

    # Size of the 4D table where we will stock the effective area and the resolution for each true energy, zenithal angle, offset and efficiency of the MCs
    TableArea = np.zeros((binEMC, binoff, binzen, bineff))
    TableBiais = np.zeros((binEMC, binoff, binzen, bineff))
    TableSigma = np.zeros((binEMC, binoff, binzen, bineff))

    # In Hap-FR: 20 true energies used for the MC simulation
    # In the .root histogram for the effective area and the resolution there are more than 20 energy bins so we take one MC simulation and we determine to which index in these histograms match the MC energy.
    filehistoarea = results.directory + "/" + results.config + "/gFixedEnergy_paris_0-8-8-8_CamOptimal_hwtrig_eff100/CollectionArea.root"
    filehistoresol = results.directory + "/" + results.config + "/gFixedEnergy_paris_0-8-8-8_CamOptimal_hwtrig_eff100/EnergyResolution.root"
    name_hist_area = "EffArea_67deg_2.5off_eff100_FixedE"
    name_hist_biais = "Resol_Biais_67deg_2.5off_eff100_FixedE"
    TFileArea = ROOT.TFile(filehistoarea)
    TFileResol = ROOT.TFile(filehistoresol)
    histoarea = TFileArea.Get(name_hist_area)
    historesol = TFileResol.Get(name_hist_biais)
    Ehistoarea = [histoarea.GetXaxis().GetBinLowEdge(i) for i in range(1, 22)]
    Ehistoresol = [historesol.GetXaxis().GetBinLowEdge(i) for i in range(1, 42)]
    ind_area = [np.where(Ehistoarea < i)[0][-1] for i in lnenMC]
    ind_resol = [np.where(Ehistoresol < i)[0][-1] for i in lnenMC]

    #Tranform ROOT histogram into numpy array and store the results for each MC simultation in the 4D table
    for (ieff, eff) in enumerate(effMC):
        print eff
        for (ioff, off) in enumerate(offMC):
            print off
            for (izen, zen) in enumerate(zenMC):
                print zen
                filehistoarea = results.directory + "/" + results.config + "/gFixedEnergy_paris_0-8-8-8_CamOptimal_hwtrig_eff" + str(
                    eff) + "/CollectionArea.root"
                filehistoresol = results.directory + "/" + results.config + "/gFixedEnergy_paris_0-8-8-8_CamOptimal_hwtrig_eff" + str(
                    eff) + "/EnergyResolution.root"
                TFileArea = ROOT.TFile(filehistoarea)
                TFileResol = ROOT.TFile(filehistoresol)
                if (izen == 0):
                    name_hist_area = "EffArea_00deg_" + str(off) + "off_eff" + str(eff) + "_FixedE"
                    name_hist_biais = "Resol_Biais_00deg_" + str(off) + "off_eff" + str(eff) + "_FixedE"
                    name_hist_sigma = "Resol_Sigma_00deg_" + str(off) + "off_eff" + str(eff) + "_FixedE"
                else:
                    name_hist_area = "EffArea_" + str(zen) + "deg_" + str(off) + "off_eff" + str(eff) + "_FixedE"
                    name_hist_biais = "Resol_Biais_" + str(zen) + "deg_" + str(off) + "off_eff" + str(eff) + "_FixedE"
                    name_hist_sigma = "Resol_Sigma_" + str(zen) + "deg_" + str(off) + "off_eff" + str(eff) + "_FixedE"
                histoarea = TFileArea.Get(name_hist_area)
                histobiais = TFileResol.Get(name_hist_biais)
                histosigma = TFileResol.Get(name_hist_sigma)
                AreaArray = hist2array(histoarea)[ind_area]
                AreaBiais = hist2array(histobiais)[ind_resol]
                AreaSigma = hist2array(histosigma)[ind_resol]
                TableArea[:, ioff, izen, ieff] = AreaArray
                TableBiais[:, ioff, izen, ieff] = AreaBiais
                TableSigma[:, ioff, izen, ieff] = AreaSigma
    outdir="output_4Dnumpyarrays"
    np.savez(outdir+"/IRF_"+results.config+".npz", TableArea=TableArea, TableBiais=TableBiais, TableSigma=TableSigma, enMC=enMC, lnenMC=lnenMC,
             zenMC=zenMC, offMC=offMC, effMC=effMC)
