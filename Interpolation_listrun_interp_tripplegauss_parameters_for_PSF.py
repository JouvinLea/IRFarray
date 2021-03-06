#! /usr/bin/env python
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy import interpolate
import math
from astropy.io.fits import Column
import sys
import os
from glob import glob
from pathlib import Path
import os
import argparse

"""
Commande a lancer pour pouvoir donner des arguments au scripts
"""


# ./Interpolation_listrun.py 'Crab_All.list' 'ash_stereo' "Prod15_4_stereo"
#./Interpolation_listrun.py 'Crab_All.list' 'ash_stereo_thsq64' "Prod15_4_stereo"

class Observation:
    """Helper functions to compute file and folder names.
    """

    # filetypes = ['events', 'aeff', 'edisp', 'psf_3gauss']
    filetypes = ['events']

    def __init__(self, obs_id, hap_config=None, telpattern=None):
        self.obs_id = obs_id
        self.hap_config = hap_config
        self.telpattern = telpattern

    @property
    def obs_group(self):
        obs_id_min = self.obs_id - (self.obs_id % 200)
        obs_id_max = obs_id_min + 199
        return obs_id_min, obs_id_max

    @property
    def _obs_group_folder(self):
        return Path('run{:06d}-{:06d}'.format(self.obs_group[0], self.obs_group[1]))

    @property
    def _obs_folder(self):
        return Path('run{:06d}'.format(self.obs_id))

    def folder(self, step=None):
        """Create folder for a given step.
        """
        if step is None:
            return self._obs_group_folder / self._obs_folder
        else:
            return Path(step) / self._obs_group_folder / self._obs_folder

    def hap_filename(self, filetype):
        """Name of FITS file generated by HAP"""
        if filetype == 'events':
            return self.folder('events') / 'run_{:07d}_{}_eventlist.fits'.format(self.obs_id, self.hap_config)
            # return self.folder('events') / 'events_{:06d}.fits.gz'.format(self.obs_id)
        elif filetype == 'aeff':
            return self.folder('irfs') / 'aeff_{:06d}.fits.gz'.format(self.obs_id)
        elif filetype == 'edisp':
            return self.folder('irfs') / 'edisp_{:06d}.fits.gz'.format(self.obs_id)
        elif filetype == 'psf_3gauss':
            return self.folder('irfs') / 'psf_3gauss_{:06d}.fits.gz'.format(self.obs_id)
        else:
            raise ValueError('Invalid {} {}'.format(filetype))

    def out_filename(self, filetype, dir, format='old'):
        """Name of FITS file in out folder"""
        filename = self.filename(filetype=filetype, format=format)
        return Path(dir) / filename

    def filename(self, filetype, format='old'):
        if format == 'old':
            TAGS = dict(
                events='events',
                aeff='aeff_2d',
                edisp='edisp_2d',
                psf_3gauss='psf_3gauss',
                psf_king='psf_king',
                psf_table='psf_table',
                background='bkg_offruns',
            )
        elif format == 'new':
            TAGS = dict(
                events='events',
                aeff='aeff',
                edisp='edisp',
                psf_3gauss='psf_3gauss',
                psf_king='psf_king',
                psf_table='psf_table',
                background='background',
            )

        tag = TAGS[filetype]
        filename = '{}_{:06d}.fits.gz'.format(tag, self.obs_id)
        return self.folder() / filename

    def mkdir(self, step):
        """Make directory (parts=True, exists_ok=True)"""
        path = self.folder(step)
        if not path.exists():
            path.mkdir(parents=True)

        return path

    def check_out_files_exist(self):
        """Check if all out files exist"""
        for filetype in self.filetypes:
            filename = self.out_filename(filetype)
            if not filename.is_file():
                log.error('MISSING: {}'.format(filename))
                return False

        return True


def gauss(x, sigma, mean):
    f = 1 / (np.sqrt(2 * math.pi) * sigma) * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))
    return f


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make the index and obs table')
    parser.add_argument('runlist', action="store",
                        help='List of run for which we want to interpolate the IRFs')
    parser.add_argument('config', action="store", help='Prod Configuration, cut we apply')
    parser.add_argument('prod', action="store", help='Prod')
    arg = parser.parse_args()
    runlist = arg.runlist
    coupure = arg.config
    prod = arg.prod
    directory = os.path.expandvars('$CALDB')
    PathListRun = directory + "/" + prod + "/" + coupure

    # Sur Lyon
    # PathListRun = directory+"/data/hess/hap-16-03_fits/"+prod+"/"+config

    RunNumber = np.loadtxt(runlist, ndmin=2)
    # Load les info sur les MCs depuis la table d'IRF ou est stocke pour toutes les nergies, zenith,offset et efficacite des MCs la valeur des la surface efficiace, du biais et sigma pour la resolution et du s1, s2, s3, A2, A3 de la tripplegauss utilisee pour fitter la psf
    PathTableIRF = os.path.expandvars('$HESSCONFIG')
    PathTablePSF = os.path.expandvars('$HESSCONFIG')

    for nrun in RunNumber[:, 0]:
        obs = Observation(int(nrun))
        informat = "old"
        namerun = str(Path(PathListRun) / obs.filename('events', format=informat))
        try:
            table = Table.read(namerun, hdu='EVENTS')
        except Exception:
            print "fits corrupted for file " + namerun
            continue
        hdurun = fits.open(namerun)
        az=hdurun[1].data["AZ"]
        az[np.where(az>180)]=az[np.where(az>180)]-360
        AZRun=az.mean()
        AltRun = hdurun[1].data["ALT"].mean()
        if ((AZRun > 90) & (AZRun < 270)):
            mode = "south"
        else:
            mode = "north"
        ZenRun = 90 - AltRun
        EffRun = hdurun[1].header["MUONEFF"] * 100
        name_config = coupure[0:3] + "_" + mode + "_" + coupure[4:]
        print(nrun)
        print(PathTableIRF + "/" + name_config + "/IRF_" + name_config + ".npz")
        IRF = np.load(PathTableIRF + "/" + name_config + "/IRF_" + name_config + ".npz")
        IRFArea = IRF["TableArea"]
        IRFSigma = IRF["TableSigma"]
        IRFBiais = IRF["TableBiais"]
        enMC = IRF["enMC"]
        lnenMC = IRF["lnenMC"]
        zenMC = IRF["zenMC"]
        effMC = IRF["effMC"]
        offMC = IRF["offMC"]

        PSF = np.load(PathTablePSF + "/" + name_config + "/PSF_triplegauss_" + name_config + ".npz")
        PSFs1 = PSF["TableSigma1"]
        PSFs2 = PSF["TableSigma2"]
        PSFs3 = PSF["TableSigma3"]
        PSFA2 = PSF["TableA2"]
        PSFA3 = PSF["TableA3"]

        binoffMC = len(offMC)
        binEMC = len(enMC)
        # binEreco=50
        binEreco = 100
        bineffarea = len(offMC) * len(enMC)
        bineffresol = len(offMC) * len(enMC) * binEreco

        # reverifier qu ils ont bien ca dans leur bin PA en low edge et upper edge
        off_low = offMC
        off_hi = offMC

        # pour les extremites prendre le milieu des bin en log
        binlnEMC = lnenMC[1:] - lnenMC[:-1]
        # Pour le premier bin en energie pour defenr le edge low du bin on prend la demilargeur du premier bin
        binlnEMClow = np.insert(binlnEMC, 0, binlnEMC[0])
        # Pour le dernier bin en energie pour defenr le edge up du bin on prend la demilargeur du dernier bin
        binlnEMCup = np.insert(binlnEMC, -1, binlnEMC[-1])
        # Retrouver
        lnEMClow = lnenMC - binlnEMClow / 2
        lnEMCup = lnenMC + binlnEMCup / 2
        E_true_low = pow(10, lnEMClow)
        E_true_up = pow(10, lnEMCup)

        # Definition de Etrue/Ereco
        # lnEreco_true=np.linspace(-1,1,binEreco)
        lnEreco_true = np.linspace(-2, 2, binEreco)
        # Le tableau d energie reco en log ont tous la meme largeur de bin donc on prend le premier
        binlnEreco_true = lnEreco_true[1] - lnEreco_true[0]
        lnE_reco_true_low = lnEreco_true - binlnEreco_true / 2
        lnE_reco_true_up = lnEreco_true + binlnEreco_true / 2
        Ereco_true = np.exp(lnEreco_true)
        E_reco_true_low = np.exp(lnE_reco_true_low)
        E_reco_true_hi = np.exp(lnE_reco_true_up)

        AreaRun = np.zeros((binoffMC, binEMC))
        ResolRun = np.zeros((binoffMC, binEreco, binEMC))

        PSFS1Run = np.zeros((binoffMC, binEMC))
        PSFS2Run = np.zeros((binoffMC, binEMC))
        PSFS3Run = np.zeros((binoffMC, binEMC))
        PSFA2Run = np.zeros((binoffMC, binEMC))
        PSFA3Run = np.zeros((binoffMC, binEMC))

        for (iEMC, EMC) in enumerate(enMC):
            for (ioff, off) in enumerate(offMC):
                # print ioff, " ", iEMC
                InterArea = interpolate.interp2d(effMC, np.cos(zenMC * math.pi / 180), IRFArea[iEMC, ioff, :, :])
                InterBiais = interpolate.interp2d(effMC, np.cos(zenMC * math.pi / 180), IRFBiais[iEMC, ioff, :, :])
                InterSigma = interpolate.interp2d(effMC, np.cos(zenMC * math.pi / 180), IRFSigma[iEMC, ioff, :, :])
                AreaRun[ioff, iEMC] = InterArea(EffRun, np.cos(ZenRun * math.pi / 180))
                BiaisRun = InterBiais(EffRun, np.cos(ZenRun * math.pi / 180))
                SigmaRun = InterSigma(EffRun, np.cos(ZenRun * math.pi / 180))
                ResolRun[ioff, :, iEMC] = gauss(lnEreco_true, SigmaRun, BiaisRun) / Ereco_true
                # etre sur que c est bien normalise
                norm = np.sum(ResolRun[ioff, :, iEMC] * (E_reco_true_hi - E_reco_true_low))

                if (np.isnan(norm)):
                    ResolRun[ioff, :, iEMC] = 0
                else:
                    ResolRun[ioff, :, iEMC] = ResolRun[ioff, :, iEMC] / norm

                ind_zen, ind_eff = np.where(PSFs1[iEMC, ioff, :, :] != -1)
                # If there is at least one simu for this offset and this energy for wich the fit works
                if (len(ind_zen) != 0):
                    zensame = np.where(ind_zen != ind_zen[0])
                    effsame = np.where(ind_eff != ind_eff[0])
                    # Il doit y avoir au moins 2 valeurs differentes en efficacite et en zenith pour que l interpolateur marche
                    if ((len(zensame[0]) != 0) & (len(effsame[0]) != 0)):
                        coord_eff = effMC[ind_eff]
                        coord_zen = zenMC[ind_zen]
                        points = (coord_eff, np.cos(coord_zen * math.pi / 180))

                        PSFS1Run[ioff, iEMC] = interpolate.griddata(points, PSFs1[iEMC, ioff, ind_zen, ind_eff],
                                                                    (EffRun, np.cos(ZenRun * math.pi / 180)),
                                                                    method='linear')
                        if np.isnan(PSFS1Run[ioff, iEMC]):
                            PSFS1Run[ioff, iEMC] = interpolate.griddata(points, PSFs1[iEMC, ioff, ind_zen, ind_eff],
                                                                        (EffRun, np.cos(ZenRun * math.pi / 180)),
                                                                        method='nearest')

                        PSFS2Run[ioff, iEMC] = interpolate.griddata(points, PSFs2[iEMC, ioff, ind_zen, ind_eff],
                                                                    (EffRun, np.cos(ZenRun * math.pi / 180)),
                                                                    method='linear')
                        if np.isnan(PSFS2Run[ioff, iEMC]):
                            PSFS2Run[ioff, iEMC] = interpolate.griddata(points, PSFs2[iEMC, ioff, ind_zen, ind_eff],
                                                                        (EffRun, np.cos(ZenRun * math.pi / 180)),
                                                                        method='nearest')

                        PSFS3Run[ioff, iEMC] = interpolate.griddata(points, PSFs3[iEMC, ioff, ind_zen, ind_eff],
                                                                    (EffRun, np.cos(ZenRun * math.pi / 180)),
                                                                    method='linear')
                        if np.isnan(PSFS3Run[ioff, iEMC]):
                            PSFS3Run[ioff, iEMC] = interpolate.griddata(points, PSFs3[iEMC, ioff, ind_zen, ind_eff],
                                                                        (EffRun, np.cos(ZenRun * math.pi / 180)),
                                                                        method='nearest')

                        PSFA2Run[ioff, iEMC] = interpolate.griddata(points, PSFA2[iEMC, ioff, ind_zen, ind_eff],
                                                                    (EffRun, np.cos(ZenRun * math.pi / 180)),
                                                                    method='linear')
                        if np.isnan(PSFA2Run[ioff, iEMC]):
                            PSFA2Run[ioff, iEMC] = interpolate.griddata(points, PSFA2[iEMC, ioff, ind_zen, ind_eff],
                                                                        (EffRun, np.cos(ZenRun * math.pi / 180)),
                                                                        method='nearest')

                        PSFA3Run[ioff, iEMC] = interpolate.griddata(points, PSFA3[iEMC, ioff, ind_zen, ind_eff],
                                                                    (EffRun, np.cos(ZenRun * math.pi / 180)),
                                                                    method='linear')
                        if np.isnan(PSFA3Run[ioff, iEMC]):
                            PSFA3Run[ioff, iEMC] = interpolate.griddata(points, PSFA3[iEMC, ioff, ind_zen, ind_eff],
                                                                        (EffRun, np.cos(ZenRun * math.pi / 180)),
                                                                        method='nearest')

                            #                else:
                            #                    PSFS1Run[ioff, iEMC] = -1
                            #                    PSFS2Run[ioff, iEMC] = -1
                            #                    PSFS3Run[ioff, iEMC] = -1
                            #                    PSFA2Run[ioff, iEMC] = -1
                            #                    PSFA3Run[ioff, iEMC] = -1
                            #            else:
                            #                PSFS1Run[ioff, iEMC] = -1
                            #                PSFS2Run[ioff, iEMC] = -1
                            #                PSFS3Run[ioff, iEMC] = -1
                            #                PSFA2Run[ioff, iEMC] = -1
                            #                PSFA3Run[ioff, iEMC] = -1

        outdir = str(Path(PathListRun) / obs.folder())
        # Ecriture des fichiers fits pour aeff, edisp et psf pour chaque observation
        # AEFF FITS FILE
        c1_area = Column(name='ENERG_LO', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_low))
        c2_area = Column(name='ENERG_HI', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_up))
        c3_area = Column(name='THETA_LO', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_low))
        c4_area = Column(name='THETA_HI', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_hi))
        c5_area = Column(name='EFFAREA', format=str(bineffarea) + 'E', unit='m2', array=np.expand_dims(AreaRun, 0))
        c6_area = Column(name='EFFAREA_RECO', format=str(bineffarea) + 'E', unit='m2', array=np.expand_dims(AreaRun, 0))
        tbhdu_area = fits.BinTableHDU.from_columns([c1_area, c2_area, c3_area, c4_area, c5_area, c6_area])
        for i in range(1, 7):
            tbhdu_area.header.comments['TTYPE' + str(i)] = 'label for field ' + str(i)
            tbhdu_area.header.comments['TFORM' + str(i)] = 'data format of field: 4-byte REAL'
            tbhdu_area.header.comments['TUNIT' + str(i)] = 'physical unit of field '

        tbhdu_area.header.set("EXTNAME", "EFFECTIVE AREA", "name of this binary table extension ")
        tbhdu_area.header.set("TDIM5", "(" + str(binEMC) + "," + str(binoffMC) + ")")
        tbhdu_area.header.set("TDIM6", "(" + str(binEMC) + "," + str(binoffMC) + ")")
        tbhdu_area.header.set("LO_THRES", -1, "TeV")
        tbhdu_area.header.set("HI_THRES", 150, "TeV")
        # tbhdu_area.header["EXTNAME"]='EFFECTIVE AREA'
        #tbhdu_area.writeto(outdir + '/aeff_2d_0' + str(int(nrun)) + '.fits', clobber=True)
        tbhdu_area.writeto(outdir + '/aeff_2d_{:06d}.fits'.format(int(nrun)), clobber=True)
        if Path(outdir + '/hess_aeff_2d_' + str(int(nrun)) + '.fits').exists():
            os.remove(outdir + '/hess_aeff_2d_' + str(int(nrun)) + '.fits')
        if Path(outdir + '/hess_aeff_2d_0' + str(int(nrun)) + '.fits').exists():
            os.remove(outdir + '/hess_aeff_2d_0' + str(int(nrun)) + '.fits')

        # EDISP FITS FILE
        c1_resol = Column(name='ETRUE_LO', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_low))
        c2_resol = Column(name='ETRUE_HI', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_up))
        c3_resol = Column(name='MIGRA_LO', format=str(binEreco) + 'E', unit='', array=np.atleast_2d(E_reco_true_low))
        c4_resol = Column(name='MIGRA_HI', format=str(binEreco) + 'E', unit='', array=np.atleast_2d(E_reco_true_hi))
        c5_resol = Column(name='THETA_LO', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_low))
        c6_resol = Column(name='THETA_HI', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_hi))
        c7_resol = Column(name='MATRIX ', format=str(bineffresol) + 'E', unit='TeV', array=np.expand_dims(ResolRun, 0))
        tbhdu_resol = fits.BinTableHDU.from_columns(
            [c1_resol, c2_resol, c3_resol, c4_resol, c5_resol, c6_resol, c7_resol])
        for i in range(1, 8):
            tbhdu_resol.header.comments['TTYPE' + str(i)] = 'label for field ' + str(i)
            tbhdu_resol.header.comments['TFORM' + str(i)] = 'data format of field: 4-byte REAL'
            #        tbhdu_resol.header.comments['TUNIT'+str(i)]='physical unit of field '

        tbhdu_resol.header.set("EXTNAME", "EDISP_2D", "name of this binary table extension ")
        tbhdu_resol.header.set("TDIM7", "(" + str(binEMC) + "," + str(binEreco) + "," + str(binoffMC) + ")")
        # tbhdu_resol.header["EXTNAME"]='EFFECTIVE RESOL'
        tbhdu_resol.writeto(outdir + '/edisp_2d_{:06d}.fits'.format(int(nrun)), clobber=True)
        #tbhdu_resol.writeto(outdir + '/edisp_2d_0' + str(int(nrun)) + '.fits', clobber=True)
        if Path(outdir + '/hess_edisp_2d_' + str(int(nrun)) + '.fits').exists():
            os.remove(outdir + '/hess_edisp_2d_' + str(int(nrun)) + '.fits')
        if Path(outdir + '/hess_edisp_2d_0' + str(int(nrun)) + '.fits').exists():
            os.remove(outdir + '/hess_edisp_2d_0' + str(int(nrun)) + '.fits')
        # PSF FITS FILE
        c1_psf = Column(name='ENERG_LO', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_low))
        c2_psf = Column(name='ENERG_HI', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_up))
        c3_psf = Column(name='THETA_LO', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_low))
        c4_psf = Column(name='THETA_HI', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_hi))

        norm = 2 * np.pi * (PSFS1Run ** 2 + PSFA2Run * PSFS2Run ** 2 + PSFA3Run * PSFS3Run ** 2)
        c5_psf = Column(name='SIGMA_1', format=str(bineffarea) + 'E', unit='deg', array=np.expand_dims(PSFS1Run, 0))
        c6_psf = Column(name='AMPL_2', format=str(bineffarea) + 'E', unit='', array=np.expand_dims(PSFA2Run, 0))
        c7_psf = Column(name='SIGMA_2', format=str(bineffarea) + 'E', unit='deg', array=np.expand_dims(PSFS2Run, 0))
        c8_psf = Column(name='AMPL_3', format=str(bineffarea) + 'E', unit='', array=np.expand_dims(PSFA3Run, 0))
        c9_psf = Column(name='SIGMA_3', format=str(bineffarea) + 'E', unit='deg', array=np.expand_dims(PSFS3Run, 0))
        c10_psf = Column(name='SCALE', format=str(bineffarea) + 'E', unit='', array=np.expand_dims(1 / norm, 0))
        tbhdu_psf = fits.BinTableHDU.from_columns(
            [c1_psf, c2_psf, c3_psf, c4_psf, c5_psf, c6_psf, c7_psf, c8_psf, c9_psf, c10_psf])
        tbhdu_psf.header.set("EXTNAME", "PSF_2D", "name of this binary table extension")
        tbhdu_psf.writeto(outdir + '/psf_3gauss_{:06d}.fits'.format(int(nrun)), clobber=True)
        #tbhdu_psf.writeto(outdir + '/psf_3gauss_0' + str(int(nrun)) + '.fits', clobber=True)
        if Path(outdir + '/hess_psf_3gauss_' + str(int(nrun)) + '.fits').exists():
            os.remove(outdir + '/hess_psf_3gauss_' + str(int(nrun)) + '.fits')
        if Path(outdir + '/hess_psf_3gauss_0' + str(int(nrun)) + '.fits').exists():
            os.remove(outdir + '/hess_psf_3gauss_0' + str(int(nrun)) + '.fits')
