import numpy as np

class FrenchMcBands(object):

    def __init__(self,i=4):
        #First digit of the run: particule type
        self.part_type=["gamma"]
        self.part_type_digit="1"
        #Second digit of the run: azimuth angle
        self.azimuth=[180,0]
        self.azimuth_digit=["4", "0"]
        #third digit of the run: MC zenithal angle
        self.zenMC = [0, 18, 26, 32, 37, 41, 46, 50, 53, 57, 60, 63, 67, 70]
        self.zenMC_digit=["00","01","02","03","04","05","06","07","08","09","10","11","12","13"]
        #fouth digit of the run: MC offset 
        self.offMC = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
        self.offMC_digit=["0","1","2","3","4","5"]
        #fith digit of the run: energy type (spectrum or fixed energy)
        self.energy_type="1"
        #sixth digit of the run: MC energy
        self.enMC = [0.02, 0.03, 0.05, 0.08, 0.125, 0.2, 0.3, 0.5, 0.8, 1.25, 2, 3, 5, 8, 12.5, 20, 30, 50, 80, 125]
        self.enMC_digit=["01", "02" , "03" , "04", "05" ,"06" , "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]

    def find_azimuth_digit(self, az):
        """
        Return the digit associated to a given azimuth angle
        """
        i=np.where(np.asarray(self.azimuth)==az)[0]
        return self.azimuth_digit[i]
    
    def find_zenMC_digit(self, zen):
        """
        Return the digit associated to a given zenithal angle
        """
        i=np.where(np.asarray(self.zenMC)==zen)[0]
        return self.zenMC_digit[i]
    
    def find_offMC_digit(self, off):
        """
        Return the digit associated to a given offset
        """
        i=np.where(np.asarray(self.offMC)==off)[0]
        return self.offMC_digit[i]
    
    def find_enMC_digit(self, E):
        """
        Return the digit associated to a given energy
        """
        i=np.where(np.asarray(self.enMC)==E)[0]
        return self.enMC_digit[i]
    
    def run_number(self,az,zen, off, E):
        nen=self.find_enMC_digit(E)
        nzen=self.find_zenMC_digit(zen)
        noff=self.find_offMC_digit(off)
        naz=self.find_azimuth_digit(az)
        return self.part_type_digit+naz+nzen+noff+self.energy_type+nen
    
    def ener_MC(self, run_number):
        i=np.where(np.array(self.enMC_digit)==run_number[6:8])
        return self.enMC[i[0][0]]
    def zen_MC(self, run_number):
        i=np.where(np.array(self.zenMC_digit)==run_number[2:4])
        return self.zenMC[i[0][0]]
    def off_MC(self, run_number):
        i=np.where(np.array(self.offMC_digit)==run_number[4])
        return self.offMC[i[0][0]]
    def MC_values(self, run_number):
        E=self.ener_MC(run_number)
        zen=self.zen_MC(run_number)
        off=self.off_MC(run_number)
        print("The run"+run_number+" match with MCs: E="+str(E)+" TeV ,zen="+str(zen)+" deg ,off="+str(off)+" deg")
