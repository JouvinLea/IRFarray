import numpy as np
import iminuit

class PSFfit(object):
    
    def __init__(self, theta2):
        self.theta2=theta2

    def nllp_triplegauss(self,s1,s2,s3,A2,A3):
        s12 = np.min([s1,s2,s3])**2
        s22 = np.median([s1,s2,s3])**2
        s32 = np.max([s1,s2,s3])**2
        norm = 2*np.pi*(s12+ np.abs(A2) * s22 + np.abs(A3) * s32)
        gaus1 = np.exp(-self.theta2/(2*s12))
        gaus2 = np.exp(-self.theta2/(2*s22))
        gaus3 = np.exp(-self.theta2/(2*s32))

        lnlike = np.log(gaus1 + A2*gaus2 + A3*gaus3) - np.log(norm)
        
        res = - lnlike.sum()

        return res

    def nllp_king(self,sig,gam):
        norm1=2*np.pi*sig**2
        norm2 = 1-1/gam
        king=(1+self.theta2/(2*gam*(sig**2)))**(-gam)
        lnlike = np.log(king) + np.log(norm2)-np.log(norm1)
        
        res = - lnlike.sum()

        return res

    def minimization(self, model_func, s1_init=0.02, s2_init=0.05,s3_init=0.08, A2_init=0.3,A3_init=0.1):
        if(model_func=="triplegauss"):
            m=iminuit.Minuit(self.nllp_triplegauss, s1=s1_init, s2=s2_init,s3=s3_init, A2=A2_init,A3=A3_init,
                                    limit_A2 = (1e-10,10.),limit_A3 = (1e-10,10.),
                                     limit_s1 = (0.005,0.3), limit_s2 = (0.005,0.6),
                                     limit_s3 = (0.02,0.7))
            m.migrad()
            s1_m=m.values['s1']
            s2_m=m.values['s2']
            s3_m=m.values['s3']
                    
            s1 = np.min([s1_m,s2_m,s3_m])
            s2 = np.median([s1_m,s2_m,s3_m])
            s3 = np.max([s1_m,s2_m,s3_m])
            return (s1, s2 , s3, m.values['A2'], m.values['A3'])        
        elif(model_func=="king"):
            m=iminuit.Minuit(self.nllp_king, sig=0.02, gam=2, limit_sig = (1e-10,10.),limit_gam = (1e-10,10.))
            #m=iminuit.Minuit(self.nllp_king, sig=0.07, gam=1.5, limit_sig = (1e-10,10.),limit_gam = (1e-10,10.))
            m.migrad()
            return (m.values['sig'],m.values['gam']) 
        else:
            print "Error: Vous n'avez pas donne de fonction..."

        
