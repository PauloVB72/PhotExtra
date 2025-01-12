
import os
import pandas as pd
import pathlib
import numpy as np
from astropy import wcs, units as u
import warnings
from astropy.utils.exceptions import AstropyWarning


config = pd.read_csv('/home/polo/Escritorio/Works/Doctorado/Code/SFHmergers/Photsfh/prm_config.csv')

def check_data():
    pass

#data = ['SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z','2MASS_K','2MASS_J','2MASS_H','GALEX_FUV','GALEX_NUV','WISE_W1','WISE_W2','WISE_W3']

#position = (130.3,20.3)

#size = 10

class parameters:
    
    def __init__(self,data,position,size):
        # data separation
        
        if isinstance(data[0],str):
            surveys = {}
            for items in data:
                surv, band = items.split('_')
                if surv not in surveys:
                    surveys[surv] = []
                
                surveys[surv].append(band)
            self.surveys = surveys
            list_bands = list(surveys.values())
        else:
            warnings
        # position 
        if isinstance(position[0],(int,float)):
            pass
        else:
            warnings
        # size
        if size >= 5.:    
                warnings
        else:
            warnings
    def survey_values(self):

        return self.surveys
    
    def check_validity(self):
        global config
        # validation of the initial parameters
        for srv in self.surveys.keys():
            if srv in np.array(config['survey']):
                for flt in self.surveys[srv]:
                    if flt in ''.join(config.loc[config['survey'] == srv, 'filters'].astype(str)):
                        pass
                    else:
                        warnings
            else:
                warnings
    
