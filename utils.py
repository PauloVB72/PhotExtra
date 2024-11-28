import numpy as np
import pandas as pd
from astropy import units as u

config = pd.read_csv('prm_config.csv')


def survey_pixel_scale(survey:str):

    global config

    pxscale = config.loc[config['survey'] == survey, 'pixel_scale']
    
    return pxscale



def folder_exists(folder_path):
    return os.path.isdir(folder_path)


def directory(path:str)-> None:
            
    mdir = path
    if isinstance(path,str):
        try: 
            os.mkdir(mdir)
        except OSError:
            print ("Creation of the directory %s failed" % mdir)
        else:
            print ("Successfully created the directory %s " % mdir)