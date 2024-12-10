import numpy as np
import pandas as pd
from astropy import units as u
import os
import ast
import requests
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground

import json
from astropy.io import fits
from astropy.wcs import WCS



path = '/home/polo/Escritorio/Works/Doctorado/Code/SFHmergers/Photsfh/prm_config.csv'

def load_data(csv_file):
    df = pd.read_csv(csv_file)
    
    df['filters'] = df['filters'].apply(ast.literal_eval)
    return df

config = load_data(path)




def survey_resolution(survey:str):

    # mejorar para valores de distintas bandas

    global config
    
    resolution = config.loc[config['survey'] == survey, 'resolution']
    
    return resolution



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

def filter_check(survey,filter_in=None):
    
    global config

    filters = config.loc[config['survey'] == survey, 'filters']

    if survey in ['PS1','SDSS','LegacySurvey','UKIDSS']:

        if filter_in == None:
            filters = filters 
            return filters
        else:
            for vals in list(filter_in):
                if vals in list(filters):
                    pass
                else:
                    message = (
                        f"filter '{vals}' is not a valid option for "
                        f"'{survey}' survey ({filters})"
                    )
                    assert vals in filters, message
            return filter_in
    else:
        for vals in filter_in:
            if vals in filters:
                pass
            else:
                message = (
                        f"filter '{vals}' is not a valid option for "
                        f"'{survey}' survey ({filters})"
                    )
                assert vals in filters, message

        return filter_in

def check_filters(survey, filters):
    global config
    valid_surveys = {entry: filters for entry, filters in zip(config["survey"], config["filters"])}
    if survey not in valid_surveys:
        return (f"Survey '{survey}' is not a valid option.")
    
    valid_filters = valid_surveys[survey]
    
    if filters is None:
        if survey in ["PS1", "SDSS", "LegacySurvey", "UKIDSS"]:
           
            return ''.join(valid_filters)
        else:
            # Para otros surveys, devuelve la lista como está
            return valid_filters

    else:
        if survey in ["PS1", "SDSS", "LegacySurvey", "UKIDSS"]:
            if isinstance(filters, str):
                if all(len(f) == 1 for f in valid_filters):  
                    user_filters = list(filters)
                else:  
                    user_filters = filters.split(',')
            else:
                user_filters = filters

            if all(f in valid_filters for f in user_filters):
                return ''.join(user_filters)
            else:

                menssage = f"Error: Filters {user_filters} are not valid in the survey '{survey}'. Valid filters: {valid_filters}"
                assert False ,menssage
        else:
             
            if isinstance(filters, str):
                if all(len(f) == 1 for f in valid_filters):  
                    user_filters = list(filters)
                else:  
                    user_filters = filters.split(',')
            else:
                user_filters = filters

            if all(f in valid_filters for f in user_filters):
                return filters
            else:
                menssage = f"Error: Filters {user_filters} are not valid in the survey '{survey}'. Valid filters: {valid_filters}"
                assert False ,menssage       




def dowload_kernel(name:str,path:str):

    hi_res = "https://www.astro.princeton.edu/~draine/Kernels/Kernels_2018/Kernel_FITS_Files/Hi_Resolution/"

    file_url = hi_res + name    #"Kernel_HiRes_BiGauss_00.5_to_GALEX_FUV.fits.gz"

    output_folder = path

    file_name = os.path.join(output_folder, file_url.split("/")[-1])
    #print(f"Descargando {file_name}...")
    response = requests.get(file_url, stream=True)
    response.raise_for_status()
    print(response.iter_content(chunk_size=8192))
    with open(file_name, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    return file_name


def bkg_sub(data:np.array,survey:str):
    # A first approximation. at future change the bkg_stimator to election 
    bkg_surveys = ["2MASS", "WISE", "VISTA", "UKIDSS"]
    
    data = data.astype(np.float64)
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

    if survey in bkg_surveys:
        return np.copy(data - bkg)
    else:
        return  data


def header_changes(hdu,ra:float,dec:float,size_img:list,survey=None):
    if survey == 'SDSS':

        header = hdu.header
    else:
        header = hdu[0].header
    header['NAXIS1'] = size_img[0]
    header['NAXIS2'] = size_img[1]
    header['CRPIX1'] = size_img[0] / 2
    header['CRPIX2'] = size_img[1] / 2
    header['CRVAL1'] = ra
    header['CRVAL2'] = dec

    return hdu



def get_data(inp_survey: str, ker_survey:str):
    with open("/home/polo/Escritorio/Works/Doctorado/Code/SFHmergers/Photsfh/kernel_survey.json", "r") as f:
        data = json.load(f)

    # Buscar el valor específico para SDSS → WISE
    for entry in data:
        if entry["From"] == inp_survey and entry["To"] == ker_survey:
            print(f"Valor SDSS→WISE: {entry['Value']}")
            return entry['Value']
        
print(survey_pixel_scale('SDSS'))