import numpy as np
import pandas as pd
from astropy import units as u
import os
import ast



path = '/home/polo/Escritorio/Works/Doctorado/Code/SFHmergers/Photsfh/prm_config.csv'

def load_data(csv_file):
    df = pd.read_csv(csv_file)
    
    df['filters'] = df['filters'].apply(ast.literal_eval)
    return df

config = load_data(path)


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

    
print(check_filters('WISE',['W1','W2','W3']))