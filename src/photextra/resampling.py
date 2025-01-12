from reproject import reproject_adaptive,reproject_exact,reproject_interp
from astropy.io import fits
from astropy.wcs import WCS

from .utils import setup_directories

import logging

class ReprojectIMG:
    def __init__(self,survey:int,ref_survey:str,name:str,path=None,data=None,method='exact'):
        """
        Initialize the ReprojectIMG class.

        Parameters:
        survey (int): The survey to reproject.
        ref_survey (str): The reference survey.
        name (str): The name of the object.
        path (str): The path to the data. Defaults to None.
        data (fits.HDUList): The data to reproject. Defaults to None.
        method (str): The resampling method to use. Defaults to 'exact'.
        """
        self.inp_surveys = survey
        self.ref_survey = ref_survey
        self.path = path
        self.name = name
        self.data = data
        self.method = method

    def get_reproject(self,**kwargs):

        if self.data is not None:
                hdu_inp = self.data

        else:
                obj_dir = setup_directories(self.name,path=self.path)['images']
                path_inp = obj_dir+'/'+str(self.inp_surveys)+'.fits'
                hdu_inp = fits.open(path_inp)

        obj_dir = setup_directories(self.name,path=self.path)['images']
        path_ref = obj_dir+'/'+str(self.ref_survey)+'.fits'
        hdu_ref = fits.open(path_ref)
        header_ref = hdu_ref[0].header
        wcs_ref = WCS(header_ref)

        self.ref_survey = self.ref_survey.split('_')[0]
        self.inp_surveys = self.inp_surveys.split('_')[0]

       # pxs_ker = float(survey_pixel_scale(self.ref_survey))
       # pxs_inp = float(survey_pixel_scale(self.inp_surveys))
        if self.method == 'adaptive':
            reprojection_data, _ = reproject_adaptive(hdu_inp, wcs_ref,**kwargs)
        elif self.method == 'exact':
            reprojection_data, _ = reproject_exact(hdu_inp, wcs_ref,**kwargs)
        elif self.method == 'interp':
            reprojection_data, _ = reproject_interp(hdu_inp, wcs_ref,**kwargs)
        else:
            logging.error("Invalid resampling method. Please choose from 'adaptive', 'exact', or 'interp'.")
            return None

        return reprojection_data




#gs = ReprojectIMG('SDSS_g','unWISE_W2','/home/polo/Escritorio/Works/Doctorado/Code/SFHmergers/images','SIT45')
#gs.get_reproject()