from reproject import reproject_adaptive,reproject_exact,reproject_interp
from utils import dowload_kernel,setup_directories
from utils import get_data
from utils import directory
from utils import folder_exists
from utils import survey_pixel_scale ,survey_resolution
from params import parameters
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve,convolve_fft
from astropy.wcs import WCS


class ReprojectIMG:
    def __init__(self,survey:int,ref_survey:str,name:str,path=None,data=None):

        self.inp_surveys = survey
        self.ref_survey = ref_survey
        self.path = path
        self.name = name
        self.data = data

    def get_reproject(self):

        if self.data is not None:
                hdu_inp = self.data
                data_inp = self.data[0]
        else:
                obj_dir = setup_directories(self.name,path=self.path)['images']
                path_inp = obj_dir+'/'+str(self.inp_surveys)+'.fits'
                hdu_inp = fits.open(path_inp)
                data_inp = hdu_inp[0].data
                header_inp = hdu_inp[0].header

        obj_dir = setup_directories(self.name,path=self.path)['images']
        path_ref = obj_dir+'/'+str(self.ref_survey)+'.fits'
        hdu_ref = fits.open(path_ref)
        data_ref = hdu_ref[0].data
        header_ref = hdu_ref[0].header
        wcs_ref = WCS(header_ref)

        self.ref_survey = self.ref_survey.split('_')[0]
        self.inp_surveys = self.inp_surveys.split('_')[0]

       # pxs_ker = float(survey_pixel_scale(self.ref_survey))
       # pxs_inp = float(survey_pixel_scale(self.inp_surveys))

        reprojection_data , footprint = reproject_exact(hdu_inp,wcs_ref)

        import matplotlib.pyplot as plt
        import numpy as np


        return reprojection_data



#gs = ReprojectIMG('SDSS_g','unWISE_W2','/home/polo/Escritorio/Works/Doctorado/Code/SFHmergers/images','SIT45')
#gs.get_reproject()