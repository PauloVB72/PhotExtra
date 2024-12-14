import os
import requests
from utils import dowload_kernel,setup_directories
from utils import get_data
from utils import directory
from utils import folder_exists
from utils import survey_pixel_scale ,survey_resolution

from params import parameters
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve,convolve_fft
from astropy.wcs import WCS

from scipy.ndimage import zoom

class ConvolutionIMG:
    def __init__(self,survey:int,ker_survey:str,name:str,path=None,data=None):

        self.inp_surveys = survey
        self.ker_survey = ker_survey
        self.path = path
        self.name = name
        self.data = data
        
    def get_kernels(self):
        
        #workdir = os.getenv("workdir", "images")
        #if self.path == None:
         #   self.path = workdir

        #if folder_exists(self.path):
         #   obj_dir = os.path.join(self.path, 'kernels') 
          #  if not folder_exists(obj_dir):
           #     directory(obj_dir)
        #else:
         #   directory(self.path)
          #  obj_dir = os.path.join(self.path, 'kernels') 
           # directory(obj_dir)
        obj_dir = setup_directories(self.name,path=self.path)['main']
        obj_dir_ker = setup_directories(self.name,path=self.path)['kernels']


        kernel_name = get_data(self.inp_surveys, self.ker_survey)
        filename = os.path.join(obj_dir, kernel_name)

        if os.path.isfile(filename):
            print(f"Kernel already exists: {filename}")
            return filename
        filename = dowload_kernel(kernel_name, obj_dir_ker)
        print(filename)
        print(f"Kernels dowloaded")
        return filename


    def get_convolve(self):
        
        obj_dir = self.get_kernels()
        print('--------------------------------')
        print(obj_dir)

        if self.data is not None:
            data_inp = self.data
        else:
            path_inp = self.path+'/'+str(self.name)+'/images/'+str(self.inp_surveys)+'.fits'
            hdu_inp = fits.open(path_inp)
            data_inp = hdu_inp[0].data
            header_inp = hdu_inp[0].header
            wcs_inp = WCS(header_inp)

        path_ker = str(obj_dir)
        hdu_ker = fits.open(path_ker)
        data_ker = hdu_ker[0].data
        header_ker = hdu_ker[0].header

        self.ker_survey = self.ker_survey.split('_')[0]
        self.inp_surveys = self.inp_surveys.split('_')[0]

        #pxs_ker = float(survey_pixel_scale(self.ker_survey))
        #pxs_inp = float(survey_pixel_scale(self.inp_surveys))
        #res = float(survey_resolution(self.ker_survey))


        astropy_conv = convolve_fft(data_inp, data_ker, nan_treatment = 'interpolate', normalize_kernel=True,
                                preserve_nan=True)
        
       # import matplotlib.pyplot as plt
       # import numpy as np

        #fig, axs = plt.subplots(1, 3, figsize=(15, 5))

       # axs[0].imshow(data_inp, cmap='gray', origin='lower',vmin=np.nanmean(data_inp)-np.nanstd(data_inp),vmax=np.nanmean(data_inp)+np.nanstd(data_inp),
        #        interpolation='nearest')

       # axs[1].imshow(data_ker, cmap='gray', origin='lower',vmin=np.nanmean(data_ker)-np.nanstd(data_ker),vmax=np.nanmean(data_ker)+np.nanstd(data_ker),
       #         interpolation='nearest')
       # axs[2].imshow(astropy_conv, cmap='gray', origin='lower',vmin=np.nanmean(astropy_conv)-np.nanstd(astropy_conv),vmax=np.nanmean(astropy_conv)+np.nanstd(astropy_conv),
        #        interpolation='nearest')       
       # plt.show() 
        return astropy_conv
    

#gs = convolveimg('SDSS_g','unWISE_W2','/home/polo/Escritorio/Works/Doctorado/Code/SFHmergers/images','SIT45')
#gs.get_convolve()

