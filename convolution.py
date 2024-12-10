import os
import requests
from utils import dowload_kernel
from utils import get_data
from utils import directory
from utils import folder_exists
from utils import survey_pixel_scale ,survey_resolution

from params import parameters
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve,convolve_fft
from astropy.wcs import WCS

from scipy.ndimage import zoom

class convolveimg:
    def __init__(self,survey:int,ker_survey:str,path:str,name:str):

        self.inp_surveys = survey
        self.ker_survey = ker_survey
        self.path = path
        self.name = name

    def get_kernels(self):
        
        workdir = os.getenv("workdir", "images")
        if self.path == None:
            self.path = workdir
        else:
            pass

        if folder_exists(self.path) is True:

            obj_dir = os.path.join(self.path, 'kernels') 
            if folder_exists(obj_dir) is True:
                pass
            else:
                directory(obj_dir)
        else:
            directory(self.path)
            obj_dir = os.path.join(self.path, 'kernels') 
            directory(obj_dir)

        kernel_name = get_data(self.inp_surveys,self.ker_survey)
        
        filename = dowload_kernel(kernel_name,obj_dir)

        print(f"Kernels dowloaded")
        return filename


    def get_convolve(self):
        
        obj_dir = self.get_kernels()

        path_inp = self.path+'/'+str(self.name)+'/'+str(self.inp_surveys)+'.fits'
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


        pxs_ker = float(survey_pixel_scale(self.ker_survey))
        pxs_inp = float(survey_pixel_scale(self.inp_surveys))
        res = float(survey_resolution(self.ker_survey))

        print(pxs_ker,pxs_inp)

        if pxs_ker != pxs_inp:
            ratio = pxs_ker/pxs_inp
            size = ratio * data_ker.shape[0]
            if round(size) % 2 == 0:
                size += 1
                ratio = size / data_ker.shape[0]
        else:
            ratio = 1.
        print(ratio)
        newkernel = zoom(data_ker, ratio) / ratio**2
        print(res)
        astropy_conv = convolve_fft(data_inp, newkernel, nan_treatment = 'interpolate', normalize_kernel=True,
                                preserve_nan=True)
        
        import matplotlib.pyplot as plt
        import numpy as np

        fig, axs = plt.subplots(1, 3, figsize=(15, 5))

        axs[0].imshow(data_inp, cmap='gray', origin='lower',vmin=np.nanmean(data_inp)-np.nanstd(data_inp),vmax=np.nanmean(data_inp)+np.nanstd(data_inp),
                interpolation='nearest')

        axs[1].imshow(newkernel, cmap='gray', origin='lower',vmin=np.nanmean(newkernel)-np.nanstd(newkernel),vmax=np.nanmean(newkernel)+np.nanstd(newkernel),
                interpolation='nearest')
        axs[2].imshow(astropy_conv, cmap='gray', origin='lower',vmin=np.nanmean(astropy_conv)-np.nanstd(astropy_conv),vmax=np.nanmean(astropy_conv)+np.nanstd(astropy_conv),
                interpolation='nearest')       
        plt.show() 
        return astropy_conv
    

gs = convolveimg('SDSS_g','WISE_W1','/home/polo/Escritorio/Works/Doctorado/Code/SFHmergers/images','SIT45')
gs.get_convolve()

