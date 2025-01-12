import os
from .utils import download_kernel,setup_directories
from .utils import get_data
from .utils import pxscale

from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve,convolve_fft
from scipy.ndimage import zoom

class ConvolutionIMG:
    def __init__(self,survey:int,ker_survey:str,name:str,path=None,hdul=None,base='high_res'):

        self.inp_surveys = survey
        self.ker_survey = ker_survey
        self.path = path
        self.name = name
        self.hdul = hdul
        self.base = base
        
    def get_kernels(self):
        
        obj_dir = setup_directories(self.name,path=self.path)['main']
        obj_dir_ker = setup_directories(self.name,path=self.path)['kernels']
        kernel_name = get_data(self.inp_surveys, self.ker_survey)
        filename = os.path.join(obj_dir, kernel_name)

        if os.path.isfile(filename):
            print(f"Kernel already exists: {filename}")
            return filename
        else:
            filename = download_kernel(kernel_name, obj_dir_ker,base=self.base)

            print(f"Kernels dowloaded")
            return filename


    def get_convolve(self):
        
        obj_dir = self.get_kernels()

        if self.hdul is not None:
            hdu_inp = self.hdul
            data_inp = hdu_inp.data
            header_inp = hdu_inp.header
            #wcs_inp = WCS(header_inp)
            
        else:
            path_inp = self.path+'/'+str(self.name)+'/images/'+str(self.inp_surveys)+'.fits'
            hdu_inp = fits.open(path_inp)
            hdu_inp = hdu_inp[0]
            data_inp = hdu_inp.data
            header_inp = hdu_inp.header
            #wcs_inp = WCS(header_inp)

        path_ker = str(obj_dir)
        hdu_ker = fits.open(path_ker)
        hdu_ker = hdu_ker[0]
        data_ker = hdu_ker.data
        header_ker = hdu_ker.header

        self.ker_survey = self.ker_survey.split('_')[0]
        self.inp_surveys = self.inp_surveys.split('_')[0]

        pxs_ker = pxscale(header_ker)
        pxs_inp = pxscale(header_inp)
        if pxs_ker != pxs_inp:
                ratio = pxs_ker/pxs_inp
                size = ratio * data_ker.shape[0]
                if round(size) % 2 == 0:
                    size += 1
                    ratio = size / data_ker.shape[0]
        else:
                ratio = 1.
                #res = float(survey_resolution(self.ker_survey))

        newkernel = zoom(data_ker, ratio) / ratio**2

        astropy_conv = convolve_fft(data_inp, newkernel, nan_treatment = 'interpolate', normalize_kernel=True, preserve_nan=True)

        return astropy_conv
    

#gs = convolveimg('SDSS_g','unWISE_W2','/home/polo/Escritorio/Works/Doctorado/Code/SFHmergers/images','SIT45')
#gs.get_convolve()

