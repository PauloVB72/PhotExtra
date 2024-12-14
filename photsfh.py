import os
import numpy as np
from astropy.io import fits
from utils import directory

from images import GetImages
from convolution import ConvolutionIMG
from resampling import ReprojectIMG

class PHOTsfh:

    def __init__(self ,name ,ra ,dec ,size , survey_filter,survey_ref, version=None, path=None,**kwargs):
        """ Photometry by SFH models"""
        self.name = name
        self.ra = ra
        self.dec = dec
        self.size = size
        self.survey_filter = survey_filter
        self.version = version
        self.path = path
        self.survey_ref = survey_ref

    def processing_img(self):

        data_cube = []

        gs = GetImages(self.name, self.ra, self.dec, self.size, self.survey_filter,versions = self.version, path = self.path)
        gs.download()

        path_inp = self.path+'images'

        for srv_inp in self.survey_filter:
            convolve_img = ConvolutionIMG(srv_inp,self.survey_ref, path_inp, self.name)
            convolved_image = convolve_img.get_convolve()
            #data_cube.append(convolved_image)

            resample_img = ReprojectIMG(srv_inp,self.survey_ref,path = None, data=convolved_image)
            resampled_image = resample_img.get_reproject()
            data_cube.append(resampled_image)

        # Convert data_cube to a 3D numpy array
        data_cube = np.array(data_cube)
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['NAME_OBJ'] = self.name
        primary_hdu.header['RA'] = self.ra
        primary_hdu.header['DEC'] = self.dec
        primary_hdu.header['FOV'] = self.size
        primary_hdu.header['MASK_TYPE'] = 'MAP'
        primary_hdu.header['CONV_TYPE'] = 'FFT'
        primary_hdu.header['REP_TYPE'] = 'EXACT'
        primary_hdu.header['SURVEYS'] = ', '.join(self.survey_filter)  # Lista de surveys

        # Save the data cube to a FITS file
        hdus = [primary_hdu]
        for i, srv_inp in enumerate(self.survey_filter):
            hdu = fits.ImageHDU(data=data_cube[i], name=f'SURVEY_{srv_inp}')
            hdus.append(hdu)

        output_file = os.path.join(output_path, 'data_cube.fits')
        hdu = fits.PrimaryHDU(data=data_cube)
        hdu.writeto(output_file, overwrite=True)

        print(f"Data cube saved to {output_file}")
        return data_cube
