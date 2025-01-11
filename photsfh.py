import os
import numpy as np
from astropy.io import fits
from utils import directory,setup_directories
from images import GetImages
from convolution import ConvolutionIMG
from resampling import ReprojectIMG
from maskstars import  ObjectsIMG
import matplotlib.pyplot as plt
from ploting import Plotting
from astropy.wcs import WCS


class PHOTsfh:

    def __init__(self ,name ,ra ,dec ,size , survey_filter,survey_ref, version=None, path=None, mask_stars = True, table = False, plot = True,
                 spike_threshold = 0.03, gradient_threshold = 0.05,
                 sep_threshold = 5, sep_deblend = 0.005, r = 3.8, **kwargs):
        """ Photometry"""
        self.name = name
        self.ra = ra
        self.dec = dec
        self.size = size
        self.survey_filter = survey_filter
        self.version = version
        self.path = path
        self.survey_ref = survey_ref
        self.mask_stars = mask_stars
        self.spike_threshold = spike_threshold
        self.gradient_threshold = gradient_threshold
        self.sep_threshold = sep_threshold
        self.r = r 
        self.sep_deblend = sep_deblend
        self.table = table
        self.plot = plot
    def processing_img(self):

        data_cube = []
        data_cube_HDU = []
        print('...............................................................................................')
        print('...............................................................................................')
        print('...............................................................................................')
        print('...............................................................................................')
        print(f"Downloading images for {self.name}")
        gs = GetImages(self.name, self.ra, self.dec, self.size, self.survey_filter,versions = self.version, path = self.path)
        gs.download()

        path_inp = self.path
        print('...............................................................................................')
        print('...............................................................................................')
        print(f"Images downloaded for {self.name}")
        survey_images = {srv: [] for srv in self.survey_filter}
        surv_reference = []


        
        for srv_inp in self.survey_filter:
            if srv_inp != self.survey_ref:
                srv_org = srv_inp.split('_')[0]
                if self.mask_stars:
                    print('...............................................................................................')
                    print('...............................................................................................')
                    print(f"Masking stars for {srv_inp}")
                    hdu= fits.open(path_inp+'/'+self.name+'/images/'+srv_inp+'.fits')[0]
                    masking = ObjectsIMG(hdu, self.ra, self.dec, self.size, spike_threshold = self.spike_threshold,
                                                  gradient_threshold = self.gradient_threshold, survey = srv_org,
                                                   sep_threshold = self.sep_threshold, sep_deblend = self.sep_deblend, r = self.r).masked()
                    print('...............................................................................................')
                    print('...............................................................................................')
                    print(f"Convolution for {srv_inp}")
                    convolve_img = ConvolutionIMG(srv_inp,self.survey_ref,self.name, path=path_inp, hdul=masking)
                    convolved_image = convolve_img.get_convolve()

                else:
                    print('...............................................................................................')
                    print('...............................................................................................')
                    print(f"Convolution for {srv_inp}")
                    convolve_img = ConvolutionIMG(srv_inp,self.survey_ref,self.name, path=path_inp)
                    convolved_image = convolve_img.get_convolve()

                hdu_conv = fits.open(path_inp+'/'+self.name+'/images/'+srv_inp+'.fits')[0].header
                #data_cube.append(convolved_image)
                print('...............................................................................................')
                print('...............................................................................................')
                print(f"Resampling for {srv_inp}")
                resample_img = ReprojectIMG(srv_inp,self.survey_ref,self.name,path = path_inp, data=(convolved_image,hdu_conv))
                resampled_image = resample_img.get_reproject()
                data_cube.append(resampled_image)
            else:
                print('...............................................................................................')
                print('...............................................................................................')
                print(f"Survey reference: {srv_inp}")
                hdu_fits_ref = fits.open(path_inp+'/'+self.name+'/images/'+srv_inp+'.fits')[0]
                masking = ObjectsIMG(hdu_fits_ref, self.ra, self.dec, self.size, spike_threshold = self.spike_threshold,
                                                  gradient_threshold = self.gradient_threshold, survey = srv_org,
                                                   sep_threshold = self.sep_threshold, sep_deblend = self.sep_deblend, r = self.r).masked()
                data_cube.append(masking.data)
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
        primary_hdu.header['MASK_STARS'] = str(self.mask_stars)
        primary_hdu.header['SURVEYS'] = ', '.join(self.survey_filter)  # Lista de surveys

        # Save the data cube to a FITS file

        hdus = [primary_hdu]

        surveys_sep = list({filt.split('_')[0] for filt in self.survey_filter})

        for i, srv_inp in enumerate(surveys_sep):
            data_surveys = []
            for j in self.survey_filter:
                if j.split('_')[0] == srv_inp:
                    hdu_fits = fits.open(path_inp+'/'+self.name+'/images/'+j+'.fits')[0]
                    data_surveys.append(hdu_fits.data)
      
            hdu = fits.ImageHDU(data=data_surveys, name=f'SURVEY_{srv_inp}',header=hdu_fits.header)
            hdu.header['SURVEY'] = srv_inp
            hdu.header['NIMAGES'] = len(data_surveys) 
            hdus.append(hdu)

        hdu_referece = fits.open(path_inp+'/'+self.name+'/images/'+self.survey_ref+'.fits')[0].header

        for i, srv_inp in enumerate(enumerate(self.survey_filter), start=1):  # start=1 para empezar desde CHAN1

            key = f'CHAN{i}'   # Generar la clave: CHAN1, CHAN2, ...
            hdu_referece[key] = srv_inp
        #hdu_referece.flush()
        last_hdu = fits.ImageHDU(data=data_cube, name='FINAL_PRODUCT',header=hdu_referece)
        hdus.append(last_hdu)

        output_file = os.path.join(self.path+'/'+self.name, 'data_cube.fits')
        fits.HDUList(hdus).writeto(output_file, overwrite=True)

        if self.plot:
            hdu = fits.open(output_file)
            data_cube = hdu[-1].data
            data_1 = []
            wcs_1 = []
            hdu_1 = []
            for i in range(1, len(hdu) - 1):
                    datas = hdu[i].data
                    wcss = WCS(hdu[i].header)
                    for j in range(len(datas)):    
                        data_1.append(datas[j])
                        wcs_1.append(wcss)
                        new_hdu = fits.ImageHDU(data=datas[j], header=hdu[i].header)
                        hdu_1.append(new_hdu)
                        
            data_2 = hdu[-1].data
            wcs_2 = WCS(hdu[-1].header)
            
            for k in range(len(data_cube)):
                path = self.path+'/'+self.name+'/img_'+str(k)
                stars, _ =  ObjectsIMG(hdu_1[k], self.ra, self.dec, self.size).objects_search()
                Plotting(data_1[k], data_2[k], hdu_1[k], wcs_2,stars, path).plot_pdf()
        else:

            pass


        print(f"Data cube saved to {output_file}")
        return data_cube

ra = 351.2577 
dec = -0.00041

ra_gal2 = 20.0108974646	
dec_gal2 = 14.3617675139

ra_gal3=123.30674	
dec_gal3 = 24.60798
size= 3

ra_gal4 = 	122.390684521
dec_gal4 = 36.9852665635


ra_m1,dec_m1 = 141.414344,	32.803871


ra_m3,dec_m3 = 166.837,	18.4329

name ='m3'
surveys_ints = ['SDSS_r','SDSS_g','SDSS_i','SDSS_u','SDSS_z','GALEX_FUV','GALEX_NUV','unWISE_W1',
                'unWISE_W2','unWISE_W3']

version = {
    'GALEX': 'DIS',
}

ph = PHOTsfh(name,ra_m3,dec_m3,size,surveys_ints,'unWISE_W3',path='/home/polo/Escritorio/Works/PHOTSFH_PRUEBAS',version=version)
ph.processing_img()
