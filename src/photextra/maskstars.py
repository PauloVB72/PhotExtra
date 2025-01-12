
from astroquery.sdss import SDSS
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
Vizier.ROW_LIMIT = -1
import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects, disk, binary_dilation
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
import sep
from astropy.wcs import WCS
import logging

def arcsec_to_deg(segundos):

  minutos = segundos / 60
  grados = minutos / 60
  return grados
    
def arcmin_to_deg(arcmin):

  grados = arcmin / 60
  return grados
    
def deg_to_arcmin(grados):

  minutos_arco = grados * 60
  return minutos_arco




class OBJECTS_IMG:
    def __init__(self, hdu, ra, dec, size, spike_threshold = 0.3, gradient_threshold = 0.05,
                survey = 'SDSS', sep_threshold = 5, sep_deblend = 0.005, r = 3.8):
        """
        Initialize the OBJECTS_IMG class.

        Parameters:
        hdu (fits.HDUList): The input image.
        ra (float): The right ascension of the object.
        dec (float): The declination of the object.
        size (float): The size of the image.
        spike_threshold (float): The threshold for spike detection. Defaults to 0.3.
        gradient_threshold (float): The threshold for gradient detection. Defaults to 0.05.
        survey (str): The survey name. Defaults to 'SDSS'.
        sep_threshold (float): The threshold for SEP detection. Defaults to 5.
        sep_deblend (float): The deblending threshold for SEP. Defaults to 0.005.
        r (float): The radius for SEP detection. Defaults to 3.8.
        """
        self.hdu = hdu
        self.ra = ra
        self.dec = dec
        self.size = size
        Vizier.ROW_LIMIT = -1
        self.spike_threshold = spike_threshold
        self.gradient_threshold = gradient_threshold
        self.survey = survey
        self.sep_threshold = sep_threshold
        self.r = r 
        self.sep_deblend = sep_deblend
        self.surveys_highresol = ['SDSS', 'LegacySurvey','SPLUS','PS1']
        3
    @staticmethod
    def query_stars(ra, dec, radius=3):
        """
        Query stars from the Tycho-2 and Gaia catalogs.

        Parameters:
        ra (float): The right ascension of the object.
        dec (float): The declination of the object.
        radius (float): The radius for the query. Defaults to 3.

        Returns:
        tuple: A tuple containing the Tycho-2 and Gaia catalogs.
        """
        radius = arcmin_to_deg(radius)*u.deg
        # Coordenadas del centro
        center = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    
        # Configurar Vizier para los catálogos Tycho-2 y Gaia
        Vizier.ROW_LIMIT = -1  # Sin límite de filas
        tycho2_catalog = "I/259/tyc2"  # Identificador del catálogo Tycho-2
        gaia_catalog = "I/345/gaia2"  # Identificador del catálogo Gaia DR2
    
        # Consultar Tycho-2s.open('/home/polo/Escritorio/Works/PHOTSFH_PRUEBAS/Prueba_1/images/SDSS_r.fits')
        tycho2_results = Vizier.query_region(center, radius=radius, catalog=tycho2_catalog)
        tycho2_table = tycho2_results[0] if tycho2_results else None
    
        # Consultar Gaia DR2
        gaia_results = Vizier.query_region(center, radius=radius, catalog=gaia_catalog)
        gaia_table = gaia_results[0] if gaia_results else None
    
        return tycho2_table, gaia_table
        
    @staticmethod
    def gal_detected(ra,dec,size=3):
        
        rad = arcmin_to_deg(size)
        ra_min = ra-rad
        ra_max = ra+rad
        dec_min = dec-rad
        dec_max = dec+rad
        query = f"""
        SELECT
            p.objid,
            p.ra,
            p.dec
        FROM PhotoObj AS p
        JOIN SpecObj AS s ON p.objid = s.bestobjid
        WHERE s.class = 'GALAXY'
            AND p.ra BETWEEN {ra_min} AND {ra_max}
            AND p.dec BETWEEN {dec_min} AND {dec_max}
        """
        
        # Ejecuta la consulta
        try:
            results = SDSS.query_sql(query)
            return results

        except Exception as e:
            print(f"Error al realizar la consulta: {e}")

        
    @staticmethod
    def find_and_mask_spikes(image, center, spike_threshold=0.3, gradient_threshold=0.05):
        """
        Identify and mask spikes radiating from the center of a star.
    
        Parameters:
            image (np.ndarray): Input image with stars and spikes.
            center (tuple): Coordinates of the star's center (x, y).
            spike_threshold (float): Intensity threshold to identify spikes.
            gradient_threshold (float): Minimum gradient change to detect spikes.
    
        Returns:
            np.ndarray: Binary mask with spikes masked.
        """
        mask = np.zeros_like(image, dtype=bool)
        x_center, y_center = center
        size = image.shape[0]
    
        # Check directions radiating from the center
        angles = np.linspace(0, 2 * np.pi, 720)  # Increase resolution for better spike detection
        for angle in angles:
            x, y = x_center, y_center
            prev_intensity = image[int(y), int(x)]
            while 0 <= x < size and 0 <= y < size:
                current_intensity = image[int(y), int(x)]
                gradient = current_intensity - prev_intensity
                if current_intensity > spike_threshold or gradient > gradient_threshold:
                    mask[int(y), int(x)] = True
                else:
                    break
                prev_intensity = current_intensity
                x += np.cos(angle)
                y += np.sin(angle)
    
        return mask

    @staticmethod
    def sep_aperture(hdu, deblend_cont = 0.005 ,threshold= 5 ):
        """
        Perform SEP aperture photometry.

        Parameters:
        hdu (fits.HDUList): The input image.
        deblend_cont (float): The deblending threshold for SEP. Defaults to 0.005.
        threshold (float): The threshold for SEP detection. Defaults to 5.

        Returns:
        np.ndarray: The SEP aperture photometry results.
        """
        wcs = WCS(hdu.header)
        data = hdu.data.byteswap().newbyteorder()
        bkg = sep.Background(data)
        objects = sep.extract(data, threshold, err=bkg.globalrms, deblend_cont=deblend_cont)
        
        return objects

    def create_star_mask(self, image, sep_threshold= 5, sep_deblend = 0.005 ,r = 3.8, spike_threshold=0.3, gradient_threshold=0.05, survey = 'SDSS'):
        """
        Create a mask for the stars in the image.

        Parameters:
        image (np.ndarray): The input image.
        sep_threshold (float): The threshold for SEP detection. Defaults to 5.
        sep_deblend (float): The deblending threshold for SEP. Defaults to 0.005.
        r (float): The radius for SEP detection. Defaults to 3.8.
        spike_threshold (float): The threshold for spike detection. Defaults to 0.3.
        gradient_threshold (float): The threshold for gradient detection. Defaults to 0.05.
        survey (str): The survey name. Defaults to 'SDSS'.

        Returns:
        np.ndarray: A binary mask indicating the positions of stars.
        """
        objetos_img = self.objects_search()
        wcs = WCS(self.hdu.header)
        estrella = []
        for k in objetos_img[0]:
            obj_x, obj_y = wcs.wcs_world2pix(k[0], k[1], 1)
            estrella.append((int(obj_x),int(obj_y)))
        estrella_x ,estrella_y= zip(*estrella)

        mask = np.zeros_like(self.hdu.data, dtype=bool)
        

        if survey in self.surveys_highresol:
            for i in range(len(estrella)):
                spike_mask = OBJECTS_IMG.find_and_mask_spikes(self.hdu.data, (estrella_x[i], estrella_y[i]),spike_threshold,gradient_threshold)
                mask |= spike_mask
                
            mask = remove_small_objects(mask, min_size=20)

        else:
            
            objects = OBJECTS_IMG.sep_aperture(self.hdu, deblend_cont = sep_deblend, threshold = sep_threshold)
            objs_coord = wcs.pixel_to_world(objects["x"], objects["y"])
            obj_cat = []

            for ra_star, dec_star in objetos_img[0]:
                
                star_coord = SkyCoord(ra_star, dec_star, unit="deg")
                
                separations = star_coord.separation(objs_coord)
                
                objs_filters = objs_coord[np.where(separations.arcsecond < 10)]
                if objs_filters :
                    obj_cat.append(np.where(separations.arcsecond < 10)[0][0])
                    
            sep.mask_ellipse(
                            mask,
                            objects["x"][obj_cat],
                            objects["y"][obj_cat],
                            objects["a"][obj_cat],
                            objects["b"][obj_cat],
                            objects["theta"][obj_cat],
                            r=r,
                        )
        return mask
    def log_star_detection(self, stars):
        """
        Log the details of detected stars.

        Parameters:
        stars (list): List of detected stars.
        """
        logging.info(f"Detected {len(stars)} stars.")
        for star in stars:
            logging.info(f"Star at RA: {star['ra']}, DEC: {star['dec']}, Magnitude: {star['magnitude']}")

        
    def masked(self):
        
        mask = self.create_star_mask(self.hdu.data,survey = self.survey, sep_threshold = self.sep_threshold, r = self.r, sep_deblend = self.sep_deblend,
                                     spike_threshold=self.spike_threshold, gradient_threshold=self.gradient_threshold)
        if self.survey in self.surveys_highresol:
            sigma_clip = SigmaClip(sigma=3.0)
            bkg_estimator = MedianBackground()
            bkg = Background2D(self.hdu.data, (50, 50), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
            data_mk = np.copy(self.hdu.data)
            background_mean = bkg.background_median
            background_rms = bkg.background_rms_median
            simulated_background = np.random.normal(loc=background_mean, scale=background_rms, size=self.hdu.data.shape)
            data_mk[mask] = simulated_background[mask]

            #data_mk[mask] = 0

        else:
            
            wcs = WCS(self.hdu.header)
            data = self.hdu.data.byteswap().newbyteorder()
            bkg = sep.Background(data)
            data_mk = np.copy(self.hdu.data)
            background_mean = bkg.back()  # Fondo medio
            background_rms = bkg.rms()  # RMS del fondo
            background_rms = np.maximum(background_rms, 0)
            simulated_background = np.random.normal(loc=background_mean, scale=background_rms, size=self.hdu.data.shape)
            data_mk[mask] = simulated_background[mask]
            #data_mk[mask] = 0

        return data_mk


    def objects_search(self):
    
        tab_gal = OBJECTS_IMG.gal_detected(self.ra,self.dec,size = self.size)
        tab_tycho, tab_gaia = OBJECTS_IMG.query_stars(self.ra,self.dec,radius = self.size)

        header = self.hdu.header
        data = self.hdu.data
        wcs = WCS(header)
        coord_pxtoworld = wcs.pixel_to_world([0,data.shape[0]], [0,data.shape[1]])
        xymin = coord_pxtoworld[0]
        xymax = coord_pxtoworld[1]
        filtered_galaxies = []
        filtered_stars = []
        
        for row in tab_gal.iterrows():
            ra, dec = row[1], row[2]

            if xymax.ra.value <= ra <= xymin.ra.value and xymin.dec.value <= dec <= xymax.dec.value:

                filtered_galaxies.append((ra, dec))
        filtered_galaxies.append((self.ra, self.dec))
        
        if tab_tycho is not None:
    
            for row in tab_tycho.iterrows():
                ra, dec = row[8], row[9]
                if xymax.ra.value <= ra <= xymin.ra.value and xymin.dec.value <= dec <= xymax.dec.value:
                   
                    filtered_stars.append((ra, dec))
        
        for row in tab_gaia.iterrows():
            ra, dec = row[0], row[2]
            
            if xymax.ra.value <= ra <= xymin.ra.value and xymin.dec.value <= dec <= xymax.dec.value:
                
                filtered_stars.append((ra, dec))
            
        galaxies_coords = SkyCoord([gal[0] for gal in filtered_galaxies], 
                                   [gal[1] for gal in filtered_galaxies], unit="deg")
        filtered_stars_result = []
    
        # Calcular separaciones angulares
        for ra_star, dec_star in filtered_stars:
            star_coord = SkyCoord(ra_star, dec_star, unit="deg")
            separations = star_coord.separation(galaxies_coords)
            
            
            if all(separations.arcsecond > 15):
                filtered_stars_result.append((ra_star, dec_star))
        
        return filtered_stars_result,filtered_galaxies
        

class ObjectsIMG(OBJECTS_IMG):

    def __init__(self, hdu, ra, dec, size, spike_threshold = 0.3, gradient_threshold = 0.05,
                survey = 'SDSS', sep_threshold = 5, sep_deblend = 0.005, r = 3.8):
        
        super().__init__(hdu, ra, dec, size, spike_threshold, gradient_threshold, survey, sep_threshold, sep_deblend, r)
        self.hdu = hdu

    
    def masked(self):
        newdata = super().masked()
        self.hdu.data = newdata
        return self.hdu



