from astroquery.sdss import SDSS
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u
Vizier.ROW_LIMIT = -1
import numpy as np
import matplotlib.pyplot as plt
from skimage.filters import threshold_otsu
from scipy.ndimage import label, generate_binary_structure, gaussian_filter
from skimage.morphology import remove_small_objects, disk, binary_dilation
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
import sep

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
    def __init__(self, hdu, ra, dec, size, spike_threshold, gradient_threshold):
        
        self.hdu = hdu
        self.ra = ra
        self.dec = dec
        self.size = size
        Vizier.ROW_LIMIT = -1
        self.spike_threshold = spike_threshold
        self.gradient_threshold = gradient_threshold
        
    @staticmethod
    def query_stars(ra, dec, radius=3):
        """
        Consulta estrellas de los catálogos Tycho-2 y Gaia en un área específica del cielo.
    
        :param ra: Ascensión recta del centro (en grados).
        :param dec: Declinación del centro (en grados).
        :param radius: Radio de búsqueda (astropy.units, por defecto 3 arcmin).
        :return: Dos tablas de astropy Table con los resultados de Tycho-2 y Gaia.
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
            # Convierte los resultados a una tabla Astropy para manejo sencillo
            #table = Table(results)
            return results
            
            # Guarda los resultados en un archivo CSV
           # table.write("galaxias_sdss.csv", format="csv", overwrite=True)
           # print("Consulta completada y resultados guardados en 'galaxias_sdss.csv'")
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

    def create_star_mask(self, image, method="otsu", intensity_threshold=None, spike_threshold=0.3, gradient_threshold=0.05,survey = 'SDSS'):
        # Accede directamente a self.objects_search
        objetos_img = self.objects_search()
        wcs = WCS(self.hdu[0].header)
        estrella = []
        for k in objetos_img[0]:
            obj_x, obj_y = wcs.wcs_world2pix(k[0], k[1], 1)
            estrella.append((int(obj_x),int(obj_y)))
        estrella_x ,estrella_y= zip(*estrella)

        
        mask = np.zeros_like(self.hdu[0].data, dtype=bool)

        if survey == 'SDSS':
            for i in range(len(estrella)):
                spike_mask = OBJECTS_IMG.find_and_mask_spikes(self.hdu[0].data, (estrella_x[i], estrella_y[i]),spike_threshold)
                mask |= spike_mask
            mask = remove_small_objects(mask, min_size=20)
            return mask
        else:
            bkg = sep.Background(self.hdu[0].data)
            background = bkg.back()  # Fondo estimado

            objects = sep.extract(self.hdu[0].data, 5, err=bkg.globalrms, deblend_cont=0.005)
        
    def masked(self):
        
        mask = self.create_star_mask(self.hdu[0].data,spike_threshold=self.spike_threshold, gradient_threshold=self.gradient_threshold)
        
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        bkg = Background2D(self.hdu[0].data, (50, 50), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        data_mk = np.copy(self.hdu[0].data)
        data_mk[mask] = bkg.background_rms_median
        return data_mk


    def objects_search(self):
    
        tab_gal = OBJECTS_IMG.gal_detected(self.ra,self.dec,size = self.size)
        tab_tycho, tab_gaia = OBJECTS_IMG.query_stars(self.ra,self.dec,radius = self.size)

        header = self.hdu[0].header
        data = self.hdu[0].data
        wcs = WCS(header)
        coord_pxtoworld = wcs.pixel_to_world([0,data.shape[0]], [0,data.shape[1]])
        xymin = coord_pxtoworld[0]
        xymax = coord_pxtoworld[1]
        print(xymin,xymax)
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
        
        #objetos = []
       # for star in range(len(filtered_stars)):     
        #    ra_star, dec_star = filtered_stars[star]  
         #   if ra_star < 100.:
          #      ra_star_ = ra_star + 360.
           # else:
        #        ra_star_ = ra_star
         #   for gal in range(len(filtered_galaxies)):
        #        
          #      ra_gal, dec_gal = filtered_galaxies[gal]
           #     
             #   if ra_gal< 100.:
            #        ra_gal_ = ra_gal + 360. 
             #   else:
              #      ra_gal_ = ra_gal
                    
              #  d1 = np.abs(ra_star_- ra_gal_)*3600
              #  d2 = np.abs(dec_star - dec_gal)*3600
              #  if d1<=15 and d2<=15:
               #     objetos.append((ra_star,dec_star))
        #rem = [filtered_stars.remove(x) for x in objetos]
        
        return filtered_stars_result