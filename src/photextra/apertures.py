
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from photutils import PetrosianPhotometry, EllipticalAperture
from sep import extract
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground

class Aperture:
    def __init__(self, image, wcs, mask=None):
        """
        Initialize the Aperture class.

        Parameters:
        image (ndarray): The input image.
        wcs (WCS): The WCS object associated with the image.
        mask (ndarray): The mask for the image. Defaults to None.
        """
        self.image = image
        self.wcs = wcs
        self.mask = mask

    def petrosian_aperture(self, threshold=0.2):
        """
        Calculate the Petrosian aperture.

        Parameters:
        threshold (float): The threshold for the Petrosian radius. Defaults to 0.2.

        Returns:
        PetrosianPhotometry: The Petrosian photometry object.
        """
        petro = PetrosianPhotometry(threshold=threshold)
        return petro(self.image, self.mask)

    def elliptical_aperture(self, x, y, a, b, theta):
        """
        Calculate the elliptical aperture.

        Parameters:
        x (float): The x-coordinate of the center of the ellipse.
        y (float): The y-coordinate of the center of the ellipse.
        a (float): The semi-major axis of the ellipse.
        b (float): The semi-minor axis of the ellipse.
        theta (float): The position angle of the ellipse.

        Returns:
        EllipticalAperture: The elliptical aperture object.
        """
        ellipse = EllipticalAperture((x, y), a, b, theta=theta)
        return ellipse

    def sep_aperture(self, threshold=5, deblend_cont=0.005):
        """
        Calculate the SEP aperture.

        Parameters:
        threshold (float): The threshold for the SEP detection. Defaults to 5.
        deblend_cont (float): The deblending threshold for SEP. Defaults to 0.005.

        Returns:
        ndarray: The SEP aperture mask.
        """
        data = self.image.byteswap().newbyteorder()
        bkg = sep.Background(data)
        objects = sep.extract(data, threshold, err=bkg.globalrms, deblend_cont=deblend_cont)
        return objects

    def map_fluxes(self, aperture):
        """
        Calculate the fluxes within the given aperture.

        Parameters:
        aperture (ndarray): The aperture mask.

        Returns:
        ndarray: The fluxes within the aperture.
        """
        fluxes = np.sum(self.image * aperture)
        return fluxes

    def calculate_aperture(self, aperture_type, **kwargs):
        """
        Calculate the aperture based on the given type.

        Parameters:
        aperture_type (str): The type of aperture to calculate. Can be 'petrosian', 'ellipse', or 'sep'.
        **kwargs: Additional keyword arguments for the aperture calculation.

        Returns:
        ndarray: The aperture mask.
        """
        if aperture_type == 'petrosian':
            return self.petrosian_aperture(**kwargs)
        elif aperture_type == 'ellipse':
            return self.elliptical_aperture(**kwargs)
        elif aperture_type == 'sep':
            return self.sep_aperture(**kwargs)
        else:
            raise ValueError("Invalid aperture type. Must be 'petrosian', 'ellipse', or 'sep'.")

    def calculate_fluxes(self, aperture):
        """
        Calculate the fluxes within the given aperture.

        Parameters:
        aperture (ndarray): The aperture mask.

        Returns:
        ndarray: The fluxes within the aperture.
        """
        return self.map_fluxes(aperture)