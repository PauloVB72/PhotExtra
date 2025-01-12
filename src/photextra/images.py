# Part of this image.py is based on https://github.com/temuller/hostphot
import os
from pathlib import Path
import logging

import time
from pyvo.dal.exceptions import DALFormatError

import glob
import copy
import shutil
import zipfile
import tarfile
import requests  # for several surveys
import numpy as np
import pandas as pd
from scipy.ndimage import rotate
# for VISTA
import re
import urllib
import montage_wrapper as montage
from astropy.io import fits
from astropy.table import Table
from astropy import wcs, units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import pyvo  # 2MASS
from pyvo.dal import sia  # DES, SkyMapper, S-Plus
from astroquery.sdss import SDSS
from astroquery.ukidss import Ukidss
from astroquery.skyview import SkyView  # other surveys
from astroquery.mast import Observations  # for GALEX
from astroquery.esa.hubble import ESAHubble  # HST
import warnings
from astropy.utils.exceptions import AstropyWarning

esahubble = ESAHubble()

from reproject import reproject_interp

from .params import parameters
from .utils import folder_exists,directory,setup_directories
from .utils import survey_pixel_scale
from .utils import check_filters
from .utils import bkg_sub
from .utils import header_changes


#repository_path = Path(Photsfh.__path__[0])



# 

class GetImages():
    valid_versions = {
        'SDSS': ['dr12', 'dr13', 'dr14', 'dr15', 'dr16', 'dr17'],
        'GALEX': ['AIS', 'MIS', 'DIS', 'NGS', 'GII'],
        'WISE': ['allwise', 'neo1', 'neo2', 'neo3', 'neo4', 'neo5', 'neo6', 'neo7'],
    }
    def __init__(self, name, ra,dec, size, surveys_init,versions = None, path = None, bkg_subtraction = True):
        """
        Initialize the GetImages class.

        Parameters:
        name (str): The name of the object.
        ra (float): The right ascension of the object.
        dec (float): The declination of the object.
        size (float): The size of the image.
        surveys_init (list): The list of surveys to download.
        versions (dict): The dictionary of versions for each survey. Defaults to None.
        path (str): The path to the data. Defaults to None.
        bkg_subtraction (bool): Whether to perform background subtraction. Defaults to True.
        """
        self.name = name
        self.position = (ra,dec)
        self.size = size
        self.versions = versions
        self.path = path
        self.bkg_subtraction = bkg_subtraction
        self.surveys = self.initialize_surveys(surveys_init, self.position, size)


    def initialize_surveys(self, surveys_init, position, size):
        try:
            parameters(surveys_init, position, size).check_validity()
            return parameters(surveys_init, position, size).survey_values()
        except Exception as e:
            warnings.warn(f'Error: {e}')
            return {}


    def download(self):
        """
        Download the images.

        Returns:
        str: The message indicating whether the download was successful.
        """
        gs = GetSurveys()
        for srv, filters in self.surveys.items():
            try:
                if self.versions is not None:
                    for survey, ver in self.versions.items():
                        if survey == srv:

                            gs.dowload_img(self.name, self.position[0], self.position[1], self.size, survey=srv, filters=filters, version=ver, path=self.path)
                        else:
                            gs.dowload_img(self.name, self.position[0], self.position[1], self.size, survey=srv, filters=filters, version=None, path=self.path)

            except Exception as e:
                warnings.warn(f'Error: {e}')
                return "Download failed"
        return "Download completed in folder."
    

class GetSurveys():
    def __init__(self,versions = None,bkg_subtraction = True):
        """
        Initialize the GetSurveys class.

        Parameters:
        versions (dict): The dictionary of versions for each survey. Defaults to None.
        bkg_subtraction (bool): Whether to perform background subtraction. Defaults to True.
        """
        self.versions = versions
        self.bkg_subtraction = bkg_subtraction

    def getimg_PS1(self, ra ,dec, size =3, filters = None):
        """
        Download the PS1 images.

        Parameters:
        ra (float): The right ascension of the object.
        dec (float): The declination of the object.
        size (float): The size of the image. Defaults to 3.
        filters (list): The list of filters to download. Defaults to None.

        Returns:
        list: The list of HDU ```python
        objects downloaded.
        """
       
        survey = "PS1"

        if filters is None:
            filters = 'grizy'

        pixel_scale = survey_pixel_scale(survey)      

        if isinstance(size, (float, int)):
            
            size_arcsec = (size * u.arcmin).to(u.arcsec).value
        else:
            size_arcsec = size.to(u.arcsec).value
        size_pixels = int(size_arcsec / pixel_scale)

        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        url = (
            f"{service}?ra={ra}&dec={dec}&size={size_pixels}&format=fits&"
            f"filters={filters}"
        )

        table = Table.read(url, format="ascii")    
        url = (
            "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
            f"ra={ra}&dec={dec}&size={size_pixels}&format=fits"
        )

        # sort filters from blue to red
        flist = ["grizy".find(x) for x in table["filter"]]
        table = table[np.argsort(flist)]

        base_url = url + "&red="
        url_list = []
        for filename in table["filename"]:
            url_list.append(base_url + filename)

        hdu_list = []
        for url, filt in zip(url_list, filters):
            hdu = fits.open(url)
            hdu_list.append(hdu)
        return hdu_list

    def getimg_SDSS(self, ra ,dec, size =3, filters = None, version=None):

        """
        Download the SDSS images.

        Parameters:
        ra (float): The right ascension of the object.
        dec (float): The declination of the object.
        size (float): The size of the image. Defaults to 3.
        filters (list): The list of filters to download. Defaults to None.
        version (str): The version of the PS1 data. Defaults to None.
        Returns:
        list: The list of HDU ```python
        objects downloaded.
        """
        survey = "SDSS"

        # check data release version
        if version is None:
            version = "dr17"

        versions = [f"dr{i}" for i in range(12, 17 + 1)]
        assert (
            version in versions
        ), f"The given version ({version}) is not a valid data release: {versions}"
        dr = int(version.replace("dr", ""))

        if isinstance(size, (float, int)):
            size_arcsec = (size * u.arcmin).to(u.arcsec)
        else:
            size_arcsec = size.to(u.arcsec)

        pixel_scale = survey_pixel_scale(survey)
        size_pixels = int(size_arcsec.value / pixel_scale)

        coords = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
        for radius in np.arange(1, 60, 1):
            ids = SDSS.query_region(
                coords, radius=radius * u.arcsec, data_release=dr
            )
            if ids is not None:
                if 'ra' in ids.colnames and 'dec' in ids.colnames:
                    break
                else:
                    ids = None

        if ids is None:
            return None

        # get the pointing closest to the given coordinates
        coords_imgs = SkyCoord(
            ra=ids["ra"].value,
            dec=ids["dec"].value,
            unit=(u.degree, u.degree),
            frame="icrs",
        )
        separation = coords.separation(coords_imgs).value

        pointing_id = np.argmin(separation)
        ids2remove = list(np.arange(len(separation)))
        del ids2remove[pointing_id]
        ids.remove_rows(ids2remove)

        # download images here
        hdu_list = SDSS.get_images(matches=ids, band=filters, data_release=dr)

        # SDSS images are large so need to be trimmed
        hdul_list = list()
    
        for hdul in hdu_list:
            
            hdulnew = montage.reproject_hdu(hdul[0], north_aligned=True)
            hdul_list.append(hdulnew)
        images_list = list()
        for hdu in hdul_list:
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", AstropyWarning)
                img_wcs = wcs.WCS(hdu.header)
            

            trimmed_data = Cutout2D(hdu.data, coords, size_pixels, img_wcs)
            hdu.data = trimmed_data.data
            hdu.header.update(trimmed_data.wcs.to_header())

        return hdul_list

    def getimg_GALEX(self, ra ,dec, size =3, filters = None, version=None):
        """
        Get GALEX images from the GALEX archive.

        Parameters:
        ra (float): Right ascension of the pointing in degrees.
        dec (float): Declination of the pointing in degrees.
        size (int): Size of the image in arcseconds.
        filters (list): List of filters to download. Default is None, which means all filters.
        version (str): Version of the GALEX data to download. Default is None, which
        means the latest version.
        Returns:
        list: The list of HDU ```python
        objects downloaded.
        """

        survey = "GALEX"

        if isinstance(size, (float, int)):
            size_arcsec = (size * u.arcmin).to(u.arcsec)
        else:
            size_arcsec = size.to(u.arcsec)

        pixel_scale = survey_pixel_scale(survey)
        size_pixels = int(size_arcsec.value / pixel_scale)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            coords = SkyCoord(
                ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs"
            )

        obs_table = Observations.query_criteria(
            coordinates=coords, radius=size_arcsec, obs_collection=["GALEX"]
        )
        obs_df_ = obs_table.to_pandas()

        if version is None:
            # starts from the survey with the highest images first
            obs_df_ = obs_df_.sort_values("t_exptime", ascending=True) # sort by exposure time, last modify: False
            projects = obs_df_.project.unique()
        else:
            # only use the survey requested by the user
            projects = [version]

        hdu_dict = {"NUV": None, "FUV": None}
        for project in projects:
            obs_df = obs_df_[obs_df_.project == project]

            for filt in filters:
                if hdu_dict[filt] is not None:
                    # if the image was already found, skip this filter
                    continue

                # get only "intensity" images
                filt_extension = {
                    "NUV": "-nd-int.fits.gz",
                    "FUV": "-fd-int.fits.gz",
                }
                # get unique image sectors
                files = []
                for file in obs_df.dataURL.values:
                    sector_info = file.split("-")[:-2]
                    file = (
                        "-".join(string for string in sector_info)
                        + filt_extension[filt]
                    )
                    if file not in files:
                        files.append(file)

                # download the FITS images
                hdu_list = []
                for file in files:
                    try:
                        hdu = fits.open(file)
                        hdu_list.append(hdu)
                    except:
                        pass

                # calculate the separation of the galaxy to the image center
                separations = []
                for hdu in hdu_list:
                    ra_img = float(hdu[0].header["RA_CENT"])
                    dec_img = float(hdu[0].header["DEC_CENT"])
                    coords_img = SkyCoord(
                        ra=ra_img,
                        dec=dec_img,
                        unit=(u.degree, u.degree),
                        frame="icrs",
                    )
                    separation = coords.separation(coords_img).value
                    separations.append(separation)

                # get the image with the galaxy closest to the center
                if len(separations) != 0:
                    id_file = np.argmin(separations)
                    hdu = hdu_list[id_file]
                else:
                    hdu = None

                hdu_dict[filt] = hdu

        hdu_list = []
        for filt in filters:
            hdu = hdu_dict[filt]
            if hdu is None:
                print('No image in this filter')
                continue  # no image in this filter
            
            # trim data to requested size
            img_wcs = wcs.WCS(hdu[0].header)
            pos = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)

            trimmed_data = Cutout2D(hdu[0].data, pos, size_pixels, img_wcs)
            hdu[0].data = trimmed_data.data
            hdu[0].header.update(trimmed_data.wcs.to_header())
            hdu_list.append(hdu)

        return hdu_list


    def getimg_WISE(self, ra ,dec, size =3, filters = None, version=None):

        """Downloads a set of GALEX fits images for a given set
        of coordinates and filters.

        Parameters
        ----------
        ra: str or float
            Right ascension in degrees.
        dec: str or float
            Declination in degrees.
        size: float or ~astropy.units.Quantity, default ``3``
            Image size. If a float is given, the units are assumed to be arcmin.
        filters: str, default ``None``
            Filters to use. If ``None``, uses ``FUV, NUV``.
        version: str, default ``None``
            Version of GALEX images. Either Deep (``DIS``), Medium (``MIS``) or
            All-Sky Imaging Survey (``AIS``), or Nearby Galaxy Survey (``NGS``) or
            Guest Investigator Survey (``GII``). If ``None``, take the image with the
            longest exposure time.

        Return
        ------
        hdu_list: list
            List with fits images for the given filters.
            ``None`` is returned if no image is found.
        """
        survey = "WISE"

        if isinstance(size, (float, int)):
            size_arcsec = (size * u.arcmin).to(u.arcsec)
        else:
            size_arcsec = size.to(u.arcsec)

        pixel_scale = survey_pixel_scale(survey)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            coords = SkyCoord(
                ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs"
            )
            # just to get the original image used
            skyview_fits = SkyView.get_images(
                position=coords,
                coordinates="icrs",
                pixels="100",
                survey="WISE 3.4",
            )
            header = skyview_fits[0][0].header
            image_line = header["HISTORY"][-3]
            used_image = image_line.split("/")[-1].split("-")[0]
            coadd_id = used_image
            coadd_id1 = coadd_id[:4]
            coadd_id2 = coadd_id1[:2]

            # for more info: https://irsa.ipac.caltech.edu/ibe/docs/wise/allwise/p3am_cdd/#sample_code
            base_url = (
                "http://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/"
            )
            coadd_url = os.path.join(coadd_id2, coadd_id1, coadd_id)
            params_url = f"center={ra},{dec}&size={size_arcsec.value}arcsec&gzip=0"  # center and size of the image

            hdu_list = []
            for filt in filters:
                i = filt[-1]
                band_url = f"{coadd_id}-w{i}-int-3.fits"
                url = os.path.join(
                    base_url, coadd_url, band_url + "?" + params_url
                )
                hdu = fits.open(url)
                hdu_list.append(hdu)

        for dt in hdu_list:
            try:
                dt[0].data
                if self.bkg_subtraction == True:
                    dt[0].data = bkg_sub(dt[0].data,survey=survey)
            except:
                pass       
        return hdu_list

   
    def getimg_unWISE(self, ra ,dec, size =3, filters = None, version="allwise"):
        """
        Get unWISE images from the IRSA website.
        Parameters
        ----------
        ra : float
        Right ascension of the target in degrees.
        dec : float
        Declination of the target in degrees.
        size : float
        Size of the image in arcseconds.
        filters : list
        List of filters to retrieve. Default is None, which means all filters.
        version : str
        Version of the unWISE catalog. Default is "allwise".
        Returns
        -------
        hdu_list : list

        List of HDU objects downloaded.
        """

        survey = "unWISE"

    
        if version is None:
            version = "allwise"
        else:
            # check validity of the version used
            neo_versions = [f"neo{i}" for i in range(1, 8)]
            all_versions = ["allwise"] + neo_versions
            assert (
                version in all_versions
            ), f"Not a valid version ({version}): {all_versions}"

        if "neo" in version:
            # these only have W1 and W2 data
            if "W3" in filters:
                filters.remove("W3")
            if "W4" in filters:
                filters.remove("W4")

        if isinstance(size, (float, int)):
            size_arcsec = (size * u.arcmin).to(u.arcsec)
        else:
            size_arcsec = size.to(u.arcsec)

        pixel_scale = survey_pixel_scale(survey)
        size_pixels = int(size_arcsec.value / pixel_scale)
        assert size_pixels <= 1024, "Maximum cutout size for unWISE is 1024 pixels"

        bands = "".join(filt[-1] for filt in filters)  # e.g. 1234

        # for more info: http://unwise.me/imgsearch/
        base_url = "http://unwise.me/cutout_fits?"
        params_url = (
            f"version={version}&ra={ra}&dec={dec}&size={size_pixels}&bands={bands}"
        )
        master_url = base_url + params_url

        response = requests.get(master_url, stream=True)
        target_file = Path(f"unWISE_images_{ra}_{dec}.tar.gz")  # current directory
        if response.status_code == 200:
            with open(target_file, "wb") as f:
                f.write(response.raw.read())

        hdu_list = []
        with tarfile.open(target_file) as tar_file:
            files_list = tar_file.getnames()
            for fits_file in files_list:
                for filt in filters:
                    if f"{filt.lower()}-img-m.fits" in fits_file:
                        tar_file.extract(fits_file, ".")
                        hdu = fits.open(fits_file)
                        hdu_list.append(hdu)
                        os.remove(fits_file)

        # sometimes, the file is not downloaded, so let's check if it exists
        if os.path.isfile(target_file) is True:
            os.remove(target_file)

        return hdu_list


    def getimg_2MASS(self, ra ,dec, size =3, filters = None):
        """
        Get 2MASS images from the GALEX archive.

        Parameters:
        ra (float): Right ascension of the pointing in degrees.
        dec (float): Declination of the pointing in degrees.
        size (int): Size of the image in arcseconds.
        filters (list): List of filters to download. Default is None, which means all filters.

        Returns:
        list: The list of HDU ```python
        objects downloaded.
        """
        
        survey = "2MASS"

        def retry_regsearch(max_retries=3, delay=5):
            for attempt in range(max_retries):
                try:
                    return pyvo.regsearch(servicetype="image", keywords=["2mass"])
                except DALFormatError as e:
                    print(f"Error en intento {attempt + 1}: {e}")
                    time.sleep(delay * (2 ** attempt))  # Backoff exponencial
            raise Exception("No se pudo conectar al servicio después de varios intentos")

        if isinstance(size, (float, int)):
            size_degree = (size * u.arcmin).to(u.degree)
        else:
            size_degree = size.to(u.degree)
        size_degree = size_degree.value

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            coords = SkyCoord(
                ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs"
            )
            twomass_services = retry_regsearch()

         #   twomass_services = pyvo.regsearch(
          #      servicetype="image", keywords=["2mass"]
          #
          #   )


            table = twomass_services[0].search(pos=coords, size=size_degree)
            twomass_df = table.to_table().to_pandas()
            twomass_df = twomass_df[twomass_df.format == "image/fits"]

            # for more info: https://irsa.ipac.caltech.edu/ibe/docs/twomass/allsky/allsky/#main
            base_url = (
                "https://irsa.ipac.caltech.edu/ibe/data/twomass/allsky/allsky"
            )

            hdu_list = []
            for filt in filters:
                if filt == "Ks":
                    filt = "K"
                band_df = twomass_df[twomass_df.band == filt]
                if len(band_df) == 0:
                    # no data for this band:
                    hdu_list.append(None)
                    continue

                # pick the largest images which is also the most "square" one
                # this works better than picking the image closest to the given coordinates
                # don't know why
                sizes = []
                tmp_hdu_list = []
                for i in range(len(band_df)):
                    fname = band_df.download.values[i].split("=")[-1]
                    hemisphere = band_df.hem.values[i]
                    ordate = band_df.date.values[i]
                    scanno = band_df.scan.values[i]
                    # add leading zeros for scanno bellow 100
                    n_zeros = 3 - len(str(scanno))
                    scanno = n_zeros * "0" + str(scanno)

                    tile_url = os.path.join(f"{ordate}{hemisphere}", f"s{scanno}")
                    fits_url = os.path.join("image", f"{fname}.gz")
                    params_url = f"center={ra},{dec}&size={size_degree}degree&gzip=0"  # center and size of the image

                    url = os.path.join(
                        base_url, tile_url, fits_url + "?" + params_url
                    )
                    try:
                        hdu = fits.open(url)
                        ny, nx = hdu[0].data.shape
                        sizes.append(nx * ny)
                        tmp_hdu_list.append(hdu)
                    except:
                        # some images might give 500 Internal Server Error
                        # because the cutout does not overlap the image
                        pass

                if len(tmp_hdu_list)==0:
                    hdu_list.append(None)
                else:
                    # pick largest image, which usually is the best
                    i = np.argmax(sizes)
                    hdu_list.append(tmp_hdu_list[i])
        
        for dt in hdu_list:
            try:
                dt[0].data
                if self.bkg_subtraction == True:
                    dt[0].data = bkg_sub(dt[0].data,survey=survey)
            except:
                pass

        return hdu_list

    def getimg_LegacySurvey(self, ra ,dec, size =3, filters = None,version=None):
  
        survey = "LegacySurvey"

        pixel_scale = survey_pixel_scale(survey)
        if isinstance(size, (float, int)):
            size_arcsec = (size * u.arcmin).to(u.arcsec).value
        else:
            size_arcsec = size.to(u.arcsec).value
        size_pixels = int(size_arcsec / pixel_scale)

        if version is None:
            version = "dr10"  # latest data release

        base_url = "https://www.legacysurvey.org/viewer/fits-cutout?"
        params = f"ra={ra}&dec={dec}&layer=ls-{version}&pixscale={pixel_scale}&bands={filters}&size={size_pixels}&invvar"
        url = base_url + params

        master_hdu = fits.open(url)
        master_header = master_hdu[0].header
        master_header_invvar = master_hdu[1].header

        hdu_list = []
        for i, filt in enumerate(filters):
            data = master_hdu[0].data[i]
            header = master_header.copy()
            header.append(('BAND', filt, ' Band - added by HostPhot'), end=True)
            hdu = fits.PrimaryHDU(data=data, header=header)
            
            header_invvar = master_header_invvar.copy()
            header_invvar.append(('BAND', filt, ' Band - added by HostPhot'), end=True)
            data_invvar = master_hdu[1].data[i]
            hdu_invvar = fits.ImageHDU(data=data_invvar, 
                                    header=header_invvar)

            hdu_list.append(fits.HDUList([hdu, hdu_invvar]))

        return hdu_list
    
    def getimg_VISTA(self, ra ,dec, size =3, filters = None,version = None):

        survey = "VISTA"
        

        if not isinstance(size, (float, int)):
            size = size.to(u.arcmin).value

        if version is None:
            version = "VHS"
        # These are final data releases - VIDEO has a DR6 and VIKING a DR5, but not 
        # if there is any difference:
        # VHSDR6: https://b2find.eudat.eu/dataset/0b10d3a0-1cfe-5e67-8a5c-0949db9d19cb
        # VIDEODR5: https://www.eso.org/sci/publications/announcements/sciann17491.html
        # VIKINGDR4: https://www.eso.org/sci/publications/announcements/sciann17289.html
        database_dict = {
            "VHS": "VHSDR6",
            "VIDEO": "VIDEODR6",
            "VIKING": "VIKINGDR5",
        }
        valid_surveys = list(database_dict.keys())
        assert (
            version in valid_surveys
        ), f"Not a valid VISTA survey: choose from {valid_surveys}"
        database = database_dict[version]

        base_url = "http://horus.roe.ac.uk:8080/vdfs/GetImage?archive=VSA&"
        survey_dict = {
            "database": database,
            "ra": ra,
            "dec": dec,
            "sys": "J",
            "filterID": "all",
            "size": size,  # in arcmin
            "obsType": "object",
            "frameType": "tilestack",
        }
        survey_url = "&".join([f"{key}={val}" for key, val in survey_dict.items()])

        url = base_url + survey_url
        results = urllib.request.urlopen(url).read()
        links = re.findall('href="(http://.*?)"', results.decode("utf-8"))

        # find url for each filter (None if not found)
        urls_dict = {filt: None for filt in filters}
        for filt in filters:
            for link in links:
                url = link.replace("getImage", "getFImage", 1)
                if f"band={filt}" in url:
                    urls_dict[filt] = url
                    break

        hdu_list = []
        for filt, url in urls_dict.items():
            if url is not None:
                hdu = fits.open(url)

                hdu[0].data = hdu[1].data
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", AstropyWarning)
                    img_wcs = wcs.WCS(hdu[1].header)
                    hdu[0].header.update(img_wcs.to_header())
                    
                    # add some keywords to the PHU
                    hdu[0].header['EXPTIME'] = hdu[1].header['EXPTIME']
                    hdu[0].header['MAGZRR'] = hdu[1].header['MAGZRR']

                    # calculate effective ZP (considering atmospheric extinction)
                    # calculate extinction first
                    airmass = (hdu[1].header['HIERARCH ESO TEL AIRM START'] + hdu[1].header['HIERARCH ESO TEL AIRM END'])/2
                    ext_coeff = hdu[1].header['EXTINCT']
                    extinction = ext_coeff*(airmass - 1)
                    # calculate effective ZP
                    zp = hdu[1].header['MAGZPT']
                    hdu[0].header['MAGZP'] = zp - extinction

                hdu_list.append(hdu)
            else:
                hdu_list.append(None)

        for dt in hdu_list:
            try:
                dt[0].data
                dt[0].data = bkg_sub(dt[0].data,survey=survey)
            except:
                pass
        return hdu_list


    def getimg_SPLUS(self, ra, dec, size, filters=None):

        survey = "SPLUS"
    
        splus_url = "https://datalab.noirlab.edu/sia/splus_dr1"

        if isinstance(size, (float, int)):
            fov = (size * u.arcmin).to(u.degree).value
        else:
            fov = size.to(u.degree).value

        svc = sia.SIAService(splus_url)
        imgs_table = svc.search(
            (ra, dec), (fov / np.cos(dec * np.pi / 180), fov), verbosity=2
        )
        if len(imgs_table) == 0:
            print(("Warning: empty table returned for " f"ra={ra}, dec={dec}"))
            return None

        imgs_df = pd.DataFrame(imgs_table)
        imgs_df = imgs_df[imgs_df.access_format.str.endswith('fits')]
        if len(imgs_df) == 0:
            return None

        url_list = []
        for filt in filters:
            filt_df = imgs_df[imgs_df.obs_bandpass==filt]
            if len(filt_df) == 0:
                url_list.append(None)
            else:
                fits_url = filt_df.access_url.values[0]  # first image
                url_list.append(fits_url)

        #global repository_path

        if url_list is None:
            return None

        hdu_list = []
        for url in url_list:
            if url is None:
                hdu_list.append(None)
            else:
                hdu = fits.open(url)

                # add zeropoint
                # file from https://splus.cloud/documentation/dr2_3
                #zps_file = repository_path.joinpath('filters', 'SPLUS', 'iDR3_zps.cat')
                zps_df = pd.read_csv('iDR3_zps.cat', sep='\\s+')

                field = hdu[0].header['OBJECT'].replace('_', '-')
                field_df = zps_df[zps_df['#field']==field]
                img_filt = hdu[0].header['FILTER']
                zp = field_df[img_filt].values[0]
                hdu[0].header['MAGZP'] = zp

                # initial EXPTIME is normalised, so it doesn't help
                hdu[0].header['EXPTIME'] = hdu[0].header['TEXPOSED']

                hdu_list.append(hdu)

        return hdu_list
    

    def getimg_UKIDSS(self, ra, dec, size, filters=None):


        survey = "UKIDSS"
        
        database = 'UKIDSSDR11PLUS'
        # programme = 'LAS'  # ['UDS', 'GCS', 'GPS', 'DXS', 'LAS']
        
        survey_pixel_scale(survey)
        if isinstance(size, (float, int)):
            size = (size * u.arcmin)

        u_obj = Ukidss(database=database)
        pos = SkyCoord(ra, dec, unit="deg")
        urls = u_obj.get_image_list(pos, waveband='all', frame_type="stack", image_width=size)

        hdu_dict = {filt:None for filt in filters}
        for filt in filters:
            for url in urls:
                # pick the first FITS image only (is this correct?)
                if f'band={filt}' in url:
                    hdu = fits.open(url)

                    # update first extension with data and WCS
                    hdu[0].data = hdu[1].data
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", AstropyWarning)
                        img_wcs = wcs.WCS(hdu[1].header)
                        hdu[0].header.update(img_wcs.to_header())
                    # add some of the keywords to the PHU
                    hdu[0].header['MAGZRR'] = hdu[1].header['MAGZRR']
                    hdu[0].header['GAIN'] = hdu[1].header['GAIN']
                    hdu[0].header['READNOIS'] = hdu[1].header['READNOIS']

                    # calculate effective ZP (considering atmospheric extinction)
                    # calculate extinction first
                    airmass = (hdu[0].header['AMSTART'] + hdu[0].header['AMEND'])/2
                    ext_coeff = hdu[1].header['EXTINCT']
                    extinction = ext_coeff*(airmass - 1)
                    # calculate effective ZP
                    zp = hdu[1].header['MAGZPT']
                    hdu[0].header['MAGZP'] = zp - extinction

                    hdu_dict[filt] = hdu
                    break                

        hdu_list = list(hdu_dict.values())

        for dt in hdu_list:
            try:
                dt[0].data
                dt[0].data = bkg_sub(dt[0].data,survey=survey)
            except:
                pass

        return hdu_list
    
    def dowload_img(self,name, ra ,dec,size =3, survey='SDSS',filters = None,version = None, path = None, overwrite=True):
        """
        Download images from the specified survey.

        Parameters:
        name (str): The name of the object.
        ra (float): The right ascension of the object.
        dec (float): The declination of the object.
        size (float): The size of the image.
        survey (str): The survey to download from.
        filters (list): The list of filters to apply.
        version (str): The version of the survey. Defaults to None.
        path (str): The path to save the images. Defaults to None.
        """
        obj_dir = setup_directories(name,path=path)['main']



        if survey == 'PS1':
            hdu_list = self.getimg_PS1(ra,dec,size=size,filters=filters)
        elif survey == 'SDSS':
            hdu_list = self.getimg_SDSS(ra,dec,size=size,filters=filters,version = None)
        elif survey == 'GALEX':
            hdu_list = self.getimg_GALEX(ra,dec,size=size,filters=filters,version = None)
        elif survey == 'WISE':
            hdu_list = self.getimg_WISE(ra,dec,size=size,filters=filters)
        elif survey == 'unWISE':
            hdu_list = self.getimg_unWISE(ra,dec,size=size,filters=filters,version = 'allwise')
        elif survey == '2MASS':
            hdu_list = self.getimg_2MASS(ra,dec,size=size,filters=filters)
        elif survey == 'LegacySurvey':
            hdu_list = self.getimg_LegacySurvey(ra,dec,size=size,filters=filters)
        elif survey == 'VISTA':
            hdu_list = self.getimg_VISTA(ra,dec,size=size,filters=filters,version = None)
        elif survey == 'SPLUS':
            hdu_list = self.getimg_SPLUS(ra,dec,size=size,filters=filters)
        elif survey == 'UKIDSS':
            hdu_list = self.getimg_UKIDSS(ra,dec,size=size,filters=filters)

        if hdu_list:
            for hdu, filt in zip(hdu_list, filters):
                if hdu is None:
                    continue  # skip missing filter/image

                else:
                    outfile = os.path.join(obj_dir+'/images', f"{survey}_{filt}.fits")
                if survey=='SDSS':
                    size_img = hdu.data.shape
                else:
                    size_img = hdu[0].data.shape
                hdu = header_changes(hdu,ra,dec,size_img,survey)
                if overwrite is True or os.path.isfile(outfile) is False:
                    hdu.writeto(outfile, overwrite=overwrite)
                else:
                    continue
  
# funcion que diga que no hay imagenes para ese survey de datos en caso de none. 




#ra = 351.2577 
#dec = -0.00041
#ra_gal2 = 20.0108974646	
#dec_gal2 = 14.3617675139

#ra_gal3=123.30674	
#dec_gal3 = 24.60798
#size= 3
#name ='SIT45'
#surveys_ints = ['SDSS_r']
#gs = GetImages(name,ra,dec,size,surveys_ints)


#gs.download()

