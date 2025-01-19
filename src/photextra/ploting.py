import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval, ImageNormalize


class Plotting:
    
    def __init__(self, data_1, data_2, wcs_1, wcs_2,stars, path):
        self.data_1 = data_1
        self.data_2 = data_2
        self.wcs_1 = wcs_1
        self.wcs_2 = wcs_2
        self.stars = stars
        self.path = path

    def plot_pdf(self, **kwargs):

        interval = ZScaleInterval()
        norm_1 = ImageNormalize(self.data_1, interval=interval)
        norm_2 = ImageNormalize(self.data_2, interval=interval)
        estrella = []
        wcs = WCS(self.wcs_1.header)
        for k in self.stars:

            obj_x, obj_y = wcs.wcs_world2pix(k[0], k[1], 1)

            estrella.append((int(obj_x),int(obj_y)))
        estrella_x ,estrella_y= zip(*estrella)
        # Graficar la imagen con los ejes en RA y DEC en grados
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4),subplot_kw={'projection': wcs})

        im1 = ax1.imshow(self.data_1, origin='lower', cmap='BuPu', norm=norm_1) 
        ax1.scatter(estrella_x, estrella_y, s=30, edgecolor='coral', facecolor='C0', lw=2)
        ax1.coords.grid(color='white', ls='dotted')
        ax1.coords[0].set_axislabel('RA (deg)')
        ax1.coords[1].set_axislabel('DEC (deg)')
        ax1.coords[0].set_major_formatter('d.ddd')
        ax1.coords[1].set_major_formatter('d.ddd')

        im2 = ax2.imshow(self.data_2, origin='lower', cmap='BuPu', norm=norm_2)
        ax2.coords.grid(color='white', ls='dotted')
        ax2.coords[0].set_axislabel('RA (deg)')
        ax2.coords[1].set_axislabel('DEC (deg)')
        ax2.coords[0].set_major_formatter('d.ddd')
        ax2.coords[1].set_major_formatter('d.ddd')

        fig.colorbar(im1, ax=ax1)
        fig.colorbar(im2, ax=ax2)
        plt.savefig(self.path+'.pdf', bbox_inches='tight')
        
        plt.close()
        
        return 