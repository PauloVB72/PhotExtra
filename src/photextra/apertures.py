

class Aperture_map:
    pass

class Aperture_ellipse:
    pass
class Aperture_segmentation:
    pass
class ApertureIMG(Aperture_map, Aperture_ellipse, Aperture_segmentation):

    def __init__(self, img, aperture = 'MAP' ):
        self.img = img
        self.aperture = aperture



    def __str__(self):
        return f"ApertureIMG(img={self.img}, x={self.x}, y={self.y}, r={self.r})"