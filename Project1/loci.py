#! /usr/bin/env python
from __future__ import print_function, division
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
## pre-defined numbers
FWHM = 8.58 # pixels, measrued from unsaturated PSF
dr = 2.0 * FWHM
NA = 300
g = 1.0
pixelScale = 0.009942 #arcsec/pixel
rmin = 0.2/pixelScale
rmax = 3.0/pixelScale # the range of annuli, inner radius: 0.3 arcsec, outer radius: 3.0
rlist = np.linspace(rmin, rmax, np.floor((rmax - rmin)/dr))
deltaR = np.sqrt(np.pi*g*NA/4) * FWHM


def getPSFSet(angle0, r, angleList, minDist = 1.0 * FWHM):
    """
    get the proper PSFSet that
    input angle in degree
    r minDist in pixel
    """
    dis = np.abs(angle0 - angleList) * np.pi/180. * r
    return np.where(dis > minDist)[0]

def lociCoeff(image0, PSFCube):
    """
    calculate the coefficient for each individual subset
    """
    A = []
    b = []
    for i in range(PSFCube.shape[0]):
        A.append([(PSFCube[i] * PSFCube[j]).sum() for j in range(PSFCube.shape[0])])
        b.append((PSFCube[i] * image0).sum())
        
    A = np.mat(A)
    b = np.mat(b)
    return np.array((np.linalg.pinv(A.T) * b.T).flat) # pinv function in numpy linear algebra module (np.linalg) use moore-penrose method to calculate the pseudo-inverse of the matrix to avoid singular matrix
    
def loci(imageCube, rotateAngle):
    nImages = imageCube.shape[0] # number of images
    imageSize = imageCube.shape[1:] # image size in 1 and 2 dimension
    xx, yy = np.meshgrid(np.arange(imageSize[0]), np.arange(imageSize[1]))
    radius = np.sqrt((xx - imageSize[0]/2.)**2 + (yy - imageSize[1]/2.)**2)
    posAngle = np.arccos((xx - imageSize[0]/2.)/radius) * 180/np.pi
    posAngle[np.where(yy - imageSize[1]/2. < 0)] = posAngle[np.where(yy - imageSize[1]/2. < 0)] + 180
    subtractedCube = imageCube.copy()
    for i in range(nImages):
        image0 = imageCube[i, :, :]
        print('Number {0} image starts to be processed...'.format(i))
        for r_i in rlist:
            phiList = np.linspace(0, 360, np.round(360/(180/np.pi/(g/2. + 2*r_i/FWHM * np.sqrt(g/np.pi/NA)))))
            PSFList = getPSFSet(rotateAngle[i], r_i, rotateAngle)

            for phiIndex, phi_i in enumerate(phiList[:-1]):
                dim1, dim2 = np.where((radius>=r_i) & (radius <= r_i + deltaR) & (posAngle >= phi_i) & (posAngle <= phiList[phiIndex + 1]))
                coeff = lociCoeff(image0[dim1, dim2], imageCube[PSFList,:,:][:, dim1, dim2])

                for coeffIndex, coeff_i in enumerate(coeff):
                    subDim1, subDim2 = np.where((radius>=r_i) & (radius <= r_i + dr) & (posAngle >= phi_i) & (posAngle <= phiList[phiIndex + 1]))
                    subtractedCube[i, subDim1, subDim2] = subtractedCube[i, subDim1, subDim2] - imageCube[PSFList[coeffIndex], subDim1, subDim2] * coeff_i

    return subtractedCube

        
if __name__ == '__main__':
    dataCubeFN = 'center_im.fits' # file name for image cube
    imCube = fits.getdata(dataCubeFN, 'primary') # use getdata function in astropy.io module to obtain data from fits file
    rotAngleFN = 'rotnth.fits'
    rotAngle = fits.getdata(rotAngleFN, 'primary')
    outFN = 'loci_test'
    lociImageCube = loci(imCube, rotAngle)
    fits.writeto(outFN + '_cube.fits', lociImageCube, clobber = True)
