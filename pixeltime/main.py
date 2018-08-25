from astropy.io import fits
import numpy as np


#defaultPath = ('/data/External/ISIMCV3_unzipped/' +
#               'NRCNRCA4-DARK-6020143529_1_1_40269_JW1_JLAB88_20160120T143723.248_20160120T145653.554/' +
#               'NRCNRCA4-DARK-60201435291_1_484_SE_2016-01-20T15h12m28.fits')


defaultPath = ('/Users/everettschlawin/Documents/jwst/test_data/cv3/full_frame_wlp8_05A4/' +
               'NRCN821WLP8FULL5-6012184421_1_481_SE_2016-01-12T20h26m57_I025.fits')

class get_refpix_tseries():
    """ A class to get the reference pixel time series """
    
    def __init__(self,path=defaultPath,inputHDU=None):
        self.nrefRows = 4
        self.nrefCols = 4
        self.nAmps = 4
        self.ampDirs = [1,-1,1,-1]
        self.fastDirs = 0
        if inputHDU is None:
            self.hdu = fits.open(path)
        else:
            self.hdu = inputHDU
        
        self.head = self.hdu[0].header
        self.data = self.hdu[0].data
        
        self.nint = self.head['NINT']
        self.ngroup = self.ngroup['NGROUP']
    
    def get_refpix_series(self,ampn=0):
        """ Get a reference pixel series for a given amplifier
        
        Parameters
        ------------
        ampn: int
            The amplifier (0-based counter) to get a time series for
        """
        if HDU[0].header['SUBARRAY'] == True:
            raise NotImplementedError
        else:
            ### length of amplifier in the fast direction
            self.ampFastSize = self.head['NAXIS1'] / self.nAmps
            self.ampStarts = np.arange(0,self.head['NAXIS1'],self.ampFastSize)
            self.ampEnds = self.ampStarts + self.ampFastSize - 1
            
            self.nRefPix = (self.nrefRows * self.head['NAXIS1'] * 2
                            + self.nrefCols * self.head['NAXSI1'] * 2
                            - self.nrefCols * self.nrefRows * 4)
            
            for oneInt in self.nint:
                thisIntStart = oneInt * self.ngroup
                oneInt = self.data[:,:,thisIntStart:thisIntStart]
                for oneGroup in self.ngroup:
                    oneImg = oneInt[:,:,oneGroup]
                    if ampDirs[ampn] == 1:
                        datBottom = oneImg[self.ampStarts[ampn]:self.ampEnds[ampn],0:self.nrefRows]
                        datTop = oneImg[self.ampStarts[ampn]:self.ampEnds[ampn],-self.nrefRows:-1]
                    else:
                        datBottom = oneImg[self.ampEnds[ampn]:self.ampStarts[ampn],0:self.nrefRows]
                        datTop = oneImg[self.ampEnds[ampn]:self.ampStarts[ampn],-self.nrefRows:-1]
                    if ampn == 0:
                        datSides = oneImg[0:self.nrefCols,self.nRefRows:-self.nRefRows]
                    elif ampn == 1:
                        datSides = oneImg[-self.nrefCols:-1,self.nRefRows:-self.nRefRows]
                    else:
                        datSides = None
                    pdb.set_trace()
        
        





