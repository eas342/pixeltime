from astropy.io import fits
import numpy as np
import pdb
import matplotlib.pyplot as plt

#defaultPath = ('/data/External/ISIMCV3_unzipped/' +
#               'NRCNRCA4-DARK-6020143529_1_1_40269_JW1_JLAB88_20160120T143723.248_20160120T145653.554/' +
#               'NRCNRCA4-DARK-60201435291_1_484_SE_2016-01-20T15h12m28.fits')


defaultPath = ('/Users/everettschlawin/Documents/jwst/test_data/cv3/full_frame_wlp8_05A4/' +
               'NRCN821WLP8FULL5-6012184421_1_481_SE_2016-01-12T20h26m57_I025.fits')

nDrop = {'RAPID':0,'BRIGHT1':1,'BRIGHT2':0,'SHALLOW2':3,'SHALLOW4':1,
         'MEDIUM2':8,'MEDIUM8':2,'DEEP2':18,'DEEP8':12}

pixelRate = 1e-5 ## seconds per pixel

class exposure():
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
        self.origdata = self.hdu[0].data
        self.data = self.avgSub(self.origdata)
        self.nint = self.head['NINT']
        
        self.ngroup = self.head['NGROUP']
        self.nframe = self.head['NFRAME']
        self.ndrop = nDrop[self.head['READOUT']]
        
    
    def avgSub(self,origData):
        """ Subtract the average frame from all frames 
        This will get rid of the average frame
        permitting more careful study of the reference pixel time series
        """
        self.avgImg = np.nanmean(origData,axis=0)
        return origData - self.avgImg
        
    
    def get_refpix_series(self,ampn=0):
        """ Get a reference pixel series for a given amplifier
        
        Parameters
        ------------
        ampn: int
            The amplifier (0-based counter) to get a time series for
        """
        if self.head['SUBARRAY'] == True:
            raise NotImplementedError
        else:
            ### length of amplifier in the fast direction
            self.tFrame = 10.73677
            nReset = 1 ## number of reset frames
            self.ampFastSize = self.head['NAXIS1'] / self.nAmps
            self.ampStarts = np.arange(0,self.head['NAXIS1'],self.ampFastSize)
            self.ampEnds = self.ampStarts + self.ampFastSize - 1
            
            self.nRefPix = (self.nrefRows * self.head['NAXIS1'] * 2
                            + self.nrefCols * self.head['NAXIS1'] * 2
                            - self.nrefCols * self.nrefRows * 4)
            
            groupSeries, groupCount, intCount = [], [], []
            for oneInt in np.arange(self.nint):
                thisIntStart = oneInt * self.ngroup
                thisIntEnd = (oneInt + 1) * self.ngroup
                oneInt = self.data[thisIntStart:thisIntEnd,:,:]
                framesBefore = (self.ngroup * self.nframe + (self.ngroup - 1) * self.ndrop + nReset) * oneInt
                
                for oneGroup in np.arange(self.ngroup):
                    oneImg = oneInt[oneGroup,:,:]
                                ## frames coadded             ## frames dropped
                    frameStart = ((self.nframe * oneGroup) + oneGroup * self.ndrop + framesBefore) * self.tFrame
                    if self.ampDirs[ampn] == 1:
                        datBottom = oneImg[0:self.nrefRows,self.ampStarts[ampn]:self.ampEnds[ampn]]
                        datTop = oneImg[-self.nrefRows:-1,self.ampStarts[ampn]:self.ampEnds[ampn]]
                        timeBottom = np.arange(0,self.ampFastSize)
                    else:
                        datBottom = oneImg[0:self.nrefRows,self.ampEnds[ampn]:self.ampStarts[ampn]]
                        datTop = oneImg[-self.nrefRows:-1,self.ampEnds[ampn]:self.ampStarts[ampn]]
                    if ampn == 0:
                        datSides = oneImg[self.nrefRows:-self.nrefRows,0:self.nrefCols]
                    elif ampn == 1:
                        datSides = oneImg[self.nrefRows:-self.nrefRows,-self.nrefCols:-1]
                    else:
                        datSides = None
                    combSeries = [np.ravel(datBottom),np.ravel(datSides),np.ravel(datTop)]
                    allRefPix = np.hstack(combSeries)
                    groupSeries.append(allRefPix)
                    groupCount.append(oneGroup)
                intCount.append(oneInt)
                
        return intCount, groupCount, groupSeries
        
    def plot_refpix(self):
        """ Plot the reference pixel time series """
        intC, groupC, groupS = self.get_refpix_series()
        for oneGroup in groupC:
            plt.plot(groupS[oneGroup],rasterized=True,label='Grp {}'.format(oneGroup))
        plt.show()
        
def make_timeImage():
    """ 
    Makes a time coordinate image for full frame images
    """
    tD = np.zeros([2048,2048])
    ampSize = 512
    lineOverheads = 12
    directions = [1,-1,1,-1]
    for oneAmp in np.arange(4):
        
        startX, endX = ampSize * oneAmp, ampSize * (1 + oneAmp)
        
        for oneRow in np.arange(2048):
            oneAmpRow = oneRow * (lineOverheads + ampSize) + np.arange(ampSize)
            if directions[oneAmp] == 1:
                tD[oneRow,startX:endX] = oneAmpRow
            else:
                tD[oneRow,startX:endX] = oneAmpRow[::-1]
    
    HDU = fits.PrimaryHDU(tD)
    HDUList = fits.HDUList(HDU)
    HDUList.writeto('pixeltime/data/full_array_timing.fits',overwrite=True)