### Pixel Time Series

This set of scripts is designed to create a time series of the pixels within an image

## Installation

```bash
git clone https://github.com/eas342/pixeltime
cd pixeltime
python setup.py install
```

## Basic Usage
Here's how to get a table of the pixel

```python
import pixeltime
import numpy as np

## here we make up fake data where the flux is the Y or X coordinate
## in reality, you could put measured data here
y, x = np.mgrid[0:2048,0:2048]
outT = pixeltime.main.all_pixels_tser(y)

```
`outT` is an astropy table with a time column and the flattend data for each output amplifier

You can plot the time series like so:

```python
import matplotlib.pyplot as plt
plt.plot(outT['time'],outT['Amp 0 Val'])
plt.show()
```
