#!/Users/jeggen/miniconda2/envs/iraf27/bin/python

# Plot the optical lightcure for OE 110 for the following instruments
# - 72"
# - 42"
# - 31"
# - SMARTS

import os
import numpy
import matplotlib
import matplotlib.pyplot

chart_title = 'OE 110: 2008-2018 R-band Photometry'
file72 = '/Users/jeggen/Optical/oe110_recal/72_inch/combined/date_mag_err.dat'
file42 = '/Users/jeggen/Optical/oe110_recal/42_inch/combined/date_mag_err.dat'
file31 = '/Users/jeggen/Optical/oe110_recal/31_inch/combined/date_mag_err.dat'
fileSMARTS = '/Users/jeggen/Optical/oe110_recal/SMARTS/combined/date_mag_err.dat'
filePol = '/Users/jeggen/Optical/oe110_recal/Polarimetry/recal/date_mag_pol_evpa.dat'

toffset = 2400000
xmin = 54600 # in MJD (time = JD - 2400000)
xmax = 58200
ymin = 21.7
ymax = 16

# extract data from files
time72,mag72,magErr72 = numpy.loadtxt(file72, unpack=True)
time42,mag42,magErr42 = numpy.loadtxt(file42, unpack=True)
time31,mag31,magErr31 = numpy.loadtxt(file31, unpack=True)
SMtime,SMmag,SMmagErr = numpy.loadtxt(fileSMARTS, unpack=True)
Poltime,Polmag,PolmagErr,PolPer,PolPerErr,PolEVPA,PolEVPAErr = numpy.loadtxt(filePol, unpack=True)

# Plot it all up
fig, plot = matplotlib.pyplot.subplots(nrows=1, ncols=1, sharex=True, figsize=(7,2))
plot.set_title(chart_title)
plot.set_ylabel('R mag.')
plot.set_xlim(xmin,xmax)
plot.set_ylim(ymin,ymax)
plot.errorbar(SMtime - toffset,SMmag,yerr=SMmagErr,marker = '.',fmt='go',label="SMARTS")
plot.errorbar(Poltime - toffset,Polmag,yerr=PolmagErr,marker = '.',fmt='mx',label='72" (Pol.)')
plot.errorbar(time72 - toffset,mag72,yerr=magErr72,marker = '.',fmt='ko',label='72" (Phot.)')
plot.errorbar(time42 - toffset,mag42,yerr=magErr42,marker = '.',fmt='ro',label='42"')
plot.errorbar(time31 - toffset,mag31,yerr=magErr31,marker = '.',fmt='bo',label='31"')
#plot.locator_params(nbins=6)
plot.legend(loc='best')
plot.set_xlabel('Time (MJD)')

matplotlib.pyplot.show()

