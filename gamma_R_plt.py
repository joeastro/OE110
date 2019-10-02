#!/Users/jeggen/miniconda2/envs/iraf27/bin/python

# Program to make a plot of Gamma-Ray and R-band data for one target.  Both 
# Will be converted to the same units and plotted on the same plot window.  

# Author: Joseph R. Eggen, Graduate Student at Georgia State U.

# Date: Oct 2012
###############################################################################

import os
import numpy
import matplotlib
import pylab

chart_title = r'OE 110: 2008-2018 R-band & $\gamma$-ray Data'
r_band_data_file = 'all_R_date_mag_err.dat'
#r_band_data_file = 'MASTER_sort_date_mag_err.dat'
#gamma_data_file = 'light_curve_Eb_composit.txt'
gamma_data_file = 'light_curve_monthly_logpar_1_115b.txt'
#gamma_data_file = 'flux_lc_fermi_no_sun_daily.dat'
#gamma_data_file = 'light_curve_f.txt'
#gamma_bin_size = 2551443.0/4
gamma_bin_size = 28*24*60*60
norm_scale = 1e-10
opt_flux_scale = 1e-5
gamma_flux_scale = 1e-8
energy_start = 100
energy_stop = 300000
xmin = 55050
xmax = 58200
xerr_bar = gamma_bin_size/(2.*60*60*24)
UL_Threshold = 9.0

# set limits for 3 different panel sets
tlength = xmax - xmin
xmin1 = xmin
#xmax1 = xmin + tlength/2
#xmin2 = xmax1
#xmax2 = xmax1 + tlength/2
#xmin3 = xmax2
#xmax3 = xmax
xmax1 = tlength/2 + xmin
xmin2 = xmax1
xmax2 = xmax

# extract necessary values from data files
g_data = numpy.loadtxt(gamma_data_file)
time_g,TS,Flux,FluxErr,Norm,NormErr,Alpha,AlphaErr,Beta,BetaErr,Eb,EbErr,UpperLim,RetCode,BinWidth = numpy.loadtxt(gamma_data_file, unpack=True)

time_o,Mag,MagErr = numpy.loadtxt(r_band_data_file, unpack=True)

# determine how many gamma-ray data points we have to work with
gamma_range=range(0,len(time_g))
x = pylab.arange(energy_start,energy_stop,((energy_stop-energy_start)/(len(gamma_range))))
# same for the optical data
optical_range=range(0,len(time_o))


# format optical data
r_flux_mJy = []
r_fluxErr_mJy = []
r_flux_MeV = []
r_fluxErr_MeV = []
time_o_mjd = []
for i in optical_range:
	r_flux_mJy.append(2941*10**(-0.4*Mag[i]))
#	r_fluxErr_mJy.append((2941**2)*(10**(-0.4*Mag[i]))*(0.4)*pylab.log(10)*MagErr[i])
	r_fluxErr_mJy.append(r_flux_mJy[i]*MagErr[i])
	r_flux_MeV.append(r_flux_mJy[i]*6.2415e-21)
	r_fluxErr_MeV.append(r_fluxErr_mJy[i]*6.2415e-21)
	time_o_mjd.append(time_o[i]-2.4e6)
	r_flux_MeV[i] = r_flux_MeV[i]*1/opt_flux_scale
	r_fluxErr_MeV[i] = r_fluxErr_MeV[i]*1/opt_flux_scale
	r_flux_mJy[i] = r_flux_mJy[i]*1/opt_flux_scale
	r_fluxErr_mJy[i] = r_fluxErr_mJy[i]*1/opt_flux_scale

# format gamma-ray data
lg_ratio = []
phot_flux = []
energy_flux = []
energy_fluxErr = []
time_g_mjd = []
UpperLimArray = []
UpperLimTime = []
FluxOnly = []
FluxOnlyErr = []
FluxOnlyTime = []
for bn in gamma_range:
	Norm[bn] = Norm[bn]*norm_scale
	NormErr[bn] = NormErr[bn]*norm_scale
	lg_ratio.append(pylab.log(x[bn]/Eb[bn]))
	phot_flux.append(Norm[bn]*(Eb[bn]/Eb[bn])**(-(Alpha[bn]+Beta[bn]*lg_ratio[bn])))
	energy_flux.append(Norm[bn]*Eb[bn])
	energy_fluxErr.append(NormErr[bn]*Eb[bn])
	# convert FERMI timests from sec since mission start (jan01,2001) to MJD
	current_time = float(((time_g[bn]+(BinWidth[bn]/2))/(60*60*24)+51910.5))
	time_g_mjd.append(current_time)
	energy_flux[bn] = energy_flux[bn]*1/gamma_flux_scale
	energy_fluxErr[bn] = energy_fluxErr[bn]*1/gamma_flux_scale
	CurrentFlux = Flux[bn]*1/gamma_flux_scale
	Flux[bn] = CurrentFlux
	CurrentFluxErr = FluxErr[bn]*1/gamma_flux_scale
	FluxErr[bn] = CurrentFluxErr
	CurrentUpperLim = UpperLim[bn]*1/gamma_flux_scale
	UpperLim[bn] = CurrentUpperLim
	if TS[bn] < UL_Threshold:
		UpperLimArray.append(CurrentUpperLim)
		UpperLimTime.append(current_time)
	else:
		FluxOnly.append(CurrentFlux)
		FluxOnlyErr.append(CurrentFluxErr)
		FluxOnlyTime.append(current_time)
UpperLimRange = range(0,len(UpperLimArray))

# ouput 3 data files containing the MJD & Flux formatted data
os.system('rm -f gamma_mjd_fluxPh_err.dat')
g_output = open('gamma_mjd_fluxPh_err.dat', 'a')
for i in gamma_range:
	gamma_data_line = time_g_mjd[i],Flux[i],FluxErr[i],0,0
	g_output.write(' '.join(map(str, gamma_data_line)) + '\n')
g_output.close

os.system('rm -f optical_mjd_fluxmJy_err.dat')
o_output = open('optical_mjd_fluxmJy_err.dat', 'a')
for i in optical_range:
	optical_data_line = time_o_mjd[i],r_flux_mJy[i],r_fluxErr_mJy[i],0,0
	o_output.write(' '.join(map(str, optical_data_line)) + '\n')
o_output.close

# plot the gamma-ray data
fig, plot = matplotlib.pyplot.subplots(nrows=4, ncols=1, sharex=False, figsize=(7,7))
plot1 = plot[0]
plot1.set_ylabel(r'$\gamma$', color='k')
plot1.set_title(chart_title)
plot1.locator_params(nbins=8)
plot1.set_ylim(-0.01,45.5)
plot1.set_xlim(xmin1,xmax1)
plot1.set_xticklabels([])

#plot1.errorbar(time_g_mjd,Flux,yerr=FluxErr,marker = '.', fmt='ko')
plot1.errorbar(FluxOnlyTime,FluxOnly,yerr=FluxOnlyErr,marker = '.', fmt='ko')
for bn in UpperLimRange:
	plot1.arrow(UpperLimTime[bn],UpperLimArray[bn],0,-0.8,width=.01,fc="r",ec="r",head_width=6.0,head_length=0.15)
	plot1.errorbar(UpperLimTime[bn],UpperLimArray[bn],marker='',fmt='ro')

# plot the optical data on the same plot, same x-axis but diff. y-axis
plot2 = plot[1]
plot2.yaxis.set_label_position('left')
plot2.errorbar(time_o_mjd,r_flux_mJy,yerr=r_fluxErr_mJy,marker = 'x', fmt='bo')
plot2.set_ylabel('R', color='k')
#plot2.set_ylim(.401,70.01)
plot2.set_ylim(-.01,40.01)
plot2.set_xlim(xmin1,xmax1)
plot2.locator_params(nbins=6)
fig.tight_layout()
#plot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

plot3 = plot[2]
plot3.set_ylabel(r'$\gamma$', color='k')
plot3.locator_params(nbins=8)
plot3.set_ylim(-0.01,45.01)
plot3.set_xlim(xmin2,xmax2)
plot3.set_xticklabels([])

plot3.errorbar(FluxOnlyTime,FluxOnly,yerr=FluxOnlyErr,marker = '.', fmt='ko')
for bn in UpperLimRange:
	plot3.arrow(UpperLimTime[bn],UpperLimArray[bn],0,-0.8,width=.01,fc="r",ec="r",head_width=6.0,head_length=0.15)
	plot3.errorbar(UpperLimTime[bn],UpperLimArray[bn],marker='',fmt='ro')

# plot the optical data on the same plot, same x-axis but diff. y-axis
plot4 = plot[3]
plot4.yaxis.set_label_position('left')
plot4.errorbar(time_o_mjd,r_flux_mJy,yerr=r_fluxErr_mJy,marker = 'x', fmt='bo')
plot4.set_ylabel('R', color='k')
#plot4.set_ylim(.401,70.01)
plot4.set_ylim(-.01,40.01)
plot4.set_xlim(xmin2,xmax2)
plot4.locator_params(nbins=6)
plot4.set_xlabel('Time (MJD)')

matplotlib.pyplot.show()
