import scipy.ndimage as ndimage
import numpy as np
import matplotlib.cm as cm
from matplotlib import pyplot as plt
#from scipy.stats import gaussian_kde
from utils import *

star_filter  = np.loadtxt('probs_l90_b40_age8.0e9_Z0.01_d20kpc.dat',delimiter = ',', usecols=(0,1,2,3,4)) 

min_probs = star_filter[:,4].min()

star_filter_min = star_filter[star_filter[:,4] == min_probs]
star_filter     = star_filter[star_filter[:,4] != min_probs]

percentis = np.percentile(star_filter[:,4],[20,40,60,80,90])

kernels = np.zeros(len(star_filter[:,4]))
kernels_values = [25,20,15,7,6,5]


for j in range(len(star_filter[:,4])):
	if (star_filter[j,4] <= percentis[0]):
		kernels[j] = kernels_values[0]


for i in range(1,len(percentis)):
	for j in range(len(star_filter[:,4])):
		if (star_filter[j,4] > percentis[i-1]) and (star_filter[j,4] <= percentis[i]):
			kernels[j] = kernels_values[i]

for j in range(len(star_filter[:,4])):
	if (star_filter[j,4] > percentis[-1]):
		kernels[j] = kernels_values[-1]
		
h, x, y, p = plt.hist2d(star_filter_min[:,1],star_filter_min[:,0], bins = 130)

gaussian  = ndimage.gaussian_filter(h, 40)


for i in range(len(kernels_values)):
	
	star_filteri = star_filter[kernels == kernels_values[i], :]

	h, x, y, p = plt.hist2d(star_filteri[:,1],star_filteri[:,0], bins = 130)
	gaussian_i  = ndimage.gaussian_filter(h, kernels_values[i])

	gaussian = gaussian+gaussian_i


gaussian = gaussian / (130*130*gaussian.sum())
print(gaussian.max())

#################################
hbg, x, y, p = plt.hist2d(star_filter_min[:,1],star_filter_min[:,0], bins = 130)

kernels_values_bg = [15,15,15,15,15,15]
gaussianbg  = ndimage.gaussian_filter(hbg, 15)


for i in range(len(kernels_values)):
	
	star_filteri = star_filter[kernels == kernels_values[i], :]

	hbg, x, y, p = plt.hist2d(star_filteri[:,1],star_filteri[:,0], bins = 130)
	gaussianbg_i  = ndimage.gaussian_filter(hbg, kernels_values_bg[i])

	gaussianbg = gaussianbg+gaussianbg_i
####################################


A_mean, A_sdev = running_mean_sdev(gaussian, nline = 40)

Mean   = gaussian-A_mean
desvio = (gaussian-A_mean)/A_sdev

print(Mean[15:115, 15:115].max())
print(Mean[15:115, 15:115].min())
print(desvio[15:115, 15:115].max())
print(desvio[15:115, 15:115].min())
view_xmin = min(star_filter[:,0])
view_xmax = max(star_filter[:,0])
view_ymin = min(star_filter[:,1])
view_ymax = max(star_filter[:,1])




#fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(100,100))

#axes[0, 0].imshow(gaussian, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#axes[0, 0].set(title='Weight > 0.7')
#plt.minorticks_on()
#plt.grid(True)

#axes[0, 1].imshow(Mean, interpolation='nearest', cmap = plt.cm.cubehelix, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#axes[0, 1].set(title='Weight < 0.2')
#plt.minorticks_on()
#plt.grid(True)

#axes[0, 2].imshow(desvio, interpolation='nearest', cmap = plt.cm.cubehelix, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#axes[0, 2].set(title='Residual')
#plt.minorticks_on()
#plt.grid(True)

#axes[1,0].imshow(gaussian, interpolation='nearest', extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#axes[1,0].set(title='Weight > 0.7')
#plt.minorticks_on()
#plt.grid(True)

#axes[1,1].imshow(gaussian, cmap = plt.cm.cubehelix, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#axes[1,1].set(title='Weight < 0.2')
#plt.minorticks_on()
#plt.grid(True)

#axes[1,2].imshow(gaussian, interpolation='nearest', cmap = plt.cm.cubehelix, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#axes[1,2].set(title='Residual')
#plt.minorticks_on()
#plt.grid(True)
##plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.05, hspace=.1)

##fig.tight_layout()
#plt.show()




plt.subplots(nrows=3, ncols=2,figsize=(18,14))
plt.subplots_adjust(left  = 0.1,right = 0.94,bottom = 0.1,top = 0.9,wspace = 0.15,hspace = 0.2)

plt.subplot2grid((2, 3), (0,0))
plt.imshow(gaussian, cmap = plt.cm.cubehelix, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
plt.yscale('linear')
plt.ylabel("latitude", fontsize=12)
plt.xlabel("longitude", fontsize=12)
plt.title('Filtro')
plt.grid(True)
plt.minorticks_on()

plt.subplot2grid((2, 3), (0, 1))
plt.imshow(Mean, cmap = plt.cm.cubehelix, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#plt.yscale('Filtro')
plt.ylabel("latitude", fontsize=12)
plt.xlabel("longitude", fontsize=12)
plt.title('Mean')
plt.grid(True)
plt.minorticks_on()

plt.subplot2grid((2, 3), (0,2))
plt.imshow(desvio, cmap = plt.cm.cubehelix, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#plt.yscale('linear')
plt.ylabel("latitude", fontsize=12)
plt.xlabel("longitude", fontsize=12)
plt.title('desvio')
plt.grid(True)
plt.minorticks_on()


plt.subplot2grid((2, 3), (1, 0))
plt.imshow(gaussian, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#plt.yscale('Mean')
plt.ylabel("latitude", fontsize=12)
plt.xlabel("longitude", fontsize=12)
plt.title('Filtro')
plt.grid(True)
plt.minorticks_on()

plt.subplot2grid((2, 3), (1, 1))
plt.imshow(Mean, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
#plt.yscale('desvio')
plt.ylabel("latitude", fontsize=12)
plt.xlabel("longitude", fontsize=12)
plt.title('Mean')
plt.grid(True)
plt.minorticks_on()




plt.subplot2grid((2, 3), (1, 2))
plt.imshow(desvio, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
plt.yscale('linear')
plt.ylabel("latitude", fontsize=12)
plt.xlabel("longitude", fontsize=12)
plt.title('Desvio')
plt.grid(True)
plt.minorticks_on()



plt.savefig('figureb.png', format='png')
plt.cla()
plt.clf()
plt.close()
