import scipy.ndimage as ndimage
import numpy as np
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from utils import *

def get_name_probs(base_name, gal_l, gal_b, age, metal, dist):
    """
    Returns th name of the probs file for coordinates gal_l, gal_b, and ssp defined by age, metal and dist
    :param base_name: probs base name file
    :param     gal_l: galactic longitude of the field in degrees
    :param     gal_b: galactic latitude of the field in degrees
    :param       age: age of the ssp in years
    :param     metal: metallicity (Z) of the ssp
    :param      dist: distance of the ssp in parsecs
    :return         : probs filename

    :type  base_name: string
    :type      gal_l: float
    :type      gal_b: float
    :type        age: float
    :type      metal: float
    :type       dist: float
    :rtype          : string
    """
    return (base_name+'_l{}_b{}_age{:.1f}e9_Z{}_d{:.0f}kpc.dat'.format(gal_l, gal_b, age/1e9, metal, dist/1e3))

def field_kernel(probs_base_name, 
                 probs_path,
                 gal_l, 
                 gal_b, 
                 age, 
                 metal, 
                 dist,
                 color = 'gr',
                 percentiles = [20,40,60,80,90],
                 kernels =  [25,20,15,7,6,5],
                 kernel_bg = 40,
                 bins = 130):
    
    # load file containing isochronal likelihoods
    filename = probs_path + get_name_probs(probs_base_name, gal_l, gal_b, age, metal, dist)
    probs_data = np.loadtxt(filename, delimiter = ',')

    star_gal_l = probs_data[:,0]
    star_gal_b = probs_data[:,1]
    star_mag_r = probs_data[:,2]
    star_mag_r_err = probs_data[:,3]
    star_mag_i = probs_data[:,4]
    star_mag_i_err = probs_data[:,5]
    star_cor_gr = probs_data[:,6]
    star_cor_gi = probs_data[:,7]

    if color == 'gr':
        star_probs = probs_data[:,8]
    elif color == 'gi':
        star_probs = probs_data[:,9]

    # Create two filters to separate values for which star_probs is minimmum from other values
    # filter_min is used to represent the background    
    filter_min = star_probs == star_probs.min()
    filter_else = star_probs != star_probs.min()
    

    # obtain percentils
    percentils = np.percentile(star_probs[filter_else], [20,40,60,80,90])

    # set kernels
    star_kernel = np.zeros(len(star_probs[filter_else]))
    for i in range(len(star_probs[filter_else])):
        if (star_probs[filter_else][i] <= percentils[0]):
            star_kernel[i] = kernels[0]

        elif (star_probs[filter_else][i] > percentils[-1]):
            star_kernel[i] = kernels[-1]

        else:
            for j in range(1, len(percentils)):
                if (star_probs[filter_else][i] > percentils[j-1]) and (star_probs[filter_else][i] <= percentils[j]):
                    star_kernel[i] = kernels[j]

    #print star_kernel
    # Background
    hbg, xbg, ybg, pbg = plt.hist2d(star_gal_l[filter_min], star_gal_b[filter_min], bins = bins)
    density_field_bg = ndimage.gaussian_filter(hbg, kernel_bg)
    
    # Density field
    density_field = density_field_bg
    
    for i in range(len(kernels)):
        kernel_filter = star_kernel == kernels[i]

        h, x, y, p = plt.hist2d(star_gal_l[filter_else][kernel_filter], star_gal_b[filter_else][kernel_filter], bins = bins)
        density_field_i = ndimage.gaussian_filter(h, kernels[i])
        density_field = density_field + density_field_i

    # Normalizar density_field
    density_field = density_field / (bins*bins*density_field.sum())

    return density_field, star_gal_b,star_gal_l

gal_l = 90
gal_b = 40
isoc_ages = [1e9,3e9,5e9,7e9,9e9,11e9,13e9]#,9e9,10e9]
isoc_metal = [0.00015,0.0005,0.0015,0.005,0.015]#,0.02,0.03]
d_seq = [20000,25000,30000,40000,45000,60000,90000,100000,140000,150000,160000]
save_path = '/home/hperottoni/Documentos/Trabalhos_andamento/caca_fantasmas/MaGIK-edicao1/canes/'

Nlines = len(isoc_ages)*len(isoc_metal)*len(d_seq)
results = np.zeros((Nlines, 11))

line = 0

for age in isoc_ages: # for each age
    for metal in isoc_metal: # for each metallicity  
        for dist in d_seq:

            results[line, 0] = gal_l
            results[line, 1] = gal_b
            results[line, 2] = age
            results[line, 3] = metal
            results[line, 4] = dist

            density_field,star_gal_b,star_gal_l = field_kernel(probs_base_name = 'probs',
                                         probs_path = '/home/hperottoni/Documentos/Trabalhos_andamento/caca_fantasmas/MaGIK-edicao1/canes/',
                                         gal_l = gal_l,
                                         gal_b = gal_b,
                                         age = age,
                                         metal = metal,
                                         dist = dist)

            A_mean, A_sdev = running_mean_sdev(density_field, nline = 40)

            Mean   = density_field-A_mean
            Desvio = (density_field-A_mean)/A_sdev

            view_xmin = min(star_gal_b)
            view_xmax = max(star_gal_b)
            view_ymin = min(star_gal_l)
            view_ymax = max(star_gal_l)

            plt.subplots(nrows=3, ncols=1,figsize=(18,14))
            plt.subplots_adjust(left  = 0.1,right = 0.94,bottom = 0.1,top = 0.9,wspace = 0.15,hspace = 0.2)

            plt.subplot2grid((1, 3), (0,0))
            plt.imshow(density_field, cmap = plt.cm.cubehelix, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
            plt.yscale('linear')
            plt.ylabel("latitude", fontsize=12)
            plt.xlabel("longitude", fontsize=12)
            plt.title('Filtro')
            plt.grid(True)
            plt.minorticks_on()

            plt.subplot2grid((1, 3), (0, 1))
            plt.imshow(Mean, cmap = plt.cm.cubehelix, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
            #plt.yscale('Filtro')
            plt.ylabel("latitude", fontsize=12)
            plt.xlabel("longitude", fontsize=12)
            plt.title('Mean')
            plt.grid(True)
            plt.minorticks_on()

            plt.subplot2grid((1, 3), (0,2))
            plt.imshow(Desvio, cmap = plt.cm.cubehelix, extent=[view_xmin, view_xmax, view_ymin, view_ymax])
            #plt.yscale('linear')
            plt.ylabel("latitude", fontsize=12)
            plt.xlabel("longitude", fontsize=12)
            plt.title('desvio')
            plt.grid(True)
            plt.minorticks_on()

            #plt.show()

            save_filename = save_path+'density_l{0}_b{1}_age{2}e9_Z{3}_d{4}kpc.png'.format(gal_l,gal_b,age/1e9,metal,dist/1000)
            plt.savefig(save_filename, format='png')
            plt.cla()
            plt.clf()
            plt.close()

            line = line + 1

            results[line, 5]= density_field.max()
            results[line, 6]= density_field.min()
            results[line, 7]= Mean.max()
            results[line, 8]= Mean.min()
            results[line, 9]= Desvio.max()
            results[line, 10]= Desvio.min()

np.savetxt(save_path+'results.dat', results, delimiter = ',')
