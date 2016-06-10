# Importing libraries
from utils import *
from scipy.interpolate import interp1d

def get_name_iso(base_name,age,metal):
	"""
	returns the name of an isochrone file for a given age and metallicity
	
	:param base_name: isochronal base name
	:type  base_name: string
	:param       age: isochronal age in yr
	:type        age: float
	:param     metal: isochronal metallicity
	:type      metal: float
	:return		: isochronal name
	:rtype		: string
	"""
	return(base_name+'_Z{}_t{:.1f}e9.dat'.format(metal, age/1e9))
		
def get_name_field(base_name,l,b):
	"""
	returns the name of a field file for given latitude and longitude
	
	:param base_name: field base name
	:type  base_name: string
	:param         l: galactic longitude in degrees
	:type	       l: float
	:param         b: galactic latitude in degrees
	:type	       b: float
	:return	        : field name
	:rtype          : string
	"""
	return(base_name+'_l{0}_b{1}.dat'.format(l, b))

##################
# Parameters  Iso#
############################################

isoc_ages = [1e9,3e9,5e9,7e9,9e9,11e9,13e9]#,9e9,10e9]
isoc_metal = [0.00015,0.0005,0.0015,0.005,0.015]#,0.02,0.03]
d_seq = [20000,25000,30000,40000,45000,60000,90000,100000,140000,150000,160000]
isoc_path = '/home/hperottoni/Documentos/Trabalhos_andamento/caca_fantasmas/MaGIK-edicao1/isocs/'
isoc_base_name = 'iso_magik'

#####################
# Parameters  fields#
############################################

gal_l = [90]
gal_b = [40]
field_path = '/home/hperottoni/Documentos/Trabalhos_andamento/caca_fantasmas/MaGIK-edicao1/canes/'
field_base_name = 'canes'

#####################
# Save to#
############################################
save_path = '/home/hperottoni/Documentos/Trabalhos_andamento/caca_fantasmas/MaGIK-edicao1/canes/'

def get_probs(isoc_path, isoc_base_name, isoc_ages, isoc_metal, d_seq, field_path, field_base_name, gal_l, gal_b, save_path):
    """
    For each star in each field, calculate the probability of the star belonging to each simple stellar population defined 
    by each age, metallicity and distance given. Saves one file for each ssp containning the obtained probabilities for all
    stars.
    
    :param       isoc_path: path to the isochrone files
    :param  isoc_base_name: isochronal base name
    :param       isoc_ages: list of all ages to be considered. Ages in yr.
    :param      isoc_metal: list of all metallicities (Z) to be considered.
    :param           d_seq: list of all distances to be considered. Distances in parsecs.
    :param      field_path: path to the field files
    :param field_base_name: field base name
    :param           gal_l: list of galactic longitudes to be considered.
    :param           gal_b: list of galactic latitudes to be considered.
    :param       save_path: path to where result files are to be saved
    
    :type        isoc_path: string
    :type   isoc_base_name: string
    :type        isoc_ages: list of floats
    :type       isoc_metal: list of floats
    :type            d_seq: list of floats
    :type       field_path: string
    :type  field_base_name: string
    :type            gal_l: list of floats
    :type            gal_b: list of floats
    :type        save_path: string
    
    obs: stars with high magnitude error's are removed.
    """
    
    for lon in gal_l: # for each given longitude
        for lat in gal_b: # for each given latitude
       
            # Load field data
            field_name = get_name_field(field_base_name, l = lon, b = lat)
            # TODO: dar um jeito de permitir que os argumentos de np.loadtxt sejam passados nos argumentos da funcao
            field_data = np.loadtxt(field_path+field_name, delimiter = ",", skiprows = 1)
             
            obj_lon = field_data[:,1]
            obj_lat = field_data[:,2]
            obj_mag_g = field_data[:,3]
            obj_mag_r = field_data[:,4]
            obj_mag_i = field_data[:,5]
            obj_mag_g_err = field_data[:,6]
            obj_mag_r_err = field_data[:,7]
            obj_mag_i_err = field_data[:,8]
             
            # Calculate colors from magnitudes
            obj_cor_gr = obj_mag_g - obj_mag_r
            obj_cor_gi = obj_mag_g - obj_mag_i
             
            # Adopting a color error
            obj_mag_r_err = obj_mag_r_err/4 # TODO: descobrir porque estamos usando isso
            obj_cor_gr_err = 0.1*obj_mag_r_err
             
            obj_mag_i_err = obj_mag_i_err/4 # TODO: descobrir porque estamos usando isso
            obj_cor_gi_err = 0.1*obj_mag_i_err
             
            # Removing stars with high magnitude error
            # TODO: fazer este passo diretamente em field_data para nao ter que escrever linha por linha
            obj_mag_filt_high_error = obj_mag_r_err < 0.2
            obj_lon = obj_lon[obj_mag_filt_high_error]
            obj_lat = obj_lat[obj_mag_filt_high_error]
            obj_mag_g = obj_mag_g[obj_mag_filt_high_error]
            obj_mag_r = obj_mag_r[obj_mag_filt_high_error]
            obj_mag_i = obj_mag_i[obj_mag_filt_high_error]
            obj_mag_g_err = obj_mag_g_err[obj_mag_filt_high_error]
            obj_mag_r_err = obj_mag_r_err[obj_mag_filt_high_error]
            obj_mag_i_err = obj_mag_i_err[obj_mag_filt_high_error]           
            obj_cor_gr = obj_cor_gr[obj_mag_filt_high_error]
            obj_cor_gi = obj_cor_gi[obj_mag_filt_high_error]
            obj_cor_gr_err = obj_cor_gr_err[obj_mag_filt_high_error]
            obj_cor_gi_err = obj_cor_gi_err[obj_mag_filt_high_error]
             
            for age in isoc_ages: # for each age
            	for metal in isoc_metal: # for each metallicity
                     	    
            	    # Load isochrone file
            	    isoc_name = get_name_iso(isoc_base_name,age,metal)
            	    isoc_data = np.loadtxt(isoc_path+isoc_name) # TODO: dar um jeito de poder passar argumentos para esta funcao
                     	    
            	    isoc_mass = isoc_data[:,2]
            	    isoc_mag_g = isoc_data[:,9]
            	    isoc_mag_r = isoc_data[:,10]
            	    isoc_mag_i = isoc_data[:,11]
                      	    
            	    # Interpolating isochrone
            	    isoc_mag_g_interp_fun = interp1d(isoc_mass, isoc_mag_g)
            	    isoc_mag_r_interp_fun = interp1d(isoc_mass, isoc_mag_r)
            	    isoc_mag_i_interp_fun = interp1d(isoc_mass, isoc_mag_i)
                      	    
            	    isoc_mass = np.linspace(isoc_mass[0], isoc_mass[-1], 10000)
            	    isoc_mag_g = isoc_mag_g_interp_fun(isoc_mass)
            	    isoc_mag_r = isoc_mag_r_interp_fun(isoc_mass)
            	    isoc_mag_i = isoc_mag_i_interp_fun(isoc_mass)
            	    
            	    # Calculating isochrone colors
            	    isoc_cor_gr = isoc_mag_g - isoc_mag_r
            	    isoc_cor_gi = isoc_mag_g - isoc_mag_i

            	    for dist in d_seq: # for each distance
            	        isoc_mag_r_dist = 5*(np.log10(dist)) - 5 + isoc_mag_r
            	        isoc_mag_i_dist = 5*(np.log10(dist)) - 5 + isoc_mag_i
            	        
            	        print(age, metal, dist)
            	        # Obtain the probabilities of each star belonging to the ssp defined by age, metal and dist using the
            	        # two different colors
            	        Probs_gr_r = P_multiple(cor_iso = isoc_cor_gr, 
            	                                cor_obs = obj_cor_gr, 
            	                                cor_err = obj_cor_gr_err, 
            	                                mag_iso = isoc_mag_r_dist,
            	                                mag_obs = obj_mag_r,
            	                                mag_err = obj_mag_r_err,
                                                m_iso = isoc_mass)
                         	                                
            	        Probs_gi_i = P_multiple(cor_iso = isoc_cor_gi, 
            	                                cor_obs = obj_cor_gi, 
            	                                cor_err = obj_cor_gi_err, 
            	                                mag_iso = isoc_mag_i_dist,
            	                                mag_obs = obj_mag_i,
            	                                mag_err = obj_mag_i_err,
                                                m_iso = isoc_mass)
                        
                        # This solves a problem with the plot
            	        Probs_gr_r_min = (Probs_gr_r[Probs_gr_r != 0]).min()
            	        for k in range(len(Probs_gr_r)):
            	            if Probs_gr_r[k] == 0: Probs_gr_r[k] = Probs_gr_r_min
            	            
            	        Probs_gi_i_min = (Probs_gi_i[Probs_gi_i != 0]).min()
            	        for k in range(len(Probs_gi_i)):
            	            if Probs_gi_i[k] == 0: Probs_gi_i[k] = Probs_gi_i_min
            	        
            	        # Normalizing the Probs
            	        Probs_gr_r = Probs_gr_r/Probs_gr_r.max()
            	        Probs_gi_i = Probs_gi_i/Probs_gi_i.max()
            	        
            	        # Preparing data to be saved
            	        save_data = np.zeros((len(Probs_gr_r), 10))
            	        save_data[:,0] = np.round(obj_lon,3)
            	        save_data[:,1] = np.round(obj_lat, 3)
            	        save_data[:,2] = np.round(obj_mag_r, 2)
            	        save_data[:,3] = np.round(obj_mag_r_err,2)
            	        save_data[:,4] = np.round(obj_mag_i,2)
            	        save_data[:,5] = np.round(obj_mag_i_err,2)
            	        save_data[:,6] = np.round(obj_cor_gr,2)
            	        save_data[:,7] = np.round(obj_cor_gi,2)
            	        save_data[:,8] = Probs_gr_r
            	        save_data[:,9] = Probs_gi_i
            	        
            	        # Saving data
            	        save_filename = save_path+'probs_l{0}_b{1}_age{2}e9_Z{3}_d{4}kpc.dat'.format(lon,lat,age/1e9,metal,dist/1000)
                        np.savetxt(save_filename, save_data, delimiter = ',')
						
                        # Plotting data
#                        plt.scatter(obj_cor_gr, obj_mag_r, c = np.log(Probs_gr_r), cmap = cm.jet, marker='o')
#                        plt.plot(isoc_cor_gr, isoc_mag_r_dist)
#                        plt.gca().invert_yaxis()
#                        plt.show()

get_probs(  isoc_path=isoc_path, 
            isoc_base_name=isoc_base_name, 
            isoc_ages=isoc_ages, 
            isoc_metal=isoc_metal, 
            d_seq=d_seq, 
            field_path=field_path, 
            field_base_name=field_base_name, 
            gal_l=gal_l, 
            gal_b=gal_b, 
            save_path=save_path)
