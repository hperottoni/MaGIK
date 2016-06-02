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
	return(base_name+'_Z{0}_t{1}e9.dat'.format(metal, age/1e9))
		
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

ages = [8e9]#,9e9,10e9]
metal = [0.01]#,0.02,0.03]
d_seq = [20000,60000]
path = '/'
base_name = 'isoc'

#####################
# Parameters  fields#
############################################

gal_l = [90]
gal_b = [40]
path_field = '/'
base_name_field = 'field'

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
    """
    
    for lon in gal_l: # for each given longitude
        for lat in gal_b: # for each given latitude
            
            # Load field data
            field_name = get_name_field(field_base_name, l = lon, b = lat)
            # TODO: dar um jeito de permitir que os argumentos de np.loadtxt sejam passados nos argumentos da funcao
            field_data = np.loadtxt(field_name, delimiter = ",", skiprows = 1)
            
            obj_lon = field_data[:,1]
            obj_lat = field_data[:,2]
            obj_mag_g = field_data[:,3]
            obj_mag_r = field_data[:,4]
            obj_mag_i = field_data[:,5]
            obj_mag_g_err = field_data[:,6]
            obj_mag_r_err = field_data[:,7]
            obj_mag_i_err = field_data[:,8]
            
for l in gal_l:

    for b in gal_b:

		field_name = 'field0.dat' 
		data = np.loadtxt(field_name, delimiter = ',',skiprows=1)


		#objid,l,b,g,r,i,psfMagErr_g,psfMagErr_r,psfMagErr_i
		
		Longitude = data[:,1]
		Latitude  = data[:,2]
		Mag_g = data[:,3]
		Mag_r = data[:,4]
		Mag_i = data[:,5]

		err_Mag_g = data[:,6]
		err_Mag_r = data[:,7]
		err_Mag_i = data[:,8]
		Cor_obs = Mag_g-Mag_r

		for i in range(len(ages)):

			for j in range(len(metal)):

				ages_i  = ages[i]
				metal_j = metal[j]

				iso_name = 'iso_12.00e9.dat'
		
				data_iso = np.loadtxt(iso_name)
				
				m_iso     = data_iso[:,2]
				mag_iso_g = data_iso[:,9]
				mag_iso_r = data_iso[:,10]
				#cor_iso   = mag_iso_g-mag_iso_r
				
				interp_mag_g = interp1d(m_iso, mag_iso_g)
				interp_mag_r = interp1d(m_iso, mag_iso_r)

				m_iso = np.linspace(m_iso[0],m_iso[-1],10000)
				mag_iso_g = interp_mag_g(m_iso)
				mag_iso_r = interp_mag_r(m_iso)	
				cor_iso   = mag_iso_g-mag_iso_r
				

				for d in d_seq:
					
					mag_g = Mag_g
					mag_r = Mag_r
					mag_i = Mag_i
					latitude  = Latitude
					longitude = Longitude					

					err_mag_g = err_Mag_g
					err_mag_r = err_Mag_r
					err_mag_i = err_Mag_i
					cor_obs   = Cor_obs


					Mag_iso_r = 5*(np.log10(d)) - 5 + mag_iso_r
					mag_err = err_mag_r
					mag_filt = mag_err < 0.2

					longitude = longitude[mag_filt]
					latitude  = latitude[mag_filt]
					mag_err = mag_err[mag_filt]
					mag_err = mag_err/2.
					mag_r = mag_r[mag_filt]
					cor_obs = cor_obs[mag_filt]
					cor_err = 0.15*mag_err


					print(ages_i,metal_j,d)

					Probs = P_multiple(cor_iso, cor_obs, cor_err, Mag_iso_r, mag_r, mag_err, m_iso, a = 2.7)

					P_min = (Probs[Probs != 0]).min()
					for k in range(len(Probs)):
					    if Probs[k] == 0: Probs[k] = P_min
					Probs = Probs/Probs.max()


					save_array = np.zeros((len(Probs), 5))
					save_array[:,0] = longitude
					save_array[:,1] = latitude			
					save_array[:,2] = mag_r
					save_array[:,3] = cor_obs
					save_array[:,4] = Probs

					save_filename = 'probs_l{0}_b{1}_age{2}e9_Z{3}_d{4}kpc.dat'.format(l,b,ages_i/1e9,metal_j,d/1000)

					np.savetxt(save_filename, save_array, delimiter = ',')


					plt.scatter(cor_obs, mag_r, c = np.log(Probs), cmap = cm.jet, marker='o')
					plt.plot(cor_iso,Mag_iso_r)
					plt.gca().invert_yaxis()
					plt.show()
