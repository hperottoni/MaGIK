from utils import *
from scipy.interpolate import interp1d
import random

def get_name_iso(base_name,age,metal):
	return(base_name+'_Z{0}_t{1}e9.dat'.format(metal, age/1e9))
		
def get_name_field(base_name,l,b):
	return(base_name+'_l{0}_b{1}.dat'.format(l, b))


##################
# Parameters  Iso#
############################################


ages = [8e9]#,9e9,10e9]
metal =[0.01]#,0.02,0.03]
d_seq = [20000,60000]
path = '/home/hperottoni/Documentos/Trabalhos_andamento/caca_fantasmas/SAB/trilegal/versao_final_simulada/'
base_name = '\\teste5'


#####################
# Parameters  fields#
############################################

gal_l = [90]
gal_b = [40]
#gal_b = [40,41]
path_field = '/home/hperottoni/Documentos/Trabalhos_andamento/caca_fantasmas/SAB/trilegal/versao_final_simulada/'
base_name_field = 'field'

for l in gal_l:
	
	for b in gal_b:

		#field_name = get_name_field(path_field+base_name_field,l,b)
		#field_name = 'new_field1.txt'
		field_name = 'out.dat'

		data = np.loadtxt(field_name, delimiter = ',')#,skiprows=1)


		#objid,l,b,g,r,i,psfMagErr_g,psfMagErr_r,psfMagErr_i


		#Gc,logAge,MeH,m_ini,logL,logT,logg,mM0,Av,m2m1,Mbol,u,g,r,i,zM,u10,g10,r10,i10,z10,l,b,Fld,dist,X,Y,Z,RG,Uv,Vv,Wv,Vr,PmRA,PmDec,Mact, errog,error,erroi,g0,r0,i0
		Mag_r = data[:,3]#[:,40]
		mask = Mag_r < 22.5
		Mag_r = Mag_r[mask]				
		Longitude = data[:,0][mask]#[:,21]
		Latitude  = data[:,1][mask]#[:,22]
		print(len(Longitude))
		Mag_g = data[:,2][mask]#[:,39]

		Mag_i = data[:,4][mask]#[:,41]

		err_Mag_g = data[:,5][mask]#[:,36]
		err_Mag_r = data[:,6][mask]#[:,37]
		err_Mag_i = data[:,7][mask]#[:,38]
		Cor_obs = data[:,8][mask]


		for i in range(len(ages)):

			for j in range(len(metal)):

				ages_i  = ages[i]
				metal_j = metal[j]

				#iso_name = get_name_iso(path+base_name,ages_i,metal_j) 
				iso_name = 'iso_12.00e9.dat'
				data_iso = np.loadtxt(iso_name)
				
				m_iso     = data_iso[:,2]
				mag_iso_g = data_iso[:,9]
				mag_iso_r = data_iso[:,10]
				mag_err_g = 0.01*(mag_iso_g)
				mag_err_r = 0.01*(mag_iso_r)

				mag_iso_g2 = random.gauss(mag_iso_g,mag_err_g)
				mag_iso_r2 = random.gauss(mag_iso_r,mag_err_r)
				cor_iso_gr2 = mag_iso_g2-mag_iso_r2

				mag_iso_r2 = 5*(np.log10(20000)) - 5 + mag_iso_r
				interp_mag_g = interp1d(m_iso, mag_iso_g)
				interp_mag_r = interp1d(m_iso, mag_iso_r)

				m_iso = np.linspace(m_iso[0],m_iso[-1],10000)
				mag_iso_g = interp_mag_g(m_iso)
				mag_iso_r = interp_mag_r(m_iso)	
				cor_iso   = mag_iso_g-mag_iso_r


				for d in d_seq:
					
					mag_g = Mag_g
					mag_r = Mag_r
					#mag_i = Mag_i
					latitude  = Latitude
					longitude = Longitude					

					err_mag_g = err_Mag_g
					err_mag_r = err_Mag_r
					#err_mag_i = err_Mag_i
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
					print Probs.max()
					print('aqui')
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
					#plt.plot(cor_iso,Mag_iso_r)
					plt.plot(cor_iso_gr2,mag_iso_r2, 'go')
					plt.gca().invert_yaxis()
					plt.show()
