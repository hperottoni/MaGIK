import numpy as np
import random
from scipy.interpolate import interp1d
import csv
import fileinput

save_name_background = 'backgroung.txt'
save_name_spp = 'ssp.txt'
save_name_out = 'out.dat'
background_name = 'field1.dat'
iso_name = 'iso_12.00e9.dat'
N = 100 #numero de objetos da populacao

#Carrega dados do Trilegal como background
data = np.loadtxt(background_name, delimiter = ',')#,skiprows=1)

Longitude = data[:,21]
Latitude  = data[:,22]
print(len(Longitude))
Mag_g = data[:,39]
Mag_r = data[:,40]
Mag_i = data[:,41]

err_Mag_g = data[:,36]
err_Mag_r = data[:,37]
err_Mag_i = data[:,38]
Cor_obs = Mag_g-Mag_r

save_array = np.zeros((len(err_Mag_g), 9))
save_array[:,0] = Longitude
save_array[:,1] = Latitude
save_array[:,2] = Mag_g
save_array[:,3] = Mag_r
save_array[:,4] = Mag_i
save_array[:,5] = err_Mag_g
save_array[:,6] = err_Mag_r
save_array[:,7] = err_Mag_i
save_array[:,8] = Cor_obs
#Salva background do Trilegal
np.savetxt(save_name_background, save_array, delimiter = ',')


#Carrega dados da isocrona

data_iso = np.loadtxt(iso_name)

mag_iso_g = data_iso[:,9]
mag_iso_r = data_iso[:,10]
mag_iso_i = data_iso[:,11]
mag_err_g = 0.01*(mag_iso_g)
mag_err_r = 0.01*(mag_iso_r)
mag_err_i = 0.01*(mag_iso_i)
mag_iso_g2 = random.gauss(mag_iso_g,mag_err_g)
mag_iso_r2 = random.gauss(mag_iso_r,mag_err_r)
mag_iso_i2 = random.gauss(mag_iso_i,mag_err_i)
cor_iso_gr2 = mag_iso_g2-mag_iso_r2
cor_iso   = mag_iso_g2-mag_iso_r2
mag_iso_r2 = 5*(np.log10(20000)) - 5 + mag_iso_r


args = np.random.randint(0, len(cor_iso_gr2), N)
#distribui aleatoariamente objetos dessa isocrona em uma regiao (coords l,b)
longitude = [round(random.uniform(132.7,132.8),10) for _ in range(N)]
#latitude = -56.5+(np.random.random(N)/8)
latitude = [round(random.uniform(-56.7,-56.5),10) for _ in range(N)]

save_array = np.zeros(((N), 9))

save_array[:,0] = longitude
save_array[:,1] = latitude
save_array[:,2] = mag_iso_g2[args]
save_array[:,3] = mag_iso_r2[args]
save_array[:,4] = mag_iso_i2[args]
save_array[:,5] = mag_err_g[args]
save_array[:,6] = mag_err_r[args]
save_array[:,7] = mag_err_i[args]
save_array[:,8] = cor_iso[args]

#Salva a ssp gerada em um arquivo
np.savetxt(save_name_spp, save_array, delimiter = ',')

#Salva background e ssp num unico arquivo.
filenames = [save_name_spp, save_name_background]

fout2 = open(save_name_out, 'w')
csvfout2 = csv.writer(fout2, delimiter=',')
for line in fileinput.input(filenames): # list of file names  as strings
    fout2.write(line)
fout2.close()
