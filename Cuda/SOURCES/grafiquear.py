import matplotlib.pyplot as plt
import numpy as np
from glob import glob

directorio = glob("TIMES/*")
totalcelulas = 0
totaltiempos = 0
meancelulas = []
meantiempos = []
for fichero in directorio:
	celulas = [0]
	tiempos = [0]
	dias = [1]
	with open(fichero, "r+") as file1:
		for line in file1:
			celulas = np.append(celulas, int(line.split()[0]) )
			tiempos = np.append(tiempos, int(line.split()[1]) )
			dias = np.append(dias, dias[-1]+1)
		
		plt.plot(dias, celulas, 'ko')
		plt.title('Tiempo de ejecución según el número de células cancerígenas')
		plt.xlabel('Dias')
		plt.ylabel('Nº de células cancerígenas')
		plt.savefig('IMAGES/tiempos.png')

#for fichero in directorio: 
#	celulas = [0]
#	tiempos = [0]
#	with open(fichero, "r+") as file1:
#		for line in file1:
#			celulas = np.append(celulas, int(line.split()[0]) )
#			tiempos = np.append(tiempos, int(line.split()[1]) )
#		totalcelulas += np.sum(celulas) 
#		totaltiempos += np.sum(tiempos)	
#		meancelulas = np.append(meancelulas, np.sum(celulas)/24)
#		meantiempos = np.append(meantiempos, np.sum(tiempos)/24)
#
#print(totalcelulas)
#print(totaltiempos)

#plt.plot(meancelulas, meantiempos, 'ko')
#plt.title('Mean of time for each day')
#plt.xlabel('Nº Tumor Cells')
#plt.ylabel('Time (ms)')
#plt.savefig('IMAGES/tiempos' + fichero[13:-7] + '.png')


#directorio = glob("TABLES/*")

#for fichero in directorio: 
#	with open(fichero, "r+") as file1:
#		long = file1.readline()
#		nueva_rejilla = np.zeros(( int(long.split()[0]), int(long.split()[0]) ))
#		for line in file1:
#			nueva_rejilla[ int(line.split()[0]) ][ int(line.split()[1]) ] = (5 + int(line.split()[2]))
			

#	plt.imshow(nueva_rejilla, vmin = 0, vmax = 75, cmap="viridis")
#	plt.title('Dia '+fichero[7:-7])
#	plt.xlabel('Eje X')
#	plt.ylabel('Eje Y')
#	plt.savefig('IMAGES'+fichero[6:-3]+'png')
	
#	plt.imshow(nueva_rejilla, vmin = 0, vmax = 1, cmap="Greys")
#	plt.title('Dia '+fichero[7:-7])
#	plt.xlabel('Eje X')
#	plt.ylabel('Eje Y')
#	plt.savefig('IMAGES'+fichero[6:-4]+'byw.png')

