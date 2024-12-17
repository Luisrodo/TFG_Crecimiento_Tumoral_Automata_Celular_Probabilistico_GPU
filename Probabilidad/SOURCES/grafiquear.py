import matplotlib.pyplot as plt
import numpy as np
from glob import glob

directorio = glob("TABLES/*")

for fichero in directorio: 
	with open(fichero, "r+") as file1:
		long = file1.readline()
		nueva_rejilla = np.zeros(( int(long), int(long) ))
		for line in file1:
			nueva_rejilla[ int(line.split()[0]) ][ int(line.split()[1]) ] = (1 + int(line.split()[2]))+1
			

	plt.imshow(nueva_rejilla, vmin = 0, vmax = 30, cmap="viridis")
	plt.title('Dia '+fichero[7:-7])
	plt.xlabel('Eje X')
	plt.ylabel('Eje Y')
	plt.savefig('IMAGES'+fichero[6:-3]+'png')
	
	plt.imshow(nueva_rejilla, vmin = 0, vmax = 1, cmap="Greys")
	plt.title('Dia '+fichero[7:-7])
	plt.xlabel('Eje X')
	plt.ylabel('Eje Y')
	plt.savefig('IMAGES'+fichero[6:-4]+'byw.png')
