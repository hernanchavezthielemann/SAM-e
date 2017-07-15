#!/usr/bin/python
#
#	-------------	Gathering and plotting TC version 2.0	----------------
#
#	By Hernan Chavez Thielemann
#	herchavezt@gmail.com
#	
#	latest version must be in PostProcessing directory
#
#	ver	2.0	postpLgv_G.py
#	ver	1.2	postpLgv.py
#	ver	1.1	postp.py
#	ver	1.0	pp.py
#
#	Atoms supported:
#			Carbon
#			Nitrogen
#			Hydrogen
#			Oxygen
#
#	Options how to call this utility:
#
#		(1)>	./check.py # And through the command prompt select the file to analyze downstream your cwd.
#		(2)>	python check.py # the same as first.
#		(3)>	python check.py filename_to_analyze # Where filename_to_analyze is the file to analyze.
#
#	Note:	if the file that you want to analyze is upstream from your current working directory (cwd)
#		just use the option (3)	putting the hole path to your data file.
#
#	----------------------------------------------------------

import os, sys
from os import walk

def taglia(f, n):
    '''Taglia il float f fino a n decimali, senza arrotondare, come round("string")'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])

def Precent2plot():
    		'''pide el % de steps a plotear'''
		percentoplot = raw_input("Percentage to be ploted (0-100, def 100) > ")
		if percentoplot=='':
			percentoplot=100
			print("\033[1A\033[42C 100")
		else: percentoplot=int(percentoplot)
		if 100>percentoplot>9: print("\033[1A\033[45C%")
		elif percentoplot<10: print("\033[1A\033[44C%")
		elif (99.9<percentoplot<100.3): print("\033[1A\033[46C%")
		else: 
			plot(">>	You dreamer, pick a number from 0 to 100 ")
			Precent2plot()
		return percentoplot

def filename2info(Filename):
	'''extract info from filename'''
	infofilename=Filename.split('_')
	if len(infofilename)>1:
		timestep=float(infofilename[1])
		lx,ly,lz=infofilename[3].split('x',2)
	if Filename=='lgct.txt':	
		timestep=float(raw_input("Timestep used (default 0.1)") or '0.1')
		lx=float(raw_input("lx used (default 83.2)") or '83.2')
		ly=float(raw_input("ly used (default 83.2)") or '83.2')
		lz=float(raw_input("lz used (default 83.2)") or '83.2')
	return (timestep,float(lx),float(ly),float(lz))

def fileseeker(path=os.getcwd()):
	'''seek & destroy'''
	files=[]
	if path==os.getcwd(): DIR='.'
	else: DIR=path
	for (root, folder, filenames) in os.walk(DIR, topdown=True):
		 for name in filenames:
      			files.append(os.path.join(root, name))	
	return files;


def wichfile(list_of_files,DoR='l'):
	'''pic a list of files and returns a selected file'''
	if list_of_files==[]: 
		sys.exit("Error!! ---- no files in list ---")
		
	elif list_of_files!=[]:
		print("\nThese are the files that you could analyze:")
		if DoR=='D': word='data'
		else: word='lgct'
		files2print=[]
		index=0
		for fs in range(len(list_of_files)):
			filel=list_of_files[fs].split('.', 1)
			# to pick something else than crap
			if len(filel)>1 and word in list_of_files[fs]:
				index+=1
				files2print.append(list_of_files[fs])
				print("		# {}.- {}").format(index, list_of_files[fs])

	# once showed - ask wich one is needed to analyze
	Selfile= raw_input("Which file do you want to analyze? (number of) #")
	if Selfile=='':
		if 'lgct.txt' in files2print: Filename='lgct.txt'
		else: Filename=files2print[-1]
	else: Filename=files2print[int(Selfile)-1]

	return Filename;

if __name__ == "__main__":
  

  tryother=True
  while tryother: #*(1) main while
	
	TCs=[]
	KRs=[]
	Steps=0
	StepsNumi=0
	StepsNumf=0


	if len(sys.argv)>1:File2analyze=sys.argv[1]
	else: File2analyze=wichfile(fileseeker())

	with open(File2analyze, 'r')  as indT:
		next(indT)
		for k_line in indT:
			line=k_line.split(';', 2)
			if StepsNumi==0: StepsNumi=line[0]
			StepsNumf=line[0]
	
	StepsNum=int(StepsNumf)-int(StepsNumi)


	with open(File2analyze, 'r')  as indT:
		

	#	timestep=float(raw_input("Timestep used (default 0.1)") or '0.1')

		timestep,lx,ly,lz=filename2info(File2analyze)
		SumT=0.0
		k=0
		i=0
		dQ=0.0
		kb=1.38064852E-023
		kcpm2J=4184/6.022E+023
		A2m=0.0000000001
		fs2s=1E-015
		s2ps=1E012

		dz=lz*A2m/2.0
		A=lx*A2m*ly*A2m
		l=72.532

		
		next(indT) 
		for k_line in indT:
			
			line=k_line.split(';', 2)
#			print (line)
			SumT=SumT+float(line[2])
			dQ=float(line[1])
			
			Steps=k*1000
			vect=range(10000,StepsNum+10000,10000)
			
			if Steps in vect:
				dT=SumT/k
				dt=timestep*Steps*fs2s
#		print ("dQ %f for %f " % (dQ, dT))
				dQ=dQ#*kcpm2J
				Kappa=dQ*dz/(dt*A*dT)

				Kapitza=(dt*A*2*dT)/dQ
				print ("TC %f for %d steps!" % (Kappa, Steps))
#		print ("TC %f for %d steps... & Kapitza %f" % (Kappa, Steps, Kapitza))
				TCs.append(Kappa)
				KRs.append(Kapitza) 
				i=i+1
			k=k+1
			
#	print (i, SumT, dQ)
	
		print ( "\n >>>>>	end result")	
	
		dT=SumT/k
		Steps=k*1000
		print (" >	Steps:		{}.").format(Steps)
	
	
		print (" >	timestep:	{} [fs].").format(timestep)
	
		dt=timestep*Steps*fs2s
		print (" >	sim time:	{} [ps].").format(dt*s2ps)
	
		dQ=dQ#*kcpm2J

		
		
		Kappa=dQ*dz/(dt*A*dT)
		Kapitza=(dt*A*dT)/dQ
		TCs.append(Kappa)
		KRs.append(Kapitza) 
	#	print ( "\n >>>>>>>	Average thermal conductivity:  {} [W/mK].\n").format(Kappa)
	#	print ( "\n >>>>>>>	Average Kapitza resistance:  {} [Km^2/W].\n").format(Kapitza)
	j=2
	avTC=0
	numtoaverage=10
	nt=numtoaverage+j
	while j < nt:
		avTC=avTC+TCs[len(TCs)-j]
		j=j+1
	avTC=avTC/numtoaverage
	print ( "\n >>>>>>>	Average thermal conductivity:  {} [W/mK].\n").format(avTC)
	
	
#------------------------------------------------------
#///////////////	Plot section	  /////////////
#------------------------------------------------------
	
	import matplotlib.pyplot as plt
	import numpy as np
	import sys



	paso=10000*timestep*fs2s*s2ps
	tmax=(Steps+10000)*timestep*fs2s*s2ps
	tmin=paso

	t1=np.arange(tmin,tmax,paso)
	# t1=range(tmin,tmax,paso) / deprecated by np.arange

	iplom= raw_input("How many plots do you want to try? (def 1) >")
	if iplom=='':
		iplom=1
		print("\033[1A\033[44C1")
	else: iplom=int(iplom)
	iplo=0
	plt.ion()
	while iplom>iplo:

		print("Graph {} of {}").format(iplo+1, iplom)
		percentoplot=float(Precent2plot())/100
		fig1 = plt.figure(iplo+1)
		
		t1aux=[]
		TCsaux=[]
		for i in range(int(len(TCs)*percentoplot)):
			t1aux.append(t1[i])
			TCsaux.append(TCs[i])
	
		TC=[avTC for i in range(len(t1aux))]
	
#                                                                          Task!
# make it variable			---- im on that
		adjust=round((max(TCs)-min(TCs))/18,3)
		v16=[round(avTC,2)-2*adjust for i in range(len(t1aux))]
		v17=[round(avTC,2)-adjust for i in range(len(t1aux))]
		v18=[round(avTC,2)+adjust for i in range(len(t1aux))]
		v19=[round(avTC,2)+2*adjust for i in range(len(t1aux))]
	
		l1, = plt.plot(t1aux, TCsaux,'-g', lw=2.5)
		l2, = plt.plot(t1aux, TC, 'r-.', lw=1)
		#l3, = plt.plot(t1aux, v16,'-k', lw=0.5)
		#l4, = plt.plot(t1aux, v17,'-k', lw=0.5)
		#l5, = plt.plot(t1aux, v18,'-k', lw=0.5)
		#l6, = plt.plot(t1aux, v19,'-k', lw=0.5)
		plt.grid(which='both',linestyle='-')
		plt.xlabel('time [ps]')
		plt.ylabel('Thermal conductivity [W/mK]')
		plt.title('Thermal conductivity per picosecond ('+str(percentoplot*100)+'% of steps)')
		iplo+=1

		plt.show()
#------------------------------------------------------
#////////////	   want to try another?   /////////////
#------------------------------------------------------
	plt.show()	
	tryother=False
	itry= raw_input("\nDo you want to try another file? (yes/no) > ")
	if itry=='' or itry=='yes' or itry=='y' or itry=='Yes' or itry=='Y' or itry=='YES' or itry=='si' or itry=='s' or itry=='Si' or itry=='SI' or itry=='S': tryother=True
	
#------------------------------------------------------
#///////////////	The End		  /////////////
#------------------------------------------------------

  print("\n		 Seems like all ends ok!\n")
