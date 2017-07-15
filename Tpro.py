#!/usr/bin/python
#
#	-------------	Gathering and plotting temperature of chunks version 1.1	----------------
#
#	By Hernan Chavez Thielemann
#	herchavezt@gmail.com
#	
#	latest version must be in in the same folder of postproc.py
#
#	ver	1.1	Tpro.py		14 03 2017
#	ver	1.0	Tpro.py		10 03 2017
#
#
#	Options how to call this utility:			## WIP
#
#		(1)>	
#		(2)>	
#		(3)>	
#
#	Note:	if the file that you want to analyze is upstream from your current working directory (cwd)
#		just use the option (3)	putting the hole path to your data file.

import os
import sys
import matplotlib.pyplot as plt
import numpy

#------------------------------------------------------
#///////	Function definitions are here	///////
#------------------------------------------------------
def fileseeker(path=os.getcwd()):
	'''seek & destroy'''
	files=[]
	if path==os.getcwd(): DIR='.'
	else: DIR=path
	for (root, folder, filenames) in os.walk(DIR, topdown=True):
		 for name in filenames:
      			files.append(os.path.join(root, name))	
	return files;


def whichfile(list_of_files,DoR='D',filetype='NE'):
	'''pic a list of files and returns a selected file'''
	if list_of_files==[]: 
		sys.exit("Error!! ---- no files in list ---")
		
	elif list_of_files!=[]:
		print("\nThese are the files that you could analyze:")
		if DoR=='D':
			word='data'
			word2=word
		elif DoR=='R':
			word='lgct'
			word2=word
			if filetype<>'NE':word2=filetype
		elif DoR=='T':
			word='temp'
			word2=word
		files2print=[]
		index=0
		for fs in range(len(list_of_files)):
			filel=list_of_files[fs].split('.', 1)
			# to pick something else than crap
			if len(filel)>1 and word in list_of_files[fs] and word2 in list_of_files[fs] :
				index+=1
				files2print.append(list_of_files[fs])
				print("		# {}.- {}").format(index, list_of_files[fs])

	# once showed - ask which one is needed to analyze
	Selfile= raw_input("Which file do you want to analyze? (number of) #")
	if Selfile=='':
		if 'lgct.txt' in files2print: Filename='lgct.txt'
		else: Filename=files2print[-1]
	else: Filename=files2print[int(Selfile)-1]

	return Filename;

def kappabytemp(tempfile2use,timestep=0.1,dQ=1.6089-16,A=70.5,lz2A=1):
	''' self explanatory '''

	temp_v=[]
	pos=[]

	with open(File2analyze, 'r')  as indata:
		chunknum=''
		firstdump=True
		index=0
		index2=0
		laststep=0
		for k_line in indata:
		   if not k_line.startswith("#"):	
			line=k_line.split(' ')
			while  line[0]=='':line.remove('');
			if chunknum=='': 
				firststep=int(line[0])
				chunknum=int(line[1])
			elif len(line)>3 and firstdump:
				temp_v.append(float(line[3]))
				pos.append(float(line[1])*lz2A)#
				if chunknum==int(line[0]): firstdump=False
			if len(line)>3 and not firstdump:
				temp_v[index]=temp_v[index]+float(line[3])
				index+=1
				if index>=chunknum: index=0; index2+=1;

			elif float(line[1])==chunknum:
				Nfreq=int(line[0])-laststep
				laststep=int(line[0])

		Steps=laststep-firststep+Nfreq		
		rango=Steps/Nfreq
	print index2

	## -------------------------------------------------------------
	#	ordering the vector
	
	beforeminflag=True
	minT=min(temp_v)
	for T in range(chunknum):
		if temp_v[T]<>minT and beforeminflag:temp_v.append(temp_v[T]/rango);
		elif temp_v[T]>=minT:
			beforeminflag=False;
			temp_v[T]=temp_v[T]/rango;

	
	aux_v=[]
	i_s=len(temp_v)-int(chunknum)
	while len(temp_v)>i_s:aux_v.append(temp_v[i_s]);i_s+=1;
	if len(aux_v)==int(chunknum): print '-- all seems ok --'
	else: print '\n\n ------something is wrong--------\n\n'	
	temp_v=aux_v
	
	## -------------------------------------------------------------
	#	obtaining slope 'm'
	i_s=0
	temp_a=[]
	pos_a=[]
	
	for T in range(chunknum):
		if T+2<chunknum:
		  m_prev=(temp_v[T]-temp_v[T-1])/(pos[T]-pos[T-1])
		  m=(temp_v[T+1]-temp_v[T])/(pos[T+1]-pos[T])	
		  if m>0:
		    m_next=(temp_v[T+2]-temp_v[T+1])/(pos[T+2]-pos[T+1])
		    coef=2.5
		    if m*2.0/coef<m_next<m*coef  :
			print m_prev,m, m_next
			if temp_a==[]:
				temp_a.append(temp_v[T])
				pos_a.append(pos[T])
			temp_a.append(temp_v[T+1])
			pos_a.append(pos[T+1])
	
	m,b = numpy.polyfit(pos_a, temp_a, 1)
	A2m=0.0000000001
	fs2s=1E-15
	dT_dz=m/A2m
	A=A*A2m*A2m
	dt=timestep*Steps*fs2s
	Kappa=dQ/(dt*A*dT_dz)
	## -------------------------------------------------------------
	#	To sum up
	print ("\n >	Steps:		{}.").format(Steps)
	print (" >	timestep:	{} [fs].").format(timestep)
	print (" >	sim time:	{} [ps].").format(dt/(fs2s*1000))
	print ( "\n >>>>>>>	Average thermal conductivity:  {:.4f} [W/mK] for {:.4f}[K/A].\n").format(Kappa,m)
	print '\ntaking in count :'
	print temp_a
	#	To plot
	temps=[m*x+b for x in pos_a]
	AverageTempTrend = plt.figure()
	l1, = plt.plot(pos, temp_v,'g*', lw=5)
	l2, = plt.plot(pos_a, temps , 'r-.', lw=4)
	plt.xlabel('Box distance ['+u'\u00c5'+']')
	plt.ylabel('Temperature [K]')
	plt.show()
	return Kappa,Steps,m,b,chunknum;
	
	#END kappabytemp


#------------------------------------------------------
#///////////////	 Main	         /////////////
#------------------------------------------------------

if __name__ == "__main__":

#--------------- Sim Setup ------------------------#
		timestep=0.5

		#dQ=1.27993479346527e-15
	#	dQ=2.42342725426779e-16#0.5
		dQ=1.60896505042396e-16
		lx=76.1488809351133
		ly=75.6683852211769
		#ly=lx
		A=lx*ly	

		lz=1
		#lz=100.696
		
#------------- Core to be copy&imported-------------#
		tryother=True
		while tryother: #*(1) main while
			tryother=False
			File2analyze=whichfile(fileseeker(),'T')
			inp= raw_input("\nSpecify the timestep in [fs]> ")
			if inp=='' or inp=='^[[A': timestep=timestep
			else: timestep=float(inp)
			inp= raw_input("\nSpecify the flux (dQ)> ")
			if inp<>'' or inp<>'^[[A': dQ=dQ
			else: dQ=float(inp)
			Kappa,Steps,m,b,FF=kappabytemp(File2analyze,timestep,dQ,A,lz)

			itry= raw_input("\nDo you want to try another file? (yes/no) > ")
			if itry=='' or itry=='yes' or itry=='y' or itry=='si' or itry=='s': tryother=True
			print("")
	
#------------------------------------------------------
#///////////////	The End		  /////////////
#------------------------------------------------------
		print 'the end'



		



