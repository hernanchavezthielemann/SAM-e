#!/usr/bin/python
#h
#h		
#hv			packmoler  version 1.0 (25 Spt 2017)
#h	----------- Managing to create a neat data file  ----------------
#h
#h	By Hernan Chavez Thielemann
#h	herchavezt@gmail.com
#h	
#
#	Ver		 Filename	  Date			Coments
#	---		----------	----------	----------------------------------------
#	1.0 	packmoler.py	26/09/2017	Child of mix2data in order to neat
#h
import os
from sys import exit
from ff_scout import bond_radius

#------------------------------------------------------
#///////	Class definitions are here	///////
#------------------------------------------------- -----

class Packmol(object):
	''' My homemade Packmol wrapper, clearly under development too '''
	
	def __init__(self, pkname="./.lib/packmol"):
		self._machine=pkname
	@property
	def machine(self):
		return self._machine
	@machine.setter
	def machine(self, value):
		self._machine = value	
	@property
	def command(self):
		self._command = self.machine+" < "
		return self._command
	
	def run_listrings(self,listofstrings):
		'''datafile maker'''
		plaintext=""
		if 	type(listofstrings)==type(["mo","del"]):
			for strings in listofstrings:
				plaintext=plaintext+"\n"+strings
		else: 
			raise ValueError("You are making some wrong code...")
		return self.run_text(plaintext)
			
	def run_text(self,plaintext,crfile="aux.inp"):
		'''datafile maker'''
		from moler import text2file
		
		text2file(crfile, plaintext)
		return self.run(self.command+crfile)
	
	def run(self,command):
		'''So self explanatory'''
		from subprocess import Popen, PIPE
		
		proc = Popen(command,stdout=PIPE,shell=True)
		print " Now Packmol is running please wait...\n",
		linesc=[]
		prev_line=''
		vec=['/', '|', u'\u00A6', '\\', '-', u'\u2013']
		
		i=0
		while True:
			print ("\033[0A\033[40C"+vec[i])
			i+=1
			if i==len(vec):
				i=0
			line = proc.stdout.readline()
			l=line.lstrip(' ')
			if linesc<>[]:
				linesc.append(line)
				if l.startswith("Running time:"):
					break
			elif l.startswith('Success') or l.startswith('Final'):
				linesc.append(prev_line)
				linesc.append(line)
			prev_line=line
		
		return linesc

#------------------------------------------------------
#///////	Function definitions are here	///////
#------------------------------------------------------
def mantain(lista,item2mantain='*'):
	for p in range(len(lista)):
		if lista[p][1]<>item2mantain: 
			del lista[p]
			lista=mantain(lista,item2mantain)
			break
	return lista

def remover(lista,item=''):
	for p in range(len(lista)):
		if lista[p]==item: 
			del lista[p]
			lista=remover(lista,item)
			break
	return lista

def pknamebuilder(matrix,filler):
	''' if there are some spare time...'''
	#the idea could be to compose a name as a function of the components

def packmoller(files, f_pars, pack_pars):
	''' makes the optimus geometry with packmol'''
	
	aux1=[]
	aux2=[]
	for p in range(len(f_pars)):
		if f_pars[p][1] is '*':	
			aux1.append([files[p],f_pars[p][0]])# this mean $file-$qty-$crap 
		else:
			aux2.append([files[p]]+f_pars[p][1]) # is filler -> this mean $file-$1-$centerplace ->[x,y,z]
	
	BoxdimX,BoxdimY,BoxdimZ,centerX,centerY,centerZ=pack_pars
	
	# ----------- "Calculus part"
	pebond=0.68
	
	cenx=float(centerX)
	ceny=float(centerY)
	cenz=float(centerZ)
	
	bd_x=float(BoxdimX)
	bd_y=float(BoxdimY)
	bd_z=float(BoxdimZ)	
	
	mindix=str(cenx-bd_x/2+ pebond)
	mindiy=str(ceny-bd_y/2+ pebond)
	mindiz=str(cenz-bd_z/2+ pebond)
	
	maxdix=str(cenx+bd_x/2- pebond)
	maxdiy=str(ceny+bd_y/2- pebond)
	maxdiz=str(cenz+bd_z/2- pebond)
	
	pkname= 'mix_box{:.2f}x{:.2f}x{:.2f}.pdb'.format(bd_x, bd_y, bd_z)
	#------------------ Header
	Header= "\n# By Hernan Chavez Thielemann\nseed 1111111111\ntolerance 1.31"
	Header+= "\n\nfiletype pdb\noutput " + pkname + "\n\n"
	#------------------ FILLERS <aux2>
	fillers_lines=''
	for f in range(len(aux2)):
		fillers_lines += "structure " + aux2[f][0]
		fillers_lines += "\n  number 1\n  center\n"
		la1= cenx - bd_x/2 + aux2[f][1]* bd_x
		la2= ceny - bd_y/2 + aux2[f][2]* bd_y
		la3= cenz - bd_z/2 + aux2[f][3]* bd_z
		fillers_lines+= "  fixed {:.4f} {:.4f} {:.4f} 0. 0. 0.".format(la1,la2,la3)
		fillers_lines+= "\n  radius 1.05\nend structure\n\n"# with this huge radius should be enough for any filler
		
	#------------------ MATRICES <aux1>
	matrix_lines="\n"
	qty_m=0
	for f in range(len(aux1)):
		
		qty_m+=aux1[f][1]
		matrix_lines+="\nstructure "+aux1[f][0]+"\n  number "+str(aux1[f][1])
		matrix_lines+="\n  inside box "
		matrix_lines+=mindix+" "+mindiy+" "+mindiz+" "+maxdix+" "+maxdiy+" "+maxdiz
		
		at_el, at_nums = pdb2atomel(aux1[f][0])# remover(.structure_file2info(aux1[f][0],1),'H')
		for j in range(len(at_el)):
			matrix_lines+= "\n  atoms"+at_nums[j]
			matrix_lines+= "\n    radius "+bond_radius(at_el[j])
		matrix_lines+= "\nend structure\n\n"
	
	#------------------ Footer
	Footer="\n\n# FILLERS number: "+str(len(aux2))
	Footer+="\n# MATRIX number: "+str(qty_m)
	Footer+="\n\nnloop 2000000\n# E basta! Viva Chile!"
	
	plaintext= Header+fillers_lines+"\n"+matrix_lines+"\n"+Footer
	
	dir_path = os.path.dirname(os.path.realpath(__file__))
	PK = Packmol(dir_path+"/.lib/packmol")
	
	print "\n"+"-"*36+"\n Packmol file created as 'test.inp'\n"+"-"*36+"\n"
	
	Pack_terminal = PK.run_text(plaintext,"test.inp")
	
	
	doneflag=False
	if len(fileseeker(dir_path,pkname))>0:
		print '#'*80 
		for x in range(len(Pack_terminal)):
			print Pack_terminal[x] ,
			if not doneflag and Pack_terminal[x].lstrip(' ').startswith('Success'):
				doneflag=True
		print '\n'
		
	return doneflag, pkname

def pdb2atomel(pubchemfile):

	atomic_el=[]	# atom element
	with open(pubchemfile, 'r')  as inpdb:
		for k_line in inpdb:
			if k_line.startswith("ATOM"):
				atomic_el.append(k_line[77])
	laux=[]
	la_nu=[]
	for x in range(len( atomic_el )):
		if atomic_el[x] == 'H':
			pass
		elif atomic_el[x] not in laux:
			laux.append(atomic_el[x])
			la_nu.append(' '+str(x+1))
		else:
			for i in range(len(laux)):
				if atomic_el[x]== laux[i]:
					la_nu[i]+=' '+str(x+1)
					break
	return laux, la_nu

def fileseeker(path= os.getcwd(), word='data'):
	'''seek data & destroy, returns a list of posible files, filtered by "word" criterion'''
	list_of_files=[]
	if path==os.getcwd():
		DIR='.'
	else:
		DIR=path
	for (root, folder, filenames) in os.walk(DIR, topdown=True):
		for name in filenames:
			list_of_files.append(os.path.join(root, name))
	files=[]
	for fs in range(len(list_of_files)):
		# Could implement "while not" and use the same
		#list_of_files with remove(list_of_files[fs])??
		if '/' in list_of_files[fs] and word in list_of_files[fs].split('/')[-1]:
			files.append(list_of_files[fs])
	if files==[]:
		sys.exit("Error!! --- no file(s) found with "+word+" criterion---")
	return files;
