#!/usr/bin/python
#h
#h		
#hv			ff_scout version 1.0 (17 Oct 2017)
#h	----------- force field handler  ----------------
#h
#h	By Hernan Chavez Thielemann
#h	herchavezt@gmail.com
#h	
#
#	Ver		 Filename	  Date			Coments
#	---		----------	----------	----------------------------------------
#	1.0 	ff_scout.py	04/10/2017	new, in order to expand the horizons
#h

from sys import exit

#------------------------------------------------------
#///////	Class definitions are here	///////
#------------------------------------------------- -----
class forcefield(object):
	
	def __init__(self,ff_file):
		self._fff=ff_file
		
#------------------------------------------------------
#////	Globals and warning definitions are here	////
#------------------------------------------------- -----
__forcefield_name__= None

def get_ffn():
	
	if __forcefield_name__ is None:
		raise Warning(" Forcefield has not been set.")
	else:
		return __forcefield_name__
	
def set_ffn(ffname):
	
	global __forcefield_name__
	if __forcefield_name__ is None:
		__forcefield_name__ = ffname
	else:
		raise Warning(" Forcefield has already been set.")
		
def wrg_y(text2color='Warning!'):
	'''  warning scheme'''
	return '\033[92m'+text2color+'\033[0m'

def wrg_r(text2color='Warning!!'):
	'''  warning scheme'''
	return '\033[93m'+text2color+'\033[0m'

#------------------------------------------------------
#///////	Function definitions are here	///////
#------------------------------------------------------

def bond_radius(at_element):
	''' radius dictionaryused by packmoler, why here?
	just to neat all regarding atoms characteristics'''
	r_dict={'C':'1.025','O':'0.915','N':'0.935'}
	rad=''
	try:
		rad=r_dict[at_element]
	except KeyError:
		rad='1.03'
	return rad

def atom_element(atom_number):
	
	atoms=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
		   'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
		   'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
		   'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
		   'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La',
		   'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
		   'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W']
	atom_element=''
	try:
		atom_element= atoms[atom_number-1]
	except IndexError:
		sys.exit("Fatal Error!! Atom number "+value+" not supported yet")
	
	return atom_element

def atom_type( atom):
	''' once we know the topology order in the molecule
	is posible to fit/stereotype that atom.
	First I thought to implement an automated method to assign
	tags based on "#templates compass" but then what happens with
	pcff or cff... this is a part that needs the user intervention'''
	if atom._type=="":
		
		if atom.atelem=='Ar':#---------------------------------------------------Ar
			atom.atype='ar'
		
		elif atom.atelem=='C':#---------------------------------------------------C
			if len(atom.bond)==1:
				atom.atype='c1o'
			elif len(atom.bond)==2:
				a1=len(atom._parent.atoms[int(atom._bondedids[0])-1].bond)
				a2=len(atom._parent.atoms[int(atom._bondedids[1])-1].bond)
				
				if a1+a2==5:# carbon in the end of a carbon nanotube
					atom.atype='c3a'
				elif atom.bond==['O','O'] or atom.bond==['S','S']:
					atom.atype='c2='
				else:
					print "Warning!!", atom.bond
					sys.exit(" weird kind of C2 don't you think")
			
			elif len(atom.bond)==3: # hybridization ---  sp2
				a1=len(atom._parent.atoms[int(atom._bondedids[0])-1].bond)
				a2=len(atom._parent.atoms[int(atom._bondedids[1])-1].bond)
				if atom.bond[0]=='O' and min(a1,a2)==1 :	# Carbonyl carbon 
					atom.atype="c3'"
				elif atom.bond[0]=='N':	
					atom.atype='c3a'
				else:
					atom.atype='c3a'
			
			elif len(atom.bond)==4: # hybridization ---  sp3
				
				if	atom.bond[0]=='O':
					aux=atom._parent.atoms[int(atom._bondedids[0])-1].bond
					if atom.bond==['O','H','H','H']:
						atom.atype=	'c41o' # methanol
					elif atom.bond==['O','C','C','H']:
						atom.atype=	'c4a'# c43o - secondary alcohols
						#if aux==['C','H']:elif aux==['C','C']:atom.atype='c4a'
					else:
						atom.atype=	'c4o'# c41o
					
				elif	atom.bond==['N','C','H','H']:	
						atom.atype=	'c4'#'c4m'
				elif	atom.bond==['N','H','H','H']:	
						atom.atype=	'c4n'#???
				elif	not 'H' in atom.bond[:-1]:
					if	not 'H' in atom.bond[-1]:
						atom.atype=	'c44'
					else:
						atom.atype=	'c43'
				else:
					atom.atype=	'c4'
			else: 
				print "Warning!!", atom.bond
				sys.exit(" weird kind of C don't you think")
		
		elif atom.atelem=='H':#---------------------------------------------------H
			if atom.bond==['N']:	atom.atype='h1n'
			elif atom.bond==['O']:	atom.atype='h1o'
			elif atom.bond==['H']:	atom.atype='h1h'
			else: 
				atom.atype='h1'
		
		elif atom.atelem=='He':#-------------------------------------------------He
				atom.atype='he'
		
		elif atom.atelem=='N':#---------------------------------------------------N
			
			if len(atom.bond)==2:
				if atom.bond==['O','O']:
					atom.atype='n2o'
				else:
					print "Warning!!", atom.bond
					sys.exit(" This nitro is not suported yet")
			elif len(atom.bond)==3:
				if atom.bond[2]=='H':
					atom.atype='n3a'# sp3 nitrogen in amines
				elif not 'H' in atom.bond: # class n33
					atom.atype='n3'
				else:
					atom.atype='n3'
			else:
				print "Warning!!", atom.bond
				atom.atype='n3'
		
		elif atom.atelem=='O':#---------------------------------------------------O
			if atom.bond==['O']:
				atom.atype='o1o'
			elif atom.bond in [['C'],['N'],['S']]:
				atom.atype='o1='
			
			elif atom.bond==['C','C']:
				ai1=len(atom._parent.atoms[int(atom._bondedids[0])-1]._bondedids)
				ai2=len(atom._parent.atoms[int(atom._bondedids[1])-1]._bondedids)
				if min(ai1,ai2)<>4:
					atom.atype='o2'
				else:
					atom.atype='o2e'
			elif atom.bond==['C','H']:
				atom.atype='o2h'
			else:
				atom.atype='o2'
		
		elif atom.atelem=='S':#---------------------------------------------------S
			if atom.bond==['O','O']:
				atom.atype='s2='
			elif atom.bond==['C']:
				atom.atype='s1='
			else:
				atom.atype='s'
				
		else:sys.exit("Fatal Error!! atom kind "+atom.atelem+" not supported yet")
	
	#return atom.atype

def atom_type_legacy( atom_l):
	''' first function _'''
	
	if atom_l._type=="":
		if atom_l.atelem=='N':
			if atom_l._bonded2==['C','C','H']:atom_l._type='N2a'
			elif atom_l._bonded2==['C','H','H']:atom_l._type='N1a'
			elif atom_l._bonded2==['C','C','C']:atom_l._type='N3a'
			else: sys.exit('Warning!! weird N ')
		elif atom_l.atelem=='O':
			if atom_l._bonded2==['C','C']:
				a1=len(atom_l._parent.atoms[int(atom_l._bondedids[0])-1].bond)
				a2=len(atom_l._parent.atoms[int(atom_l._bondedids[1])-1].bond)
				if min(a1,a2)==4:atom_l._type='O2e'
				elif  min(a1,a2)==3:atom_l._type='O2'
				else: sys.exit("Warning!! weird'O ")
		elif atom_l.atelem=='C':
			if len(atom_l._bonded2)==2:
				atom_l._type='C3a'
			elif len(atom_l._bonded2)==3:
				a1=len(atom_l._parent.atoms[int(atom_l._bondedids[0])-1].bond)
				a2=len(atom_l._parent.atoms[int(atom_l._bondedids[1])-1].bond)
				a3=len(atom_l._parent.atoms[int(atom_l._bondedids[2])-1].bond)
				if atom_l._bonded2==['C','C','C']:atom_l._type='C3a'
				elif min(a1,a2,a3)==2:atom_l._type='C3o'
				elif max(a1,a2,a3)==4:atom_l._type='C33'
				elif min(a1,a2,a3)==1:atom_l._type='C3a'
				elif atom_l._bonded2==['N','C','C']:atom_l._type='C3n'
				else:
					print "Warning!!", atom_l._bonded2
					sys.exit(" weird kind of C4 don't you think")
			
			elif atom_l._bonded2==['O','C','H','H']:atom_l._type='C4o'# oxygen conected carbon
			elif atom_l._bonded2==['O','C','C','H']:atom_l._type='C4e'# ether carbon
			elif atom_l._bonded2==['C','C','C','C']:atom_l._type='C44'
			elif atom_l._bonded2==['C','H','H','H']:atom_l._type='C4'
			elif atom_l._bonded2==['C','C','H','H']:atom_l._type='C42'
			elif atom_l._bonded2==['N','C','H','H']:atom_l._type='C4m'
			elif atom_l._bonded2==['N','H','H','H']:atom_l._type='C4n'#???
			elif atom_l._bonded2==['C','C','C','H']:atom_l._type='C43'
				
			else: 
				print "Warning!!", atom_l._bonded2
				sys.exit(" weird kind of C4 don't you think")
				
		elif atom_l.atelem=='H':
			if atom_l._bonded2==['N']:atom_l._type='H1n'
			elif atom_l._bonded2==['C']:
				bnds=atom_l._parent.atoms[int(atom_l._bondedids[0])-1].bond
				if len(bnds)==3:atom_l._type='H1n'
				elif len(bnds)==4:atom_l._type='H1'
				else: sys.exit("Warning!! weird H , really unexpected "+str(len(bnds))) 
			else: sys.exit("Warning!! weird H ") 
		else:sys.exit("Fatal Error!! Atom kind "+atom_l.atelem+" not supported yet")
	return atom_l._type	

def atom_charge(atom):
	''' search the bond increments in the force field file,
	otherwise computes the partial charge with gasteiger-marsili
	or some other method to be implemented'''
	
	''' is needed to implement a good way to handle unknown atoms...
	here i'm assuming that the pair x1 - x2 will be in the b_i list'''
	if atom.charge==None:
		Mol=atom._parent# get_norm_type is a function of parent.... 
		#	maybe can be putted in atom class poluting its beauty 
		# bond increments -  considered under bond tagging
		attys=[Mol.get_norm_type(atom.aid, 'bonds')]
		id_b=atom._bondedids
		for fid_b in id_b:
				attys.append(Mol.get_norm_type(fid_b, 'bonds'))
		#print attys
		
		charge=0	
		bn_i_flag=False
		with open(get_ffn(), 'r')  as inff:
			for k_line in inff:
				if bn_i_flag:
					
					if '#' in k_line and not k_line.startswith('#'):
						k_line=k_line.split('#')[0]
					
					aux1=k_line.split()
					info_aux=(aux1[len(aux1)-4:])
					if info_aux==[] or len(info_aux)<=1:
						break
					elif attys[0]==info_aux[0]:
						
						for x in range(len(attys))[1:]:
							if info_aux[1]==attys[x]:
								#print 'here 1', attys
								charge+=float(info_aux[2])
								if 'h'==attys[x][0]:
									atom._parent.atoms[int(atom._bondedids[x-1])-1].charge=info_aux[3]
					elif  attys[0]==info_aux[1]:
						
						for x in range(len(attys))[1:]:
							if info_aux[0]==attys[x]:
								#print 'here 2', attys
								charge+=float(info_aux[3])
								if 'h'==attys[x][0]:
									atom._parent.atoms[int(atom._bondedids[x-1])-1].charge=info_aux[2]
			
				elif k_line.startswith("#bond_increments"):
					bn_i_flag=True
					inff.next()
		atom.charge=charge

def get_equivalence():
	''' gathering the atom equivalence from forcefield,
	and stores it in an efficient dictionary list'''
	bnd_dct={}
	agl_dct={}
	tor_dct={}
	oop_dct={}
	flag1=False
	flag2=False
	with open(get_ffn(), 'r')  as inff:
		for k_line in inff:
			info_aux=(k_line.split()[2:])
			if flag2 and len(info_aux)<=1:
				break
			elif flag2:
				bnd_dct[info_aux[0]]=info_aux[2]
				agl_dct[info_aux[0]]=info_aux[3]
				tor_dct[info_aux[0]]=info_aux[4]
				oop_dct[info_aux[0]]=info_aux[5]
	
			elif flag1 and k_line.startswith("!---"):
				flag2=True
			elif k_line.startswith("#equivalence"):
				flag1=True
				inff.next()
				
	listofdicts	=	[bnd_dct,agl_dct,tor_dct,oop_dct]
	return listofdicts

def atom_mass(atom):
	''' made with atom instead of atom kind, just for beauty,
	at the end would be easier just search masses by atom type'''
	atty=atom.atype
	if atom.mass==None:
		flag1=False
		flag2=False
		with open(get_ffn(), 'r')  as inff:
			for k_line in inff:
				info_aux=(k_line.split()[2:4])
				if flag2 and len(info_aux)<=1:
					print 'Warning! mass for: {} not found'.format(atty)
					break
				elif flag2:
					if info_aux[0]==atty:
						atom.mass=float(info_aux[1])
						break
				elif flag1 and k_line.startswith("!---"):
					flag2=True
				elif k_line.startswith("#atom_types"):
					flag1=True

def typeatom_mass(atty):
	''' just search masses by atom type'''
	flag1=False
	flag2=False
	mass=0
	with open(get_ffn(), 'r')  as inff:
		for k_line in inff:
			info_aux=(k_line.split()[2:4])
			if flag2 and len(info_aux)<=1:
				print '--- Warning! mass for: {} not found'.format(atty)
				break
			elif flag2:
				if info_aux[0]==atty:
					mass= float(info_aux[1])
					break
			elif flag1 and k_line.startswith("!---"):
				flag2=True
			elif k_line.startswith("#atom_types"):
				flag1=True
	return mass

def get_lj96(atom_type_list):
	''' search all atom type coeffs for non bonded interaction '''
	flag1= False
	flag2= False
	atty_list= atom_type_list[:]
	eps_list= []
	r_list= []
	info= []
	with open(get_ffn(), 'r')  as inff:
		for k_line in inff:
			if flag2:
				if '#' in k_line and not k_line.startswith('#'):
					k_line= k_line.split('#')[0]
				elif len(k_line)<2 or k_line.startswith('#'):
					break
				aux= k_line.split()
				try:
					if aux[-3] in atty_list:
						info.append(aux[-3:])
						atty_list.remove(aux[-3])
				except IndexError: 
					print aux[-3:]
			elif flag1 and k_line.startswith("!---"):
				flag2= True
			elif k_line.startswith("#nonbond(9-6)"):
				flag1= True
	
	info.sort()
	for x in range(len(info)):
		eps_list.append(info[x][2])
		r_list.append(info[x][1])

	return eps_list, r_list

def get_bond_ks(bond_type_list):
	''' search all bond type coeffs '''
	flag1= False
	flag2= False
	bndty_list = bond_type_list[:]
	info= []
	with open(get_ffn(), 'r')  as inff:
		for k_line in inff:
			if flag2:
				if '#' in k_line and not k_line.startswith('#'):
					k_line= k_line.split('#')[0]
				elif len(k_line)<2 or k_line.startswith('#'):
					break
				
				aux= k_line.split()
				bnds1='{}-{}'.format(aux[-6],aux[-5])
				bnds2='{}-{}'.format(aux[-5],aux[-6])
				if bnds1 in bndty_list:
					info.append([bnds1]+aux[-4:])
					bndty_list.remove(bnds1)
				elif bnds2 in bndty_list:
					info.append([bnds2]+aux[-4:])
					bndty_list.remove(bnds2)
			
			elif flag1 and k_line.startswith("!---"):
				flag2= True
			elif k_line.startswith("#quartic_bond"):
				flag1= True
	
	info.sort()
	R0_list, K2_list, K3_list, K4_list= [], [], [], []
	nw_bty_list=[]
	for x in range(len(info)):
		nw_bty_list.append(info[x][0])
		R0_list.append(info[x][1])
		K2_list.append(info[x][2])
		K3_list.append(info[x][3])
		K4_list.append(info[x][4])
		
	for bty in range(len(bond_type_list)):
		if bond_type_list[bty]<>nw_bty_list[bty]:
			print "--- Fatal Error!! bond kind "+bond_type_list[bty]+" not found  ---"
			print  nw_bty_list[bty]
			exit(" Include that coefficients to avoid this exit")
		
	return R0_list, K2_list, K3_list, K4_list

def get_angle_ks(angle_type_list):
	''' search all bond type coeffs '''
	flag1= False
	flag2= False
	angty_list=angle_type_list[:]
	info=[ ['',0,0,0,0]	for x in range(len(angty_list))]
	
	with open(get_ffn(), 'r')  as inff:
		for k_line in inff:
			if flag2:
				if '#' in k_line and not k_line.startswith('#'):
					k_line= k_line.split('#')[0]
				elif len(k_line)<2 or k_line.startswith('#'):
					break
				
				aux= k_line.split()
				angle1='{}-{}-{}'.format(aux[-7], aux[-6], aux[-5])
				angle2='{}-{}-{}'.format(aux[-5], aux[-6], aux[-7])
				if angle1 in angty_list:
					i= angty_list.index( angle1)
					info[i]= [angle1]+ aux[-4:]
				
				elif angle2 in angty_list: # disordered in the force field
					i= angty_list.index( angle2)
					info[i]= [angle2]+ aux[-4:]
			
			elif flag1 and k_line.startswith("!---"):
				flag2= True
			elif k_line.startswith("#quartic_angle"):
				flag1= True
		
	Theta0_list, K2_list, K3_list, K4_list= [], [], [], []
	for x in range(len(info)):
		Theta0_list.append(info[x][1])
		K2_list.append(info[x][2])
		K3_list.append(info[x][3])
		K4_list.append(info[x][4])
		
		if info[x]==['',0,0,0,0]:
			print "--- "+wrg_r()+" angle kind miss --- "
			print angle_type_list[x]+" coeffs filled with 0! \n"
			
	return Theta0_list, K2_list, K3_list, K4_list

def get_angle_bb_Ks(angle_type_list, bond_radius):	#	BondBond Coeffs
	''' search all bond type coeffs '''
	flag1= False
	flag2= False
	angty_list = angle_type_list[:]
	info=[ ['',0,0,0]	for x in range(len(angty_list))]
	with open(get_ffn(), 'r')  as inff:
		for k_line in inff:
			if flag2:
				if '#' in k_line and not k_line.startswith('#'):
					k_line= k_line.split('#')[0]
				elif len(k_line)<2 or k_line.startswith('#'):
					break
					
				aux= k_line.split()
				angle1='{}-{}-{}'.format(aux[-4], aux[-3], aux[-2])
				angle2='{}-{}-{}'.format(aux[-2], aux[-3], aux[-4])
				if angle1 in angty_list: # same order
					
					radds= get_radii([aux[-2], aux[-3], aux[-4]], bond_radius)
					i= angty_list.index( angle1)
					info[i]= [angle1, aux[-1]]+radds
					
				elif angle2 in angty_list: # inverted order
					
					radds= get_radii([aux[-4], aux[-3], aux[-2]], bond_radius)
					i= angty_list.index( angle2)
					info[i]= [angle2, aux[-1]]+[radds[1], radds[0]]
					print info[i], ' - twist warning!'
				
			elif flag1 and k_line.startswith("!---"):
				flag2= True
			elif k_line.startswith("#bond-bond "):
				flag1= True

	Kbb_list, Rb_list, Rbp_list= [], [], []
	for x in range(len(info)):
		Kbb_list.append(info[x][1])
		Rb_list.append(info[x][2])
		Rbp_list.append(info[x][3])
		
		if info[x]==['',0,0,0]:
			print "--- "+wrg_r()+" angle (bond-bond) kind miss --- "
			print angle_type_list[x]+" coeffs filled with 0! \n"
		
	return Kbb_list, Rb_list, Rbp_list

def get_radii(attys, bond_radius_l):
	''' a list of ordered atoms is given to get their radius '''
	bnd_list = [bond_radius_l[x][0] for x in range(len(bond_radius_l))]
	pairs=[]
	for i in range(len(attys)-1):
		pairs.append([attys[0+i],attys[1+i]])
	#print pairs, bond_radius_l
	radiis=[]
	for bn in range(len(pairs)):
		bnds1='{}-{}'.format(pairs[bn][0], pairs[bn][1])
		bnds2='{}-{}'.format(pairs[bn][1], pairs[bn][0])
		if bnds1 in bnd_list:
			radiis.append(bond_radius_l[bnd_list.index(bnds1)][1])
		elif bnds2 in bnd_list:
			radiis.append(bond_radius_l[bnd_list.index(bnds2)][1])
	
	if len(radiis)<>len(pairs):
		print '--- '+wrg_r()+' radius mismatch', radiis
		
	return radiis

def get_angle_ba_Ks(angle_type_list):
	''' search all bond type coeffs '''
	flag1= False
	flag2= False
	angty_list = angle_type_list[:]
	info=[ ['',0,0,0,0]	for x in range(len(angty_list))]
	with open(get_ffn(), 'r')  as inff:
		for k_line in inff:
			if flag2:
				if '#' in k_line and not k_line.startswith('#'):
					k_line= k_line.split('#')[0]
				elif len(k_line)<2 or k_line.startswith('#'):
					break
				
				aux= k_line.split()
				angle1='{}-{}-{}'.format(aux[-5], aux[-4], aux[-3])
				angle2='{}-{}-{}'.format(aux[-3], aux[-4], aux[-5])
				if angle1 in angty_list: # same order
					i= angty_list.index( angle1)
					info[i]= [angle1]+ aux[-2:]
					
				elif angle2 in angty_list: # inverted order
					i= angty_list.index( angle2)
					info[i]= [angle2]+ aux[-2:]
					print info[i], ' - twist warning!'
			
			elif flag1 and k_line.startswith("!---"):
				flag2= True
			elif k_line.startswith("#bond-angle"):
				flag1= True
	
	Kb_theta_list, Kbp_theta_list= [], []
	for x in range(len(info)):
		Kb_theta_list.append(info[x][1])
		Kbp_theta_list.append(info[x][2])
		
		if info[x]==['',0,0,0,0]:
			print "--- "+wrg_r()+" angle(bond-angle) kind miss --- "
			print angle_type_list[x]+" coeffs filled with 0! \n"
			
	return Kb_theta_list, Kbp_theta_list

def get_dihedral_ks(dihe_type_list):
	''' search all bond type coeffs '''
	flag1= False
	flag2= False
	dihe_list=dihe_type_list[:]
	info=[ ['',0,0,0,0,0,0]	for x in range(len(dihe_list))]
	
	with open(get_ffn(), 'r')  as inff:
		for k_line in inff:
			if flag2:
				if '#' in k_line and not k_line.startswith('#'):
					k_line= k_line.split('#')[0]
				elif len(k_line)<2 or k_line.startswith('#'):
					break
				
				aux= k_line.split()
				dihe1='{}-{}-{}-{}'.format(aux[-10], aux[-9], aux[-8], aux[-7])
				dihe2='{}-{}-{}-{}'.format(aux[-7], aux[-8], aux[-9], aux[-10])
				if dihe1 in dihe_list:
					i= dihe_list.index( dihe1)
					info[i]= [dihe1]+ aux[-6:]
				
				elif dihe2 in dihe_list: # disordered in the force field
					i= dihe_list.index( dihe2)
					info[i]= [dihe2]+ aux[-2:]+ aux[-4:-2]+aux[-6:-4]
					
					print wrg_y()+' Torsion ',
					print info[i][0]+' - twist warning!'
					print '> ', info[i][1:]
			
			elif flag1 and k_line.startswith("!---"):
				flag2= True
			elif k_line.startswith("#torsion_3"):
				flag1= True
		
	V1_list, Phi01_list, V2_list, Phi02_list, V3_list, Phi03_list= [], [], [], [], [], []
	for x in range(len(info)):
		V1_list.append(info[x][1])
		Phi01_list.append(info[x][2])
		V2_list.append(info[x][3])	
		Phi02_list.append(info[x][4])
		V3_list.append(info[x][5])	
		Phi03_list.append(info[x][6])
		
		if info[x]==['',0,0,0,0,0,0]:
			print "--- "+wrg_r()+" dihedral kind miss --- "
			print dihe_type_list[x]+" coeffs filled with 0! \n"
			
	return V1_list, Phi01_list, V2_list, Phi02_list, V3_list, Phi03_list

def get_dihedral_mbt_ks(dihe_type_list, bond_radius):
	''' search all #middle_bond-torsion_3 type coeffs '''
	flag1= False
	flag2= False
	dihe_list=dihe_type_list[:]
	info=[ ['',0,0,0,0]	for x in range(len(dihe_list))]
	
	with open(get_ffn(), 'r')  as inff:
		for k_line in inff:
			if flag2:
				if '#' in k_line and not k_line.startswith('#'):
					k_line= k_line.split('#')[0]
				elif len(k_line.split())<2 or k_line.startswith('#'):
					break
				
				aux= k_line.split()
				f1= [aux[-7], aux[-6], aux[-5], aux[-4]]
				f2= [aux[-4], aux[-5], aux[-6], aux[-7]]
				dihe1='{}-{}-{}-{}'.format(*f1)
				dihe2='{}-{}-{}-{}'.format(*f2)
				
				if dihe1 in dihe_list:
					rads= get_radii(f1, bond_radius)
					i= dihe_list.index( dihe1)
					info[i]= [dihe1]+ aux[-3:]+ [rads+[0]]
				
				elif dihe2 in dihe_list: # disordered in the force field
					rads= get_radii(f2, bond_radius)
					i= dihe_list.index( dihe2)
					info[i]= [dihe2]+ aux[-3:]+ [rads+[1]]
					print wrg_y()+' Middle_bond-torsion ',
					print info[i][0]+' - twist warning!'
					print '> ', info[i][1:]
			
			elif flag1 and k_line.startswith("!---"):
				flag2= True
			elif k_line.startswith("#middle_bond-torsion_3"):
				flag1= True
		
	F1_list, F2_list, F3_list, R0_list= [], [], [], []
	for x in range(len(info)):
		F1_list.append(info[x][1])
		F2_list.append(info[x][2])
		F3_list.append(info[x][3])	
		R0_list.append(info[x][4])
		
		if info[x]==['',0,0,0,0]:
			print "--- "+wrg_r()+" dihedral (middle_bond) kind miss --- "
			print dihe_type_list[x]+" coeffs filled with 0! \n"
			
	return F1_list, F2_list, F3_list, R0_list

def get_dihedral_ebt_ks(dihe_type_list, radds_list):
	''' search all #end_bond-torsion_3 type coeffs '''
	flag1= False
	flag2= False
	dihe_list=dihe_type_list[:]
	info=[ ['',0,0,0,0,0,0,0,0]	for x in range(len(dihe_list))]
	
	with open(get_ffn(), 'r')  as inff:
		for k_line in inff:
			if flag2:
				if '#' in k_line and not k_line.startswith('#'):
					k_line= k_line.split('#')[0]
				elif len(k_line.split())<2 or k_line.startswith('#'):
					break
				aux= k_line.split()
				j=0
				m=1
				if len(aux)<12:
					j=3
					m=2
				
				dihe1='{}-{}-{}-{}'.format(aux[-10+j], aux[-9+j], aux[-8+j], aux[-7+j])
				# j really not needed for dihe2... but for completeness
				dihe2='{}-{}-{}-{}'.format(aux[-7+j], aux[-8+j], aux[-9+j], aux[-10+j])
				
				if dihe1 in dihe_list:
					i= dihe_list.index( dihe1)
					rd= [radds_list[i][0], radds_list[i][2]]
					info[i]= [dihe1]+ aux[-6+j:]*m+ rd
				
				elif dihe2 in dihe_list: # disordered in the force... field (drums)
					i= dihe_list.index( dihe2)
					if radds_list[i][3]:
						rd= [radds_list[i][0], radds_list[i][2]]
					else:
						rd= [radds_list[i][2], radds_list[i][0]]
					info[i]= [dihe2]+ aux[-3:]+ aux[-6:-3]+ rd
					print wrg_y()+' End_bond-torsion ',
					print info[i][0]+' - twist warning!'
					print '> ', info[i][1:]
			
			elif flag1 and k_line.startswith("!---"):
				flag2= True
			elif k_line.startswith("#end_bond-torsion_3"):
				flag1= True
		
	F1d_list, F2d_list, F3d_list, R0d_list= [], [], [], []
	F1i_list, F2i_list, F3i_list, R0i_list= [], [], [], []
	for x in range(len(info)):
		try:
			F1d_list.append(info[x][1])
			F2d_list.append(info[x][2])
			F3d_list.append(info[x][3])
			F1i_list.append(info[x][4])
			F2i_list.append(info[x][5])
			F3i_list.append(info[x][6])	
			R0d_list.append(info[x][7])
			R0i_list.append(info[x][8])
		except IndexError:
			print info[x]
		
		if info[x]==['',0,0,0,0,0,0,0,0]:
			print "--- "+wrg_r()+" dihedral (end_bond) kind miss --- "
			print dihe_type_list[x]+" coeffs filled with 0! \n"
	_2return = F1d_list, F2d_list, F3d_list, F1i_list, F2i_list, F3i_list
	_2return += R0d_list, R0i_list
	return _2return
