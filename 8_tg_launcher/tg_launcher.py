#!/usr/bin/python

import os
from sys import exit, argv
from numpy import std, median

class In_scr(object):
    '''to order the in.script, diferents'''
    
    def __init__(self,filein,crad='10',thermo1='1000',tims='1'):
        self._filein=filein
        self._th=thermo1
        self._cor=crad
        self._timestep=tims
        
        self._header=''
        self._body=''
        self.setheader()
        
    def setheader(self):
        Hdr="\nunits real\nboundary p p p\natom_style full"
        Hdr+=("\nbond_style class2"
              +"\nangle_style class2"
              +"\ndihedral_style class2"
              +"\nimproper_style class2"+"\n"
              +"\npair_style lj/class2/coul/long "+self._cor)
              
        if 'data' in self._filein:
            Hdr+="\nread_data "+self._filein+'\n\n'
        elif 'restart' in self._filein:
            Hdr+="\nread_restart "+self._filein+'\n\n'
        else:
            exit('Error!! -- file format -- < just restart allowed >')
            
        Hdr+= ( "\n"+"kspace_style ewald 1e-6"
               +"\n"+"neighbor 1.5 bin"
               +"\n"+"pair_modify tail yes"
               +"\n"+"special_bonds lj/coul 0.0 0.0 0.5"+"\n"
               +"\n"+"thermo "+self._th
               +"\n"+"timestep "+self._timestep+"\n"
               +"\n"+"variable mass_g equal mass(all)/6.02e23"
               +"\n"+"variable vol_cm3 equal vol*1.0e-24"
               +"\n"+"variable mydensity equal v_mass_g/v_vol_cm3"+"\n"+"\n"
               +"\n"+"thermo_style custom step elapsed temp etotal lx ly lz"
                    +" vol density press"+"\n")

        self._header=Hdr
    
    def setbody(self, values):
        
        T_start, T_end = values
        print values
        rat = 200/300.0 #[steps/(1000 K)]
        
        T_end_s = str(T_end)
        # first go to that temp
        bdy = ( "fix 1 all npt temp "+str(T_start)+" "+T_end_s
               +" 100 aniso 1 1 1000"+"\n"
               +"run "+str(int(abs(T_end-T_start)*rat)*1000)+"\n"
               +"unfix 1"+"\n\n"
              )
        # secon stabilize
        bdy += ( "fix 1 all npt temp "+T_end_s+" "+T_end_s
                +" 100 aniso 1 1 1000"+"\n"
                +"run 200000"+"\n"
                +"unfix 1"+"\n\n"
               )
        # third: run to extract data
        bdy += ( "fix 1 all npt temp "+T_end_s+" "+T_end_s
                +" 100 aniso 1 1 1000"+"\n"
                +"fix av_rho all ave/time 1 40 1000 v_mydensity"+"\n"
                +"thermo_style custom step elapsed temp etotal lx ly lz"
                +"vol f_av_rho density press"+"\n"
                +"run 200000"+"\n"
                +"write_restart restart.npt_K"+"\n"
               )
        self._body= bdy
        return int(abs(T_end-T_start)*rat)*1000+400000
    def write_script(self,crfile="in.aux"):
        '''datafile maker'''
        out_dfile = open(crfile,"w")
        out_dfile.write (self._header+self._body)
        out_dfile.close()
        return crfile

class Slurm(object):
    '''to order the in.script, diferents'''
    
    def __init__(self,na,steps=600000,lmsname="lmp",cors=1):
        self._machine=lmsname
        self._cor=str(cors)
        self._simname=na
        self._steps=steps
        
        self._script=''
        self.set_script()
        
    def set_script(self):
        
        rata = 210/200000.0 # minutos / steps
        time = int(self._steps*rata) # minutos
        horas = int(time/60.0)
        minutos = time - horas*60
        
        Sc=(
            "\n"+"#!/bin/bash"+
            "\n"+"#SBATCH --job-name=jn_"+self._simname+"_tg"+
            "\n"+"#SBATCH --partition=global"+
            "\n"+"#SBATCH --time="+str(horas)+":"+str(minutos)+":00"+
            "\n"+"#SBATCH --nodes=1"+
            "\n"+"#SBATCH --ntasks-per-node="+self._cor+
            "\n"+"#SBATCH --mem-per-cpu=1300M"+"\n"+
            "\n"+"module purge"+
            "\n"+"module load lammps/11Aug17"+
            "\n"+"cd $SLURM_SUBMIT_DIR"+
            "\n"+'CASE_IN="in.Tg"'+
            "\n"+'CASE_OUT="out.sortofoutput"'+"\n"+
            "\n"+"mpirun -np "+self._cor+' '+self._machine
                +" -log none -echo screen -in $CASE_IN > $CASE_OUT \n\n"
            )
        self._script=Sc
    @property
    def lammps(self):
        return self._machine
    @lammps.setter
    def lammps(self, value):
        self._machine = value
        
    def write_sbatch(self,crfile="j.sbatch"):
        '''datafile maker'''
        out_dfile = open(crfile,"w")
        out_dfile.write (self._script)
        out_dfile.close()
        return crfile
    
#------------------------------------------------------
#///////    Function definitions coming up      ///////
#------------------------------------------------------

def make_tg_run(T_h, T_l, T_point_num, relaxed_file, name='0xl',lms='lammps'):
    ''' the core part that organizes everything'''
    
    distance = T_h-T_l
    step = distance/((T_point_num-1)*1.0)
    listof_temps=[T_h]
    for x in range(T_point_num-1):
        listof_temps.append(round(listof_temps[-1]-step, 2))
    if verbose:
        print len(listof_temps)
        print listof_temps
        
    try:
        os.makedirs("tg")
    except OSError:
        exit('Woops!')
    for sn in range(len(listof_temps)):
        os.makedirs(os.path.join("tg", str(int(listof_temps[sn]))))
        
    #lms = Lammps(lmsname = "lmp_mpi", cpus="32")
    
    in_rth = In_scr('../.'+relaxed_file) #lms,
    tim_st = in_rth.setbody([300, T_h])
    temp_s = str(int(T_h))
    in_rth.write_script('./tg/'+temp_s+'/in.Tg')
    in_slu = Slurm(temp_s, tim_st, lms, 1)
    in_slu.write_sbatch('./tg/'+temp_s+'/j.sbatch')
    
    for sn in range(T_point_num)[1:]:
        in_rth = In_scr('../'+str(int(listof_temps[sn-1]))+'/restart.npt_K') #lms,
        tim_st = in_rth.setbody([listof_temps[sn-1], listof_temps[sn]])
        temp_s = str(int(listof_temps[sn]))
        in_rth.write_script('./tg/'+temp_s+'/in.Tg')
        in_slu = Slurm(temp_s, tim_st, lms, 1)
        in_slu.write_sbatch('./tg/'+temp_s+'/j.sbatch')
    
    
def sort_npt_filelist(_filelist):
    tosort=[]
    for n_f in range(len(_filelist)):
        tosort.append([_filelist[n_f].split('restart.npt')[-1].rstrip('ps'),
                       _filelist[n_f]])
    tosort.sort()
    for n_f in range(len(_filelist)):
        _filelist[n_f] = tosort[n_f][1]
    return _filelist

def fileseeker(path=os.getcwd(),word='data'):
    '''seek data & destroy, returns a list of posible files,
    filtered by "word" criterion'''
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
        # Could implement "while not" and use the same list_of_files with remove(list_of_files[fs])
        if '/' in list_of_files[fs] and word in list_of_files[fs].split('/')[-1]:
            files.append(list_of_files[fs])
    if files==[]:
        exit("Error!! --- no file(s) found with "+word+" criterion---")
    return files;
#------------------------------------------------------
#///////////////    Main          /////////////
#------------------------------------------------------

if __name__ == "__main__":
    
    arguments=argv
    help_args= ''' HELP\n python tg_launcher.py <+args>\n
                  arguments: [tg_launcher.py] [1] [2] [3] [4]
                                
                                [1] : T_high
                                [2] : T_low (default = 273 K)
                                [3] : T_points (default = 20)
                                [3] : verbose (default = 0)
                Example:
                
                        >python tg_launcher.py 600
                        
                        >will output Tg for the folowing temperatures:
                        > 600.0; 582.8; 565.6; 548.4; 531.2;
                        > 514.0; 496.7; 479.5; 462.3; 445.1;
                        > 427.9; 410.7; 393.5; 376.3; 359.1;
                        > 341.9; 324.6; 307.4; 290.2; 273.0
                '''
    T_low = 273
    T_points = 20
    verbose = 1
    not_help_flag = True
    if len(arguments)<2:
        exit('Error!! -- minimum input needed -- ')
    elif len(arguments)>=2:
        if arguments[1]=='-h':
            print help_args
            not_help_flag = False
        else:
            T_high = int(arguments[1])
    if not_help_flag:
        if len(arguments)>2:
            T_low = int(arguments[2])
        if len(arguments)>3:
            T_points = int(arguments[3])
        if len(arguments)>4:
            verbose = int(arguments[4])
            
        listoffiles = fileseeker(word='restart.npt')
        relaxed_file = sort_npt_filelist(listoffiles)[-1]
        
        if verbose:
            print listoffiles
            print relaxed_file
            
        make_tg_run(T_high, T_low, T_points, relaxed_file, '0xl')
    
    
    
