#!/usr/bin/python
#
#hv            
__version__ = 'ch version 1.1 (24 Jan 2018)'.split()[2]
#h    -------------    Sorting , checking and crosslinking     ----------------
#h
#h    By Hernan Chavez Thielemann
__author__ = 'Hernan Chavez Thielemann <hchavezthiele at gmail dot com>'
#
#    ver    1.1        24 01 2018

#    Note:    if the file that you want to analyze is upstream from your current working directory (cwd)
#        just use the option (3)    putting the whole path to your data file.

from os import getcwd, path, walk
from os.path import join
from sys import exit, argv
import matplotlib.pyplot as plt
from numpy import polyfit, pi

#------------------------------------------------------
#///////    Function definitions are here    ///////
#------------------------------------------------------
def fileseeker(path=getcwd()):
    '''seek & destroy'''
    files=[]
    if path==getcwd(): DIR='.'
    else: DIR=path
    for (root, folder, filenames) in walk(DIR, topdown=True):
         for name in filenames:
                  files.append(join(root, name))    
    return files;


def whichfile(list_of_files,DoR='D',filetype='NE'):
    '''pic a list of files and returns a selected file'''
    if list_of_files==[]: 
        exit("Error!! ---- no files in list ---")
        
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
                print("        # {}.- {}").format(index, list_of_files[fs])

    # once showed - ask which one is needed to analyze
    Selfile= raw_input("Which file do you want to analyze? (number of) #")
    if Selfile=='':
        if 'lgct.txt' in files2print: Filename='lgct.txt'
        else: Filename=files2print[-1]
    else: Filename=files2print[int(Selfile)-1]

    return Filename;


def kappabytrend(trendfile2use,timestep,lx,ly,lz):

    TCs=[]
    Steps=0
    with open(trendfile2use, "r") as indT:
        
        SumT=0.0
        k=0
        i=0
        dQ=0.0
        
        kb=1.38064852E-023
        kcpm2J=4184/6.022E+023
        A2m=0.0000000001
        fs2s=1E-15
        #==============================================
        A=lx*A2m*ly*A2m
        #=============================================0
    
    
        next(indT) 
    #    while k<kmax :
        for k_line in indT:
            
            line=k_line.split(';', 2)
            SumT=SumT+float(line[2])
            dQ=float(line[1])/2.
            
            Steps=k*1000
            vect=range(10000,4000000,10000)
            
            if Steps in vect:
                dT=SumT/k
                dt=timestep*Steps*fs2s
                dQ=dQ*kcpm2J
                print lz, A/(A2m*A2m), dT, dQ
                dz=lz*A2m
                
                Kappa=dQ*dz/(dt*A*dT)
                
                print ("TC %f for %d steps..." % (Kappa, Steps))
                TCs.append(Kappa)
                i=i+1
            k=k+1
        
        print ( "\n >>>>>    end result")    
        dT=SumT/k
        Steps=k*1000
        print (" >    Steps:        {}.").format(Steps)
        print (" >    timestep:    {} [fs].").format(timestep)
        dt = timestep*Steps*fs2s
        print (" >    sim time:    {} [ps].").format(dt*10**12)
    
        dQ=dQ*kcpm2J
        dz=lz*A2m
        
        Kappa=dQ*dz/(dt*A*dT)
        TCs.append(Kappa) 
        print ( "\n >>>>>>>    Average thermal conductivity:  {} [W/mK].\n").format(Kappa)
    
    
    
    fig1 = plt.figure()
    t1=range(10000,(Steps+10000),10000)
    #print len(TCs)
    #print len(t1)
    j=2
    avTC=0
    
    
    while j < 22:
        avTC=avTC+TCs[len(TCs)-j]
        j=j+1
    
    avTC=avTC/20
    
    
    print (avTC)
    
    TC=[avTC for i in range(len(t1))]
    v16=[0.16 for i in range(len(t1))]
    v17=[0.17 for i in range(len(t1))]
    v18=[0.18 for i in range(len(t1))]
    v19=[0.19 for i in range(len(t1))]
    
    l1, = plt.plot(t1, TCs,'-g', lw=2.5)
    l2, = plt.plot(t1, TC, 'r-.', lw=1)
    l3, = plt.plot(t1, v16,'-k', lw=0.5)
    l4, = plt.plot(t1, v17,'-k', lw=0.5)
    l5, = plt.plot(t1, v18,'-k', lw=0.5)
    l6, = plt.plot(t1, v19,'-k', lw=0.5)
    plt.grid(which='both',linestyle='-')
    plt.xlabel('Steps')
    plt.ylabel('Thermal conductivity [W/mK]')
    plt.title('Thermal conductivity per time step @0.1fs')
    
    plt.show()
    
    #while j < 10000:
    #    print (j)
    #    j=j+1
    #plt.close()

#------------------------------------------------------
#///////////////     Main             /////////////
#------------------------------------------------------

if __name__ == "__main__":

#--------------- Sim Setup ------------------------#
        timestep = 1
        A= pi*(4.215**2-2.815**2)
        lx = A**0.5
        lz = 95
#------------- Core to be imported-------------#
        tryother=True
        while tryother: #*(1) main while
            tryother=False
            File2analyze=whichfile(fileseeker(),'R')
            inp= raw_input("\nSpecify the timestep in [fs]> ")
            if inp=='' or inp=='^[[A': timestep=timestep
            else: timestep=float(inp)
            print '\nSpecify the length ['+u'\u00c5'+']> ',
            inp= raw_input()
            if inp=='' or inp=='^[[A': lz=lz
            else: lz=float(inp)
            
            print "\nSpecify the side of frontal area ["+u'\u00c5'+']> ',
            inp= raw_input()
            if inp=='' or inp=='^[[A': lx=lx
            else: lx=float(inp)
                
            print "\nSpecify the side of frontal area ["+u'\u00c5'+']> ',
            inp= raw_input()
            if inp=='' or inp=='^[[A': ly=lx
            else: ly=float(inp)
                
            kappabytrend(File2analyze,timestep,lx,ly,lz)

            itry= raw_input("\nDo you want to try another file? (yes/no) > ")
            if itry=='' or itry=='yes' or itry=='y' or itry=='si' or itry=='s': tryother=True
            print("")
    
#------------------------------------------------------
#///////////////    The End          /////////////
#------------------------------------------------------
        print 'Ciao'

