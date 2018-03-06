#!/usr/bin/python
#h
#h        <    Carbon handler        >
#h
#hv            
__version__ = 'ch version 1.1 (06 Feb 2018)'.split()[2]
#h    -------------    Sorting , checking and crosslinking     ----------------
#h
#h    By Hernan Chavez Thielemann
__author__ = 'Hernan Chavez Thielemann <hchavezthiele at gmail dot com>'
#h    
#h    latest version must be in xlink directory
#
#    Ver    Filename   Date            Coments
#    ---    --------  ----------    ----------------------------------------
#    1.1     ch        06/02/2018
#    1.0     ch        20/06/2017    
#
__help__ = '''
    Atoms supported:
            Carbon C
            
    Options of this utility from command line:
    
        -h      This option some day will give you a helpfull usage overview. 
        
        -x1     To build a carbon nano tube, the user will be prompted to insert 
                the needed parameters  
                
        -x2     To build a graphene sheet , the user will be prompted to insert 
                the needed parameters  
'''

#h        End Help
#------------------------------------------------------
#///    Packages and globals definitions are here   ///
#------------------------------------------------------
import os
import sys    
import math

#------------------------------------------------------
#///////    Class definitions coming up            ///////
#------------------------------------------------------
def cntfixer(coords):
    ''' do nothing by some time till I will fix it '''
    h=bndlengh/2


def cntcreator(m=5,lenght=150,center=[0,0,0],bndlengh=1.418,atelem='C'):
    ''' so self explanatory at the moment just armchair'''
    p=   (3)**(0.5)*0.5*bndlengh # cos
    h=   bndlengh/2              # sin
    phi= 2*math.pi/(3*m) ## radians
    
    cynuml= (lenght//p)
    cynumh= cynuml+1
    if (lenght-cynuml*p)<(cynumh*p-lenght):
        cynum=int(cynuml)
    else:
        cynum=int(cynumh)
        
    r_i = bndlengh/phi
    H_rad_cnt = bndlengh /(2*math.tan(phi/2))
    
    print r_i, H_rad_cnt
    alpha=math.pi/m # ROTATIONAL ANGLE, MAYBE BEST TO CHANGE IT TO RHO
    theta=alpha*2
    
    Ra = [[math.cos(alpha), -math.sin(alpha)], [math.sin(alpha), math.cos(alpha)]]
    #rotational matrix alpha
    Rt = [[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]]
    #rotational matrix theta
    cos_alpha= math.cos(alpha)
    sin_alpha= math.sin(alpha)
    cos_theta= math.cos(theta)
    sin_theta= math.sin(theta)
    
    atoms=[]
    atoms.append([h,-H_rad_cnt,0])
    atoms.append([-h,-H_rad_cnt,0])
    for y in range(m)[1:]:
        X1=atoms[-2][0]*cos_theta-atoms[-2][1]*sin_theta
        Y1=atoms[-2][0]*sin_theta+atoms[-2][1]*cos_theta
        atoms.append([X1,Y1,0])
        X1=atoms[-2][0]*cos_theta-atoms[-2][1]*sin_theta
        Y1=atoms[-2][0]*sin_theta+atoms[-2][1]*cos_theta
        atoms.append([X1,Y1,0])
    
    for z in range(cynum)[1:]:
        
        X1=atoms[-2][0]*cos_alpha-atoms[-2][1]*sin_alpha
        Y1=atoms[-2][0]*sin_alpha+atoms[-2][1]*cos_alpha
        atoms.append([X1,Y1,z*p])
        
        X1=atoms[-2][0]*cos_alpha-atoms[-2][1]*sin_alpha
        Y1=atoms[-2][0]*sin_alpha+atoms[-2][1]*cos_alpha
        atoms.append([X1,Y1,z*p])
        
        for y in range(m)[1:]:
            
            X1=atoms[-2][0]*cos_theta-atoms[-2][1]*sin_theta
            Y1=atoms[-2][0]*sin_theta+atoms[-2][1]*cos_theta
            atoms.append([X1,Y1,z*p])
            
            X1=atoms[-2][0]*cos_theta-atoms[-2][1]*sin_theta
            Y1=atoms[-2][0]*sin_theta+atoms[-2][1]*cos_theta
            atoms.append([X1,Y1,z*p])
            
            
    atoms, caps_at = firstcap( atoms, m,bndlengh,  H_rad_cnt)
    atoms = endcap( atoms, m,bndlengh,  caps_at)
        
        
    txtpdb="HEADER\nTITLE    Built with cnt handler"
    txtpdb=txtpdb+"\nREMARK   by Hernan Chavez Thielemann <hchavezthiele@gmail.com>"
    txtpdb=txtpdb+"\nREMARK   If this doesn't work, use a high-level tool like VMD to create it"
    txtpdb=txtpdb+"\nREMARK"
    
    Id=0

    ic=0
    sId='0'
    for at in range(len(atoms)):
        
        
        Id=+1
        
        cx='%8.3f'%(atoms[at][0]+center[0])
        cy='%8.3f'%(atoms[at][1]+center[1])
        cz='%8.3f'%(atoms[at][2]+center[2])
        while len(cx)<8:cx=" "+cx
        while len(cy)<8:cy=" "+cy
        while len(cz)<8:cz=" "+cz
            
        sId=str(at+1)
        sId=str(at+1-10000*((at+1)//10000))

        while len(sId)<5:sId=" "+sId
        atyp=atelem    
        while len(atyp)<3:atyp=atyp+' '
            
        txtpdb=txtpdb+"\nATOM  "+sId+'  '+atyp+" CNT     1    "+cx+cy+cz+"  1.00  0.00           "+atelem+'  '
        
    txtpdb=txtpdb+"\nEND"
    
    L=str(int(lenght))
    while len(L)<3:L='0'+L
    
    text2file("CNT"+str(m)*2+L+'.pdb',txtpdb)
def firstcap( atoms, m_n, bndlengh, H_rad):
    
    rho = math.pi/m_n
    theta = rho*2

    cos_theta= math.cos(theta)
    sin_theta= math.sin(theta)
    
    alpha_pen = math.pi/5.0 # alpha pentagon
    
    Pd = bndlengh*(2*math.cos(2*alpha_pen)+1)
    Ri1 = bndlengh/( 2*math.sin( math.atan(
                                math.sin( rho) /(Pd/bndlengh +math.cos(rho)))))
    H_ri1 = bndlengh * (Pd/bndlengh +math.cos(rho))/(2 * math.sin( rho))
    
    print H_rad,  H_ri1
    f_dr_1 = H_rad - H_ri1
    di_1 = ( bndlengh*math.sin(2*alpha_pen))
    print di_1
    fst_p = (di_1**2 - f_dr_1**2)**0.5
    
    cap_atoms =[]
    
    cap_atoms.append([Pd/2,-H_ri1,-fst_p])
    cap_atoms.append([-Pd/2,-H_ri1,-fst_p])
    
    for y in range(m_n)[1:]:
        X1=cap_atoms[-2][0]*cos_theta-cap_atoms[-2][1]*sin_theta
        Y1=cap_atoms[-2][0]*sin_theta+cap_atoms[-2][1]*cos_theta
        cap_atoms.append([X1,Y1,-fst_p])
        X1=cap_atoms[-2][0]*cos_theta-cap_atoms[-2][1]*sin_theta
        Y1=cap_atoms[-2][0]*sin_theta+cap_atoms[-2][1]*cos_theta
        cap_atoms.append([X1,Y1,-fst_p])
    
    gama_1 = math.acos(f_dr_1/di_1)
    Ph = bndlengh*math.tan(2*alpha_pen)/2
    Hr2p = H_rad - Ph * f_dr_1/di_1
    fspp = Ph *math.sin(gama_1)
    
    cap_atoms.append([0,-Hr2p,-fspp])
    for y in range(m_n)[1:]:
        X1=cap_atoms[-1][0]*cos_theta-cap_atoms[-1][1]*sin_theta
        Y1=cap_atoms[-1][0]*sin_theta+cap_atoms[-1][1]*cos_theta
        cap_atoms.append([X1,Y1,-fspp])
    
    
    gama_2 = math.pi/7.5
    Z_a_g2 = fspp+1.418*math.sin(gama_2)
    R_a_g2 = Hr2p-1.418*math.cos(gama_2)
    
    cap_atoms.append([0,-R_a_g2,-Z_a_g2])
    for y in range(m_n)[1:]:
        X1=cap_atoms[-1][0]*cos_theta-cap_atoms[-1][1]*sin_theta
        Y1=cap_atoms[-1][0]*sin_theta+cap_atoms[-1][1]*cos_theta
        cap_atoms.append([X1,Y1,-Z_a_g2])
        
        
    
    
    return list(reversed(cap_atoms)) + atoms , cap_atoms

def endcap( atoms, m_n,bndlengh,  caps_at):
    
    finat = 2*m_n
    Zl = atoms[-1][-1]
    for y in range(len(caps_at)): 
        caps_at[y] = caps_at[y][:-1]+[ Zl - caps_at[y][-1]]
        
    return atoms +caps_at
        
    
    
    
    
def sheetcreator(width=50,lenght=100, center=[0,0,0],bndlengh=1.418,atelem='C'):
    ''' so self explanatory '''
    b= bndlengh
    p= bndlengh* (3)**(0.5)/2   #sin of 60
    h= bndlengh* 1/2            #cos of 60
    
    repnum_l= int('{:.0f}'.format((lenght-2*b)/(3*b)))
    repnum_w= int('{:.0f}'.format((width-p)/(2*p)))
    
    atoms=[]
    atoms.append([  p, 0,   center[2]])     #  *
    atoms.append([  0, h,   center[2]])     #     *
    atoms.append([  0, h+b, center[2]])     #     *
    atoms.append([  p, 2*b, center[2]])     #  *
    
    for uno in range(repnum_l+1):
        for rep in range(repnum_w):# grow in x
            for i in range(4):
                atoms.append([  atoms[-4][0]+2*p, atoms[-4][1],   center[2]])
        for i in range(4): #grow in y 3b once
            atoms.append([atoms[i][0], atoms[-4][1]+3*b,   center[2]])
            
    atoms=atoms[:-4]
    pdb_file(atoms, molname='GSH')

def pdb_file(atoms, atelem='C', molname='CNT', filename='' ):
    
    txtpdb="HEADER\nTITLE    Built with cnt handler"
    txtpdb+="\nREMARK   by Hernan Chavez Thielemann <hchavezthiele@gmail.com>"
    txtpdb+="\nREMARK   If this doesn't work, use a high-level tool like VMD to create it"
    txtpdb+="\nREMARK   "+molname+' kind'
    
    ic=0
    sId='0'
    for at in range(len(atoms)):
        cx='%8.3f'%(atoms[at][0])
        cy='%8.3f'%(atoms[at][1])
        cz='%8.3f'%(atoms[at][2])
        sId=str(at+1)
        sId=str(at+1-10000*((at+1)//10000))
        atyp= atelem   
        while len(atyp)<3:
            atyp+=' '
            
        while len(sId)<5:
            sId=" "+sId
            
        txtpdb+="\nATOM  "+sId+'  '+atyp+" "+molname+"     1    "+cx+cy+cz+"  1.00  0.00           "+atelem+'  '
        
    txtpdb+= "\nEND"
    if filename=='':
        filename= molname+'_'+str(len(atoms))+atelem+'.pdb'
    text2file(filename, txtpdb)
    msg=" File pdb created as "+filename
    print( "\n"+"-"*(len(msg)+1)+"\n"+msg+"\n"+"-"*(len(msg)+1)+"\n")
    
    
    
def fileseeker(path=os.getcwd(),word='data'):
    '''seek data & destroy, returns a list of posible files, filtered by "word" criterion'''
    list_of_files=[]
    if path==os.getcwd(): DIR='.'
    else: DIR=path
    for (root, folder, filenames) in os.walk(DIR, topdown=True):
        for name in filenames:
            list_of_files.append(os.path.join(root, name))
    files=[]
    for fs in range(len(list_of_files)):
        # Could implement "while not" and use the same list_of_files with remove(list_of_files[fs])
        if '/' in list_of_files[fs] and word in list_of_files[fs].split('/')[-1]:
            files.append(list_of_files[fs])
    if files==[]:sys.exit("Error!! --- no file(s) found with "+word+" criterion---")
    return files;

def whichfile(list_of_files,word="."):
    '''pic a list of files, ask and returns a selected file'''
    if list_of_files==[]:sys.exit("Error!! --- no file(s) in list ---")
    elif list_of_files!=[]:
        print("\nThese are the available files:")
        files2print=[]
        index=0
        for fs in range(len(list_of_files)):
            if len(list_of_files[fs].split('.', 1))>1 and word in list_of_files[fs]:
                index+=1
                files2print.append(list_of_files[fs])
                print("        # {}.- {}").format(index, list_of_files[fs])
    # once showed - ask which one is needed to analyze
    Selfile= raw_input("Which file do you want to open? (number of) #")
    if Selfile=='':
        if 'lgct.txt' in files2print: Filename='lgct.txt'# Legacy parameter
        else: Filename=files2print[-1]
    elif int(Selfile)>len(files2print):
        print "\nThe number entered is out of range, please try again.."
        Filename=whichfile(list_of_files)
    else: Filename=files2print[int(Selfile)-1]
    return Filename;
    
    # once showed - ask which one is needed to analyze
    Selfile= raw_input("Which file do you want to analyze? (number of) #")
    if Selfile=='':
        if 'lgct.txt' in files2print: Filename='lgct.txt'
        else: Filename=files2print[-1]
    else: Filename=files2print[int(Selfile)-1]
    return Filename;

def text2file(filename,text):
    '''Make a file called "filename" with the "text" inside'''
    out_tfile = open(filename,"w")
    out_tfile.write (text)
    out_tfile.close()
#------------------------------------------------------
#///////////////    Main          /////////////
#------------------------------------------------------

if __name__ == "__main__":
    
    '''creates or fix a carbon nanotube'''
    print ('')
    
    arguments=sys.argv
    
    if len(sys.argv)>1:
        if sys.argv[1]=='-x1':
            chi= raw_input("Insert the chirality m:")
            dist= raw_input("Insert the lenght l:")
            xco= raw_input("Insert the x coord:")
            yco= raw_input("Insert the y coord:")
            zco= raw_input("Insert the z coord:")
            coord=[float(xco),float(yco),float(zco)]
            cntcreator(int(chi),float(dist),coord,bndlengh=1.418,atelem='C')
            
        elif sys.argv[1]=='-x2':
            wid= raw_input("Insert the width w:")
            dist= raw_input("Insert the lenght l:")
            xco= raw_input("Insert the x coord:")
            yco= raw_input("Insert the y coord:")
            zco= raw_input("Insert the z coord:")
            coord=[float(xco),float(yco),float(zco)]
            sheetcreator(int(wid),float(dist), center=coord)

        elif sys.argv[1]=='-h':
            print __help__
          
