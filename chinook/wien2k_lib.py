# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 11:31:20 2022

@author: Andrew
"""
import numpy as np

# import other chinook files to help create correct data structures
import chinook.klib as klib
import chinook.TB_lib as TBlib
import chinook.orbital as orblib
import chinook.electron_configs as eclib
import chinook.build_lib as buildlib


# setting up fundamental constants
uBohr = 5.29e-11
uAngs = 1e-10
Ry_to_eV = 13.606 # 1Ry = 13.606eV


def create_kobject(case):
    '''
    Function to create kobjects using the klib library.
    '''
    kpoints,kpoints_linear = get_kpoints_from_spaghettiene(case)
    labels,kpoints_brk_line = get_labels_from_klist(case)
    
    kobj = klib.kpath(kpoints)#need to pass in points
    
    # updates the kobject now that we've instantiated it with our points
    kobj.kcut = kpoints_linear
    kobj.kcut_brk = kpoints_linear[kpoints_brk_line]
    kobj.labels = labels
    
    return kobj
    
    
    


def get_kpoints_from_spaghettiene(case):
    '''
    Function to extract kpoints from case.spaghetti_ene file. Does this by
    going through the file until it reaches the second band, whilst extracting
    the kvalues
    '''
    
    kpoints = []
    kpoints_linear = []
    
    print('get_kpoints_from_spaghettiene: Opening ' + case + '.spaghetti_ene')
    
    f = open(case + '.spaghetti_ene','r')
    
    nband = 0
    for line in f: # will go through the first set of lines initially, getting the k values
        if 'bandindex' in line:
            nband += 1
            if nband > 1: # if we've been through the first band entirely, break from the loop
                break
        else:
            vals = [v for v in line.split(' ') if v] # splits line, removes whitespace
            kx = float(vals[0])
            ky = float(vals[1])
            kz = float(vals[2])
            
            kpoints.append([kx,ky,kz]) # are already given in absolute units!
            kpoints_linear.append(float(vals[3]))
    
    f.close()
    print('get_kpoints_from_spaghettiene: Closing {}.spaghetti_ene'.format(case))

    # formats the arrays into numpy arrays
    kpoints = np.array(kpoints)
    kpoints_linear = np.array(kpoints_linear)

    return kpoints,kpoints_linear




def get_bands_from_spaghettiene(case):
    '''
    Function to extract bands from case.spaghetti_ene file and return them in
    the appropriate unmpy array format
    '''
    
    bands = []
    numbands = 0

    print('get_bands_from_spaghettiene: Opening {}.spaghetti_ene'.format(case))
    f = open('{}.spaghetti_ene'.format(case),'r') # opens the file for reading
    
    bandE = []
    for line in f: # loops through all the lines in the file
    
        if 'bandindex' in line or not line: # implies we've reached a new band, or end of file
            numbands += 1
            # need to append current band data onto the bands data structure
            if bandE != []: # as long as the current band isn't empty (i.e. upon initialisation)
                bands.append(bandE)
                bandE = []
                
        else: # if the line simply has numbers, etc
            energies = [e for e in line.split(' ') if e] # removes all empty strings from list, incase line has weird spacings
            bandE.append(float(energies[4])) # returning the float strips the trailing \n
            
    f.close() # closes the spaghettiene file to avoid file corruption
    print('get_bands_from_spaghettiene: Closing {}.spaghetti_ene'.format(case))
    
    # we add 1 to the number of bands at the end of the file, so we subtract it here
    numbands = numbands - 1
    
    # converts bands, kpoints, to numpy arrays
    bands = np.transpose(np.array(bands))
    
    return bands,numbands




def get_labels_from_klist(case):
    '''
    Function to get the labels and cut points from the case.klist_band file.
    Returns the text labels and the linenumbers of the breaks (changes in direction)
    '''
    
    labels = []
    kpoints_brk_line = []
    numline = -1
    
    print('get_labels_from_klist: Opening {}.klist_band'.format(case))
    f = open('{}.klist_band'.format(case),'r')
    
    for line in f:
        numline += 1
        if line[0] != ' ': #if the line starts with characters, indicating a special point
            line = [e for e in line.split(' ') if e] # splits the line and removes whitespace
            if 'END' not in line[0]:
                labels.append(line[0])
                kpoints_brk_line.append(numline)
    
    f.close()
    print('get_labels_from_klist: Closing {}.klist_band'.format(case))

    # adjust labels so they display like Latex labels
    for i,lab in enumerate(labels):
        if lab == 'GAMMA': # include any other special points
            lab = '\\Gamma'
        labels[i] = '${}$'.format(lab)
    
    return labels,kpoints_brk_line



def change_kobject_TB(case,kobj,TB):
    '''
    A function to change a TB object and a kobject to allow us to load in our 
    own data generated by WIEN2K. 
    '''
    # uses the above functions to load in the kpoints, labels, etc
    kpoints,kpoints_linear = get_kpoints_from_spaghettiene(case)
    labels,kpoints_brk_line = get_labels_from_klist(case)
    
    
    # updates the kobject
    kobj.kpts = kpoints
    kobj.kcut = kpoints_linear
    kobj.kcut_brk = kpoints_linear[kpoints_brk_line]
    kobj.labels = labels
    
    TB.Kobj = kobj
    
    # uses the functions above to load in the bands and put them into the TB object
    bands,numbands = get_bands_from_spaghettiene(case)
    TB.Eband = bands  
    # need to change the TB.basis object to have a length of numbands
    #print(TB.basis)
    TB.basis = [ii for ii in range(numbands)]
    
    return kobj,TB
        
    
def get_PLV_from_struct(case):
    '''
    A function to extract the PLVs from the case.struct file.
    The function assumes the PLVs are given on the 4th line of the file.
    '''
    # load in the struct file
    f = open('{}.struct'.format(case),'r')
    fcont = f.read().split('\n') # load the files contents as lines into fcont
    f.close()
    print('get_PLV_from_struct: Data loaded from {}.struct'.format(case))
    
    # the third line contains info on the mode of calculation: relative or absolute coords, and lattice vectors in bohr or angstrom
    line = fcont[2]
    
    try:
        elem = [s for s in line.split(' ') if s] # split line and remove whitespace
        # calc_mode not needed if we're just extracting the PLVS
        #calc_mode = elem[2].split('=')[1] # gets the calc mode, RELA or (ABS?)
        unit = elem[3].split('=')[1] # gets the unit type, bohr or (ang?)
    except: # by default, if file structure is different, go for defaults
        print('get_PLV_from_struct: except triggered, setting default units')
        unit = 'bohr'
        calc_mode = 'RELA'
    print(unit)
    
    # the fourth line gives the a,b,c,alpha,beta,gamma values
    # we can use these to obtain the lattice vectors
    line = fcont[3]
    elem = [s for s in line.split(' ') if s] # split line and remove whitespace
    a = float(elem[0])
    b = float(elem[1])
    c = float(elem[2])
    # angles are converted to radians
    degToRad = np.pi / 180
    alpha = float(elem[3]) * degToRad
    beta = float(elem[4]) * degToRad
    gamma = float(elem[5]) * degToRad
    
    # all values are rounded to 8 deimal places, before being premultiplied by their magnitude
    vec_a = a * np.array([1,0,0])
    vec_b = b * np.array([np.round(np.cos(gamma),decimals=8),np.round(np.sin(gamma),decimals=8),0])
    
    l = np.round(np.cos(beta) , decimals=8)
    m = np.round(( np.cos(alpha) - (np.cos(gamma) * np.cos(beta)) ) / ( np.sin(gamma) ),decimals = 8)
    n = np.round(np.sqrt( 1 - l**2 - m**2 ),decimals=8)
    vec_c = c* np.array([l,m,n])
    
    avec = np.array([vec_a,vec_b,vec_c])
    # need to perform unit change from bohr to angstrom in this.
    if unit == 'bohr':
        avec = avec / uBohr
    avec = avec * uAngs
    
    return avec,a,b,c
    
def get_z_from_file():
    '''
    Function to load in the atomic_mass.txt file and put its contents into a dictionary
    '''
    f = open('chinook/atomic_mass.txt','r')
    fcont = [line for line in f.read().split('\n') if line]
    f.close()
    Zdict = {}
    for line in fcont:
        elem = [e for e in line.split() if e]
        Zdict[elem[1]] = int(elem[0])
    return Zdict
    

    
def get_atoms_from_struct(case,avec):
    '''
    Function for extracting the numebr of atoms in the basis and their positions,
    using the case.struct file.
    This will then be used to set up the dictionary for use in the basis object
    '''
    
    # need to open the file and extract its contents
    f = open('{}.struct'.format(case),'r')
    fcont = f.read().split('\n') # split the file into lines
    f.close()
    print('get_atoms_from_struct: Extracted data from {}.struct'.format(case))
    
    line = fcont[2] # extract the third line, to determine the calculation type for the atomic positions
    if 'MODE OF CALC' in line:
        # If we have the correct line, we take the mode of calculation, rela or abs
        calc_mode = [e for e in line.split(' ') if e][2].split('=')[1].lower()
    
    Zdict = get_z_from_file()
    
    atoms = []
    Z = {}
    pos = []
    fcont = fcont[4:]  # take the 5th line onwards, as first four are setup
    for i,line in enumerate(fcont):
        if 'ATOM' in line: # if the word atom is in the line, we've found a non-equivelent atom
            j = 0
            line = line[4:] # takes off the characters 'ATOM' from the line
            while 'LOCAL ROT MATRIX' not in line: # keeps going for all the atoms of one type
                elem = [e for e in line.split(' ') if e] # splits the line by spacing and removes whitespaces
                atom_index = elem[0][1:-1] # removes the minus sign from the front and trailing colon
                
                # need to extract the X, Y and Z values: (rounded to 8 dp)
                # these values are relative to the lattice vectors
                x = np.round(float(elem[1].split('=')[1]),decimals=8)
                y = np.round(float(elem[2].split('=')[1]),decimals=8)
                z = np.round(float(elem[3].split('=')[1]),decimals=8)
    
                # NEED TO LOOK INTO WHAT TO DO IF NOT RELATIVE!!!!
                if calc_mode == 'rela':
                    atom_pos = x*avec[0] + y*avec[1] + z*avec[2]
                else:
                    print('get_atoms_from_struct: calc_mode not RELA, send help!')
                    return
                
                atoms.append(atom_index)
                pos.append(atom_pos)
                
                j = j+ 2
                line = fcont[i+j]
            # once the line with 'LOCAL ROT MATRIX' has been reached, need to go back one line and extract atom type
            line = fcont[i+j-1]
            elem = [e for e in line.split(' ') if e]
            
            atom_type = elem[0] # extracts the string for the atom type
            Z[atom_index] = Zdict[atom_type] # get_Z is placeholder for getting Z from atom string.
            
            i = i+j
            
    return atoms,Z,pos
            
               
def get_orbitals_from_atoms(Z):
    '''
    Function to get lists with all the required orbiatals for atoms taken from struct
    Z is the dictionary used to create the basis
    '''
    orbitals_Z = {}
    
    # This dictionary from orbital.py has a list of all the names for the allowed orbitals
    configOrbits = {'s':[],'p':[],'d':[],'f':[]}
    orbitalSet = orblib.projdict.keys()
    for orbital in orbitalSet:
        if orbital[0] == '0':
            configOrbits['s'].append(orbital)
        elif orbital[0] == '1':
            configOrbits['p'].append(orbital)
        elif orbital[0] == '2':
            configOrbits['d'].append(orbital)
        elif orbital[0] == '3':
            configOrbits['f'].append(orbital)
    
    
    for i,z in Z.items(): # in Z, key is index/species, value is Z
        z_orbitals = []
        
        config = eclib.get_con(eclib.filename, z)
        config = eclib.shield_split(config) # splits the configuration into nl components
        
        for state in config:
            state_n = state[0]
            
            state_l = state[1]
            
            for orbital in configOrbits[state_l]: # for each orbital associated with the l value
                z_orbitals.append('{}{}'.format(state_n,orbital))

        orbitals_Z[i] = z_orbitals
        
    return orbitals_Z

def create_basisObject(case,avec):
    '''
    Function to create a basis object from all the information in the struct file
    Will currently run a non-spin calculation case
    '''
    # get the basis vectors, get the Z
    atoms,Z,pos = get_atoms_from_struct(case,avec)
    
    z_orbitals = get_orbitals_from_atoms(Z)
    
    orbs = []
    for a in atoms:
        orbs.append(z_orbitals[a])
    
    spin = {'bool':False}
    
    basis = {
        'atoms':atoms,
        'Z':Z,
        'orbs': orbs,
        'pos': pos,
        'spin': spin}
                
    return buildlib.gen_basis(basis),spin
            
                
            
def load_qtl(case):
    '''
    Loads in the relative charge contributions from the bands from the case.qtl file
    
    After the header (which describes all the available orbitals), the file is
    designated as follows:
        
        For n atoms, there are n+1 lines per point along the k-path, per band:
            First point dictates energy
            Second value dictates to which atom (n+1 => interstitial)
            Then each value after corresponds to charge from orbital described in header.
    '''
    
    f = open('{}.qtl'.format(case),'r')
    
    #----- MIGHT NEED TO FIND APPROACH TO LOAD IN BITS OF FILE AT A TIME, memory intensive
    flines = [l for l in f.read().split('\n') if l] # reads contents of file and splits by new line
    f.close()
    
    linei = 1 # skip first line as is just a description
    line = flines[linei]
    
    E_F = 0
    while 'JATOM' not in line:
        if 'FERMI ENERGY' in line:
            E_F = float( [e for e in line.split('=') if e][2] )
        linei += 1
        line = flines[linei]
    
    # now we've reached line with JATOM in it, we need to see what orbitals are available for each atom
    n_atoms = 0
    orbitals = {}
    while 'JATOM' in line:
        n_atoms += 1
        elements = line.split(',')
        #atom_n = int([e for e in elements[0].split(' ') if e][1]) # ensures that the JATOM tag is used proeperly
        atom_n = n_atoms # This should match the JATOM tag. If qtl file exists where thats not the case, use line above.
        elements[0] = 'tot'
        
        # for each element in orbitals, stores array for string orbital names
        orbitals[atom_n] = elements
        '''# debug print statements
        print('load_qtl: atom_n = {}'.format(atom_n))
        print('load_qtl: atom_n orbitals = {}'.format(elements))
        print('load_qtl: orbitals[atom_n] = {}'.format(orbitals[atom_n]))
        '''
        # cycles through all the lines where JATOM is present
        linei += 1
        line = flines[linei]
        
    # at the end, includes the final element in the orbitals dict, for interstitial charge
    orbitals[0] = ['tot'] # uses 0 because easier later, see modulo calcs
    
    #print('load_qtl: orbitals = {}'.format(orbitals))
    
    bands = {} # will store the energies of each band
    QTL = {} # will store the partial charges of each band
    bandnum = 1
    index_in_band = 0
    # -------- SEPERATE CASES FOR FIRST BAND AND SUBSEQUENT BANDS --------
    # incase not at line where BAND 1 is written
    while 'BAND' not in line:
        linei += 1
        line = flines[linei]
    
    bands[1] = [] # empty array for energies
    QTL[1] = orbitals # the same size as orbitals, thus, for each atom, each orbital, we can create a new list
    for atom in orbitals.keys():
        print('load_qtl: first band: atom = {}'.format(atom))
        print('load_qtl: first band: QTL[{}][{}] = {}'.format(bandnum,atom,QTL[bandnum][atom]))
        for i,o in enumerate(QTL[bandnum][atom]):
            QTL[bandnum][atom][i] = [] # rather than string for orbital name, replaced with array for QTL at energy for band
   
    
    band_start_line = linei # this will be used to determine which array values are being appendeed to (atom)
    linei += 1
    line = flines[linei]
    print('load_qtl: Loading band 1')
    while 'BAND' not in line: # assuming more than one band in file, will go until band 2
        elements = [e for e in line.split(' ') if e]
        atom = (linei - band_start_line)%(n_atoms + 1) # 0 for interstitial
        line_qtl = elements[2:]
        for i,q in enumerate(line_qtl):
            QTL[1][atom][i].append(float(q)) # append line's qtl to atoms orbital-qtl list
            # QTL[band][atom][orbital][point along band] = (qtl of that atom's particluar orbital, at that point along band)
        if not atom: # if we're on the line for interstitial atoms
            bands[1].append(float(elements[0]))
            
        linei += 1
        line = flines[linei]
    
    bands[1] = np.array(bands[1]) # sets the energy for each band as a nupy array, for later calculation reasons
    
    #for each list in QTL, turns it into a numpy array:
    for atom in QTL[1].keys():
        QTL[1][atom] = np.array(QTL[1][atom])
    
    #----- NOW FIRST BAND EXTRACTED, can get all future variables of correct size first
    # should save on processing power
    
    try:
        while line: # keeps going until an empty line is reached 
            if 'BAND' in line:
                # changes the current bandnum variable to that given in the line
                bandnum = int([e for e in line.split(' ') if e][1])
                index_in_band = 0
                band_start_line = linei
                print('load_qtl: Loading band {}'.format(bandnum))
                bands[bandnum] = np.zeros(np.size(bands[1]))
                QTL[bandnum] = QTL[1]
                
            else: # for each line in the band
                elements = [e for e in line.split(' ') if e]
                atom = (linei - band_start_line)%(n_atoms + 1) # 0 for interstitial
                line_qtl = elements[2:]
                for i,q in enumerate(line_qtl):
                    QTL[bandnum][atom][i][index_in_band] = float(q) # append line's qtl to atoms orbital-qtl list
                    # QTL[band][atom][orbital][point along band] = (qtl of that atom's particluar orbital, at that point along band)
                if not atom: # if we're on the line for interstitial atoms
                    bands[bandnum][index_in_band] = float(elements[0])
                    index_in_band +=1
                
            linei += 1
            line = flines[linei]
    except:
        print('load_qtl: exception for linei={}'.format(linei))
    print('load_qtl: All bands loaded.')
    
    for b in bands.keys():
        bands[b] = bands[b] * Ry_to_eV # as bands[j] should be np.array, can do simple multiplication
    
    E_F = E_F * Ry_to_eV
    
    return QTL,bands,orbitals,E_F
            
    
    
    
    