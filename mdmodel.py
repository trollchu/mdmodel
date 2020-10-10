import numpy as np
import pandas as pd
import os
import re
import math
import itertools
"""
index 'id' columns ['type','charge','x','y','z','vx','vy','vz']
"""


class mdmodel:
# --------------------------------------------------------------------

    def __init__(self,name='cqz'):
        self.name = name
        self.natom = 0
        self.nbond = 0
        self.nangle = 0
        self.ndihedral = 0
        self.nimproper = 0
        self.natomtype = 0
        self.nbondtype = 0
        self.nangletype = 0
        self.ndihedraltype = 0
        self.nimpropertype = 0
        # self.box = [[-80,80],[-80,80],[-80,80]]
        self.box=[]
        self.mass = []
        self.mass_index = False

    
    #----------------------------------------------
    # read molecular topos
    #
    #
    def read_lmp_data(self,filename,atom_style='charge',\
                    read_velocities = True):
        
        with open(filename) as f:
            tex = f.readlines()
            
        # read number of atoms, atom types, box, masses
        index = 1
        v_index = 0
        for i in tex:
            index+=1
            i=i.strip()
            if re.match(r'\d+\s+atoms',i):
                self.natom = int(re.match(r'\d+',i).group())
            if re.match(r'\d+\s+bonds',i):
                self.nbond = int(re.match(r'\d+',i).group())
            if re.match(r'\d+\s+angles',i):
                self.nangle = int(re.match(r'\d+',i).group())
            if re.match(r'\d+\s+dihedrals',i):
                self.ndihedral = int(re.match(r'\d+',i).group())
            if re.match(r'\d+\s+impropers',i):
                self.nimproper = int(re.match(r'\d+',i).group())

            if re.match(r'\d+\s+atom\stypes',i):
                self.natomtype = int(re.match(r'\d+',i).group())
            if re.match(r'\d+\s+bond\stypes',i):
                self.nbondtype = int(re.match(r'\d+',i).group())
            if re.match(r'\d+\s+angle\stypes',i):
                self.nangletype = int(re.match(r'\d+',i).group())
            if re.match(r'\d+\s+dihedral\stypes',i):
                self.ndihedraltype = int(re.match(r'\d+',i).group())
            if re.match(r'\d+\s+improper\stypes',i):
                self.nimpropertype = int(re.match(r'\d+',i).group())

            if re.match(r'\s*[\+\-\w\.]+\s+[\+\-\w\.]+\s+[xyzhilo]+',i):
                self.box.append(i.split()[:2])

            if re.match(r'\s*Masses\s*',i):
                self.mass_index=index

            if re.match(r'\s*Atoms\s*',i):
                x_index=index


            if re.match(r'\s*Velocities\s*',i):
                v_index=index

            if re.match(r'\s*Bonds\s*',i):
                bond_index=index            
            if re.match(r'\s*Angles\s*',i):
                angle_index=index   
            if re.match(r'\s*Dihedrals\s*',i):
                dihedral_index=index
            if re.match(r'\s*Impropers\s*',i):
                improper_index=index         
                      
        # save atoms x,y,z,vx,vy,vz   
        # print(self.ntype)    
        # print(tex) 
        if atom_style == 'charge':
            col_name=['type','charge','x','y','z']
            ncol = 6
        elif atom_style == 'atomic':
            col_name=['type','x','y','z']
            ncol = 5
        elif atom_style == 'full':
            col_name=[ 'molid' ,'type','charge','x','y','z']
            ncol = 7
        else:
            print('unsupported atom_style, (only charge and atomic)')

        
        # print('Reading lammps data:{}'.format(filename))
        self.data = pd.read_csv(filename,skiprows=x_index,header=None,index_col=0,sep='\s+',nrows=self.natom,usecols=range(ncol))
        self.data.columns=col_name

        self.data.type = self.data.type.astype(int)
        self.data['id'] = self.data.index
        # sort data by atom id
        self.data = self.data.sort_index()

        if v_index>0 and read_velocities:
            self.velocities = pd.read_csv(filename,skiprows=v_index,header=None,index_col=0,nrows=self.natom,sep='\s+')
            self.velocities.columns = ['vx','vy','vz']
            self.data[['vx','vy','vz']]=self.velocities[['vx','vy','vz']]
            self.velocities['id'] = self.velocities.index


        try:
            self.bonds = pd.read_csv(filename,skiprows = bond_index, header= None,index_col=0,nrows=self.nbond,sep='\s+')
        except:
            pass

        try:
            self.angles = pd.read_csv(filename,skiprows = angle_index, header= None,index_col=0,nrows=self.nangle,sep='\s+')
        except:
            pass

        try:
            self.dihedrals = pd.read_csv(filename,skiprows = dihedral_index, header= None,index_col=0,nrows=self.ndihedral,sep='\s+')
        except:
            pass

        try: 
            self.impropers = pd.read_csv(filename,skiprows = improper_index, header= None,index_col=0,nrows=self.nimproper,sep='\s+')
        except:
            pass
            
        if self.mass_index:
            self.mass = pd.read_csv(filename,skiprows = self.mass_index, header= None,nrows=self.natomtype,sep='\s+')
        

        
    def read_gmx_gro(self,filename):
        with open(filename) as fo:
            fo.readline()
            self.natom = int(fo.readline())
        
        self.data = pd.read_csv(filename,skiprows=2,nrows=self.natom,header=None,sep='\s+')
        self.data.columns=['resi_name','atom_name','id','x','y','z','vx','vy','vz']
        self.data[['x','y','z']]=self.data[['x','y','z']]*10
        self.data[['vx','vy','vz']]=self.data[['vx','vy','vz']]*10/1000 #nm/ps to A/fs
        
        def extract_atom_type(x):
            y = re.findall(r'[A-z]+',x)[0]
            return y
        
        self.data['type'] = self.data.atom_name.apply(extract_atom_type)
        def process_row(row):
            match = re.match('C\d+',str(row))
            if match:
                return 1
            else:
                return 2
        ## 补全信息
        self.natomtype = 2
        self.mass = [['1', '12.011'], ['2', '1.008']]
        self.data['charge'] =0
        self.data['placehold1'] = 0
        self.data['placehold2'] =0
        self.data['placehold3'] = 0
        self.data['type']=self.data['atom_name'].apply(process_row)
    
    def read_pdb(self,filename,natom):
        # with open(filename) as fo:
        #     ff = fo.readlines()
        self.natom = natom
        
        self.data = pd.read_csv(filename,skiprows=3,nrows=self.natom,header=None,sep=r'\s+',usecols=range(8))
        self.data.columns=['resi_name','id','atom_name','mol_name','mol_id','x','y','z']

        def process_row(row):
            match = re.match('C\d+',str(row))
            if match:
                return 1
            else:
                return 2
        ## 补全信息
        self.natomtype = 2
        self.mass = [['1', '12.011'], ['2', '1.008']]
        self.data['charge'] =0

        self.data['type']=self.data['atom_name'].apply(process_row)

    def read_gjf(self, filename):
        with open(filename) as fo:
            text = fo.readlines()
        pos = []
        for i in text:
            # if re.match(r'\s* [\w\(\=\)]*\s[\w\-\.]+',i):
            if re.match(r'[\w]+ \s+ [\w\-\.]+ \s+ [\w\-\.]+ \s+ [\w\-\.]+',i.strip()):
                pos.append(i.split())

        self.data = pd.DataFrame(pos,dtype = float)

        self.data.columns=['atom_name','x','y','z']
        self.data['charge'] = 0
        self.data['id']=self.data.index+1
        def f(df):
    
            if    df.atom_name =='C':
                return 1
            else:
                return 2
        self.data['type'] = self.data.apply(f,axis=1)


    def read_lmp_dump(self, filename, timestep=0,column_name= None):
        with open(filename) as f:
            data  = f.readlines()
        self.box=[[-80,80],[-80,80],[-80,80]]
        self.ntimestep=data.count('ITEM: TIMESTEP\n')
        # get number of atoms
        with open(filename) as f:
            while True:
                ff = f.readline()
                if ff == 'ITEM: TIMESTEP\n':
                    self.timestep  = int(f.readline().strip('\n'))                
                if ff == 'ITEM: NUMBER OF ATOMS\n':
                    self.natom  = int(f.readline().strip('\n'))
                    break

        lines_of_header = 9 # default for lammps dump file

        index = timestep*(lines_of_header+self.natom)+lines_of_header


        tex = data[index:index+self.natom]
        if not column_name:
            column_name = data[lines_of_header-1].strip('ITEM: ATOMS ').strip().split()

        self.data = pd.read_csv(filename,skiprows=index,nrows=self.natom,header=None,sep='\s+')
        self.data.index = self.data[0]
        self.data.columns=column_name
        # sort data by atom id
        self.data = self.data.sort_index()



    def read_xyz(self, filename):
        with open(filename) as fo:
            ff = fo.readlines()
        natom = int(ff[0].strip())
        self.natom = natom
        
        self.data = pd.read_csv(filename,skiprows=2,nrows=self.natom,header=None,sep=r'\s+')
        self.data.columns=['type','x','y','z']
        self.data['id'] = np.arange(natom)+1



#### 
##### write

    def write_xyz(self, filename):
        with open(filename,'w+') as f:
            f.writelines('{}\ntype x y z\n'.format(self.natom))
        
        self.data.loc[:,['atom_name','x','y','z']].to_csv(filename,sep=' ',header=False,index=False,mode='a',float_format='%f')

    def write_xyzr(self, filename,radius):

        for i,r in enumerate(radius):
            self.data.loc[self.data.type==i+1,'radius']=r
        self.data.loc[:,['atom_name','x','y','z','radius']].to_csv(filename,sep=' ',header=False,index=False,mode='w',float_format='%f')


    def write_gjf(self,filename,frame):
        # filename = 'test'
        memory = 8
        core = 6
        theory = 'm062x/6-311G(d,p)'
        # ccsd(t)/cc-pvtz  B3LYP/6-311+G*

        fout=open('%s.gjf' %(filename),'w')
        fout.write("%")
        fout.write("rwf=%s.rwf\n" %(filename))
        fout.write("%")
        fout.write("NoSave\n")
        fout.write("%")
        fout.write("chk=%s.chk\n" %(filename))
        fout.write("%")
        fout.write("mem=%sGB\n" %(memory))
        fout.write("%")
        fout.write("nprocshared=%s\n" %(core))
        
        if frame:
            fout.write("# B3LYP/6-311+G*  counterpoise=2\n")
        else:
            fout.write("# B3LYP/6-311+G*\n")
        fout.write("\n")
        fout.write("benzene+ar")
        fout.write("\n")
        fout.write("\n")
        fout.write("0 1 0 1 0 1\n")
        # frame = '1'

        if frame:
            print('hello  you')
            pos = self.data[['atom_name','frame','x','y','z']].values
            for temp in pos:
                fout.write(" {}(fragment={}) {} {}  {} \n".format(temp[0],temp[1],temp[2],temp[3],temp[4]))

        else:
            pos = self.data[['atom_name','x','y','z']].values
            print('hello')
            for temp in pos:
                fout.write(" {} {} {}  {} \n".format(temp[0],temp[1],temp[2],temp[3]))



        fout.write("\n")
        fout.write("\n")
        fout.write("\n")
        fout.write("\n")
        fout.close()
    
    def write_lmp_dump(self,filename):
        #lmp=pd.DataFrame()
        self.data['id']=self.data['id'].astype('int')
        self.data['type']=self.data['type'].astype('int')
        #lmp[['atom_index','atom_type','charge','x','y','z']]=self.data[['atom_index','type_num','charge','x','y','z']]
        with open(filename,'w+') as f:
            f.writelines('ITEM: TIMESTEP\n')
            f.writelines('{}\n'.format(self.timestep))
            f.writelines('ITEM: NUMBER OF ATOMS\n')
            f.writelines('{}\n'.format(self.natom))
            f.writelines('ITEM: BOX BOUNDS pp pp pp\n')
            for bhl in self.box:
                f.writelines('{} {}\n'.format(bhl[0],bhl[1]))
            f.writelines('ITEM: ATOMS {}\n'.format(' '.join(self.data.columns.tolist())))
        self.data.to_csv(filename,sep=' ',header=False,index=False,mode='a',float_format='%f')
   
    def write_lmp_data(self,filename,atom_style='charge',velocity=True):

        #  1. write header
        with open(filename,'w+') as f:
            f.writelines('LAMMPS data file via md_model.write_lmp_data\n\n')
            f.writelines('{} atoms\n'.format(self.natom))

            if self.nbond:
                f.writelines('{} bonds\n'.format(self.nbond))
            if self.nangle:
                f.writelines('{} angles\n'.format(self.nangle))
            if self.ndihedral:
                f.writelines('{} dihedrals\n'.format(self.ndihedral))
            if self.nimproper:
                f.writelines('{} impropers\n'.format(self.nimproper))


            f.writelines('{} atom types\n'.format(self.natomtype))

            if self.nbondtype:
                f.writelines('{} bond types\n'.format(self.nbondtype))
            if self.nangletype:
                f.writelines('{} angle types\n'.format(self.nangletype))
            if self.ndihedraltype:
                f.writelines('{} dihedral types\n'.format(self.ndihedraltype))
            if self.nimpropertype:
                f.writelines('{} improper types\n\n'.format(self.nimpropertype))

            try:            
                f.writelines('{} {} xlo xhi\n'.format(self.box[0][0],self.box[0][1]))
                f.writelines('{} {} ylo yhi\n'.format(self.box[1][0],self.box[1][1]))
                f.writelines('{} {} zlo zhi\n\n'.format(self.box[2][0],self.box[2][1]))
            except:
                self.box=[[-80,80],[-80,80],[-80,80]]
                f.writelines('{} {} xlo xhi\n'.format(self.box[0][0],self.box[0][1]))
                f.writelines('{} {} ylo yhi\n'.format(self.box[1][0],self.box[1][1]))
                f.writelines('{} {} zlo zhi\n\n'.format(self.box[2][0],self.box[2][1]))

        # 2. write 'Masses' part  
        try:
            z = self.mass   
            with open(filename,'a') as f:
                f.writelines('\nMasses\n\n')
            pd.DataFrame(self.mass).to_csv(filename,sep=' ',header=False,index=False,mode='a')
        except:
            print('No mass data')


            

        # 3.  write 'Atoms' part

        if atom_style == 'charge':
            # columns=['id','type','charge','x','y','z','placehold1','placehold2','placehold3']
            
            columns=['id','type','charge','x','y','z']
        elif atom_style == 'atomic':
            columns=['id','type','x','y','z']
        elif atom_style == 'full':
            columns=[ 'id','molid' ,'type','charge','x','y','z']
        else:
            print('unsupported atom_style, (only charge and atomic)')
            raise


        atoms=self.data[columns]
        with open(filename,'a') as f:
            f.writelines('\nAtoms #{}\n\n'.format(atom_style))
        atoms.to_csv(filename,sep=' ',header=False,index=False,mode='a',float_format='%f')
        # 4. write 'Velocity' part
        
        if velocity:
            z = self.velocities   # 为了验证velocity是否存在
            with open(filename,'a') as f:
                f.writelines('\nVelocities\n\n')
            v=self.velocities[['id','vx','vy','vz']]

            v.to_csv(filename,sep=' ',header=False,index=False,mode='a')

           # print('No velocities data')

        # 5. write bond angle dihedral improper part
        
        try:
            z = self.bonds   
            with open(filename,'a') as f:
                f.writelines('\nBonds\n\n')
            pd.DataFrame(self.bonds).to_csv(filename,sep=' ',header=False,index=True,mode='a')
        except:
            pass
            #print('No bonds data')
        try:
            z = self.angles   
            with open(filename,'a') as f:
                f.writelines('\nAngles\n\n')
            pd.DataFrame(self.angles).to_csv(filename,sep=' ',header=False,index=True,mode='a')
        except:
            pass
            #print('No angles data')   
        try:
            z = self.dihedrals   
            with open(filename,'a') as f:
                f.writelines('\nDihedrals\n\n')
            pd.DataFrame(self.dihedrals).to_csv(filename,sep=' ',header=False,index=True,mode='a')
        except:
            pass
            #print('No dihedrals data')
        try:
            z = self.impropers   
            with open(filename,'a') as f:
                f.writelines('\nImpropers\n\n')
            pd.DataFrame(self.impropers).to_csv(filename,sep=' ',header=False,index=True,mode='a')
        except:
            pass



        ### operation functions
    def cal_r(self):
        self.data['r'] = self.data.loc[:,['x','y','z']].apply(np.linalg.norm,axis=1)

    def get_atom_name(self,atom_name):
        def _f(x):
            return atom_name[x]

        self.data['atom_name'] = self.data.type.apply(_f)


        # fast read dump
        #2019-06-04

    def read_dump_np(self, filename,n_frames=None,atom_slice=None):
        self._frame_index = 0
        # track which line we're on. this is not essential, but its useful
        # when reporting errors to the user to say what line it occured on.
        self._line_counter = 0
        self._fh = open(filename, 'r')
        self._filename = filename
        if n_frames is None:
            frame_counter = itertools.count()
        else:
            frame_counter = range(n_frames)

        all_other = []
        for _ in frame_counter:
            try:
                frame_other = self._read()
                if atom_slice is not None:
                    frame_other = frame_other[atom_slice, :]
            except _EOF:
                break

            all_other.append(frame_other)
        all_other = np.array(all_other)
        
        self.xyz = all_other
        return all_other
    def _read(self,atom_slice=None):
        """Read a single frame. """

        # --- begin header ---
        first = self._fh.readline()  # ITEM: TIMESTEP
        if first == '':
            raise _EOF()
        self._fh.readline()  # timestep
        self._fh.readline()  # ITEM: NUMBER OF ATOMS
        a= self._fh.readline()
        # print("a{}a".format(a))
        self._n_atoms = int(a)  # num atoms

        box_header = self._fh.readline().split()  # ITEM: BOX BOUNDS
        self._line_counter += 5
        if len(box_header) == 9:
            lengths, angles = self.parse_box('triclinic')
        elif len(box_header) == 6:
            lengths, angles = self.parse_box('orthogonal')
        else:
            raise IOError('lammpstrj parse error on line {0:d} of "{1:s}". '
                          'This file does not appear to be a valid '
                          'lammpstrj file.'.format(
                    self._line_counter,  self._filename))

        column_headers = self._fh.readline().split()[2:]  # ITEM: ATOMS ...
        if self._frame_index == 0:
            # Detect which columns the atom index, type and coordinates are.
            columns = {header: idx for idx, header in enumerate(column_headers)}

            try:
                self._atom_index_column = columns['id']
        #         # self._atom_type_column = columns['type']
        #         # self._xyz_columns = [columns[keywords[0]], columns[keywords[1]], columns[keywords[2]]]
        #         self._other_columns = range(5,len(column_headers))   # cqz
            except KeyError:
                raise IOError("Invalid .lammpstrj file. Must contain 'id', "
                              "'type', 'x*', 'y*' and 'z*' entries.")
        self._line_counter += 4
        # --- end header ---

        # xyz = np.empty(shape=(self._n_atoms, 3))
        # types = np.empty(shape=self._n_atoms, dtype='int')
        other = np.empty(shape=(self._n_atoms, len(column_headers)))  ## cqz
        # --- begin body ---
        
        for _ in range(self._n_atoms):
            line = self._fh.readline()
            if line == '':
                raise _EOF()
            split_line = line.split()
            try:
                atom_index = int(split_line[self._atom_index_column])
                # types[atom_index - 1] = int(split_line[self._atom_type_column])
                # xyz[atom_index - 1] = [float(split_line[column]) for column in self._xyz_columns]
                # # print(s)
                other[atom_index - 1] = split_line   # cqz
            except Exception:
                raise IOError('lammpstrj parse error on line {0:d} of "{1:s}". '
                              'This file does not appear to be a valid '
                              'lammpstrj file.'.format(
                        self._line_counter,  self._filename))
            self._line_counter += 1
        # --- end body ---
        # print(self._frame_index)
        self._frame_index += 1
        return other

    def parse_box(self, style):
        """Extract lengths and angles from a frame.

        Parameters
        ----------
        style : str
            Type of box, 'triclinic' or 'orthogonal'.

        Returns
        -------
            lengths : ndarray
            angles : ndarray

        Notes
        -----
        For more info on how LAMMPS defines boxes:
        http://lammps.sandia.gov/doc/Section_howto.html#howto_12
        """
        box = np.empty(shape=(3, 2))
        if style == 'triclinic':
            factors = np.empty(3)
            for i in range(3):
                line = self._fh.readline().split()
                box[i] = line[:2]
                factors[i] = line[2]
            xy, xz, yz = factors

            xlo = box[0, 0] - np.min([0.0, xy, xz, xy+xz])
            xhi = box[0, 1] - np.max([0.0, xy, xz, xy+xz])
            ylo = box[1, 0] - np.min([0.0, yz])
            yhi = box[1, 1] - np.max([0.0, yz])
            zlo = box[2, 0]
            zhi = box[2, 1]

            lx = xhi - xlo
            ly = yhi - ylo
            lz = zhi - zlo

            a = lx
            b = np.sqrt(ly**2 + xy**2)
            c = np.sqrt(lz**2 + xz**2 + yz**2)
            alpha = np.arccos((xy*xz + ly*yz) / (b*c))
            beta = np.arccos(xz / c)
            gamma = np.arccos(xy / b)

            lengths = np.array([a, b, c])
            angles = np.degrees(np.array([alpha, beta, gamma]))
        elif style == 'orthogonal':
            box[0] = self._fh.readline().split()  # x-dim of box
            box[1] = self._fh.readline().split()  # y-dim of box
            box[2] = self._fh.readline().split()  # z-dim of box
            lengths = np.diff(box, axis=1).reshape(1, 3)[0]  # box lengths
            angles = np.empty(3)
            angles.fill(90.0)
        return lengths, angles



    def write(self, xyz, cell_lengths, cell_angles=None, types=None, unit_set='real'):
        """Write one or more frames of data to a lammpstrj file.

        Parameters
        ----------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms to write. By convention,
            the lengths should be in units of angstroms.
        cell_lengths : np.ndarray, dtype=np.double, shape=(n_frames, 3)
            The lengths (a,b,c) of the unit cell for each frame. By convention,
            the lengths should be in units of angstroms.
        cell_angles : np.ndarray, dtype=np.double, shape=(n_frames, 3)
            The angles (\alpha, \beta, \gamma) defining the unit cell for
            each frame. (Units of degrees).
        types : np.ndarray, shape(3, ), dtype=int
            The numeric type of each particle.
        unit_set : str, optional
            The LAMMPS unit set that the simulation was performed in. See
            http://lammps.sandia.gov/doc/units.html for options. Currently supported
            unit sets: 'real'.
        """
        if not self._mode == 'w':
            raise ValueError('write() is only available when file is opened '
                             'in mode="w"')

        xyz = ensure_type(xyz, np.float32, 3, 'xyz', can_be_none=False,
                shape=(None, None, 3), warn_on_cast=False,
                add_newaxis_on_deficient_ndim=True)
        cell_lengths = ensure_type(cell_lengths, np.float32, 2, 'cell_lengths',
                can_be_none=False, shape=(len(xyz), 3), warn_on_cast=False,
                add_newaxis_on_deficient_ndim=True)
        if cell_angles is None:
            cell_angles = np.empty_like(cell_lengths)
            cell_angles.fill(90)
        cell_angles = ensure_type(cell_angles, np.float32, 2, 'cell_angles',
                can_be_none=False, shape=(len(xyz), 3), warn_on_cast=False,
                add_newaxis_on_deficient_ndim=True)
        if not types:
            # Make all particles the same type.
            types = np.ones(shape=(xyz.shape[1]))
        types = ensure_type(types, np.int, 1, 'types', can_be_none=True,
                shape=(xyz.shape[1], ), warn_on_cast=False,
                add_newaxis_on_deficient_ndim=False)

        # TODO: Support other unit sets.
        if unit_set == 'real':
            self.distance_unit == 'angstroms'
        else:
            raise ValueError('Unsupported unit set specified: {0}.'.format(unit_set))

        for i in range(xyz.shape[0]):
            # --- begin header ---
            self._fh.write('ITEM: TIMESTEP\n')
            self._fh.write('{0}\n'.format(i))  # TODO: Write actual time if known.
            self._fh.write('ITEM: NUMBER OF ATOMS\n')
            self._fh.write('{0}\n'.format(xyz.shape[1]))
            self.write_box(cell_lengths[i], cell_angles[i], xyz[i].min(axis=0))
            # --- end header ---

            # --- begin body ---
            self._fh.write('ITEM: ATOMS id type xu yu zu\n')
            for j, coord in enumerate(xyz[i]):
                self._fh.write('{0:d} {1:d} {2:8.3f} {3:8.3f} {4:8.3f}\n'.format(
                    j+1, types[j], coord[0], coord[1], coord[2]))
            # --- end body ---

class _EOF(IOError):
    pass