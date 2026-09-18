import numpy as np
import subprocess, os
def convert_coordinate(coordinate, v1=None, v2=None, v3=None, v4=None, \
                       *, mode=0, nthreads=1, tmp_dir=None):
    # -----------------------------------------------------
    if (isinstance(mode, str)):
        if   mode == 'xyz_to_lon_lat_r': 
            mode= 0
        elif mode == 'lon_lat_r_to_xyz':
            mode= 1
        if   mode == 'xyz_to_lon2_lat2_r': 
            mode= 2
        elif mode == 'lon2_lat2_r_to_xyz':
            mode= 3
        elif mode == 'lon_lat_r_to_lon2_lat2_r':
            mode= 4
        elif mode == 'lon2_lat2_r_to_lon_lat_r':
            mode= 5
        else: Exception('Something is wrong with mode')
    # -----------------------------------------------------
    if 'dtype' not in dir(coordinate): coordinate=np.array(coordinate,'f8')
    if coordinate.shape[-1] != 3: raise Exception('Something is wrong with coordinate')

    r4flag = coordinate.dtype != 'f8'
    if r4flag and coordinate.dtype != 'f8': coordinate=np.array(coordinate,'f4')
    # -----------------------------------------------------
    present1= v1 is not None
    present2= v2 is not None
    present3= v3 is not None
    present4= v4 is not None

    if present1:
        if 'dtype' not in dir(v1): v1=np.array(v1, coordinate.dtype)
        if coordinate.shape != np.array(v1).shape: raise Exception('Something is wrong with v1')
    if present2:
        if 'dtype' not in dir(v2): v2=np.array(v2, coordinate.dtype)
        if coordinate.shape != np.array(v2).shape: raise Exception('Something is wrong with v2')
    if present3:
        if 'dtype' not in dir(v3): v3=np.array(v3, coordinate.dtype)
        if coordinate.shape != np.array(v3).shape: raise Exception('Something is wrong with v3')
    if present4:
        if 'dtype' not in dir(v4): v4=np.array(v4, coordinate.dtype)
        if coordinate.shape != np.array(v4).shape: raise Exception('Something is wrong with v4')
    # -----------------------------------------------------
    cdir = os.getcwd()+os.sep 
    if tmp_dir is not None: 
        if tmp_dir[-1] != os.sep : tmp_dir=tmp_dir+os.sep
    else: tmp_dir= cdir+'tmpFastQSL'+os.sep

    old_tmp_dir=os.path.exists(tmp_dir)
    if not old_tmp_dir: os.makedirs(tmp_dir, exist_ok=True)
    # -----------------------------------------------------
    with open(tmp_dir+'head.bin','wb') as file: 
        file.write(np.array([mode, nthreads, r4flag],'i4'))
        file.write(np.array(coordinate.size,'i8'))
    with open(tmp_dir+'coordinate.bin','wb') as file: file.write(coordinate)
    if present1: 
        with open(tmp_dir+'v1.bin','wb') as file:
            file.write(np.array(v1, dtype=coordinate.dtype, order='C'))
    if present2: 
        with open(tmp_dir+'v2.bin','wb') as file: 
            file.write(np.array(v2, dtype=coordinate.dtype, order='C'))
    if present3: 
        with open(tmp_dir+'v3.bin','wb') as file: 
            file.write(np.array(v3, dtype=coordinate.dtype, order='C'))
    if present4: 
        with open(tmp_dir+'v4.bin','wb') as file: 
            file.write(np.array(v4, dtype=coordinate.dtype, order='C'))
    # -----------------------------------------------------
    # please specify the path
    os.chdir(tmp_dir)
    subprocess.run(r'/path/of/convert_coordinate.x', shell=True)
    os.chdir(cdir)
    # -----------------------------------------------------
    with open(tmp_dir+'coordinate_out.bin','rb') as file: 
        coordinate_out=np.fromfile(file, dtype=coordinate.dtype).reshape(coordinate.shape)
    if present1:
        with open(tmp_dir+'v1out.bin','rb') as file: 
            v1out=np.fromfile(file, dtype=coordinate.dtype).reshape(coordinate.shape)
    if present2:
        with open(tmp_dir+'v2out.bin','rb') as file: 
            v2out=np.fromfile(file, dtype=coordinate.dtype).reshape(coordinate.shape)
    if present3:
        with open(tmp_dir+'v3out.bin','rb') as file: 
            v3out=np.fromfile(file, dtype=coordinate.dtype).reshape(coordinate.shape)
    if present4:
        with open(tmp_dir+'v4out.bin','rb') as file: 
            v4out=np.fromfile(file, dtype=coordinate.dtype).reshape(coordinate.shape)
    # -----------------------------------------------------
    if old_tmp_dir:
        os.remove(tmp_dir+'head.bin')
        os.remove(tmp_dir+'coordinate.bin')
        os.remove(tmp_dir+'coordinate_out.bin')
        if present1: 
            os.remove(tmp_dir+'v1.bin')
            os.remove(tmp_dir+'v1out.bin')
        if present2: 
            os.remove(tmp_dir+'v2.bin')
            os.remove(tmp_dir+'v2out.bin')
        if present3: 
            os.remove(tmp_dir+'v3.bin')
            os.remove(tmp_dir+'v3out.bin')
        if present4: 
            os.remove(tmp_dir+'v4.bin')
            os.remove(tmp_dir+'v4out.bin')
    else: os.rmdir(tmp_dir)
    # -----------------------------------------------------
    convert_out= coordinate_out
    if present1: convert_out=(convert_out,).__add__((v1out,))
    if present2: convert_out=   convert_out.__add__((v2out,))
    if present3: convert_out=   convert_out.__add__((v3out,))
    if present4: convert_out=   convert_out.__add__((v4out,))

    return convert_out