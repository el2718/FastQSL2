"""
coordinate: the array of coordinate to be tranformed, 
its dimesions can be (3) or (n1, 3) or (n2, n1, 3) or (n3, n2, n1, 3) 

v1, v2, v3: the array of vectors to be tranformed, 
their dimesions should be the same as coordinate

mode:
0 or 'xyz_to_lon_lat_r'
1 or 'lon_lat_r_to_xyz'
2 or 'lon_lat_r_to_lon2_lat2_r' ; not finished
3 or 'lon2_lat2_r_to_lon_lat_r' ; not finished

usage:
coordinate_out = trans_coordinate(coordinate, v1, mode='lon_lat_r_to_xyz')
coordinate_out, v1out, v2out = trans_coordinate(coordinate, v1, v2, mode='lon_lat_r_to_xyz')
"""

import numpy as np
def convert_coordinate(coordinate, v1=None, v2=None, v3=None, *, mode=0):
    # -----------------------------------------------------
    if (isinstance(mode, str)):
        if   mode == 'xyz_to_lon_lat_r': 
            mode= 0
        elif mode == 'lon_lat_r_to_xyz':
            mode= 1
        # elif mode == 'lon_lat_r_to_lon2_lat2_r':
        #     mode= 2
        # elif mode == 'lon2_lat2_r_to_lon_lat_r':
        #     mode= 3
        else: Exception('Something is wrong with mode')

    if mode == 0: two_pi=np.pi*2.
    # -----------------------------------------------------
    coor_shape=np.array(coordinate).shape
    if coor_shape[-1] != 3: raise Exception('Something is wrong with coordinate')

    coor_dim=np.array(coordinate).ndim
    ngrid=coordinate.size//3
    coor1d=np.array(coordinate).reshape((ngrid,3))
    coor1d_out=coor1d.copy()
    coor_out=coor1d[0,:].copy()
    # -----------------------------------------------------
    present1= v1 is not None
    present2= v2 is not None
    present3= v3 is not None

    if present1: 
        v_shape=np.array(v1).shape
        if coor_dim   != np.array(v1).ndim : raise Exception('Something is wrong with v1')
        if coor_shape != np.array(v1).shape: raise Exception('Something is wrong with v1')
        v1_1d=np.array(v1).reshape((ngrid,3))
        v1_out1d=v1_1d.copy()
    if present2: 
        v_shape=np.array(v2).shape
        if coor_dim   != np.array(v2).ndim : raise Exception('Something is wrong with v2')
        if coor_shape != np.array(v2).shape: raise Exception('Something is wrong with v2')
        v2_1d=np.array(v2).reshape((ngrid,3))
        v2_out1d=v2_1d.copy()
    if present3: 
        v_shape=np.array(v3).shape
        if coor_dim   != np.array(v3).ndim : raise Exception('Something is wrong with v3')
        if coor_shape != np.array(v3).shape: raise Exception('Something is wrong with v3')
        v3_1d=np.array(v3).reshape((ngrid,3))
        v3_out1d=v3_1d.copy()
    # -----------------------------------------------------
    for i in range(ngrid):
        coor_in=coor1d[i,:]
        if mode == 0: # 'xyz_to_lon_lat_r'
            coor_out[2]=np.linalg.norm(coor_in)
            coor_out[1]=np.arcsin(coor_in[2]/coor_out[2])

            if coor_in[0] == 0. and coor_in[1] == 0.: 
                coor_out[0]=0.
            else: 
                cos_tmp=coor_in[0]/np.linalg.norm(coor_in[0:2])
                if      cos_tmp >  1.: coor_out[0]=0.
                elif    cos_tmp < -1.: coor_out[0]=np.pi
                elif coor_in[1] >= 0.: coor_out[0]=np.acos(cos_tmp)
                else :                 coor_out[0]=two_pi-np.acos(cos_tmp)
            if present1:
                sin01=np.sin(coor_out[0:2])
                cos01=np.cos(coor_out[0:2])
        elif mode == 1: # 'lon_lat_r_to_xyz'
            sin01=np.sin(coor_in[0:2])
            cos01=np.cos(coor_in[0:2])
            coor_out=coor_in[2]*np.array([cos01[1]*cos01[0], cos01[1]*sin01[0], sin01[1]])
        coor1d_out[i,:]=coor_out
        # -----------------------------------------------------
        if present1:
            # stack e_lon, e_lat, e_r for 'xyz_to_lon_lat_r'
            matrix=np.array([[         -sin01[0],           cos01[0],       0.],\
                             [-sin01[1]*cos01[0], -sin01[1]*sin01[0], cos01[1]],\
                             [ cos01[1]*cos01[0],  cos01[1]*sin01[0], sin01[1]]])
            if mode == 1: matrix=matrix.transpose()
            v1_out1d[i,:]= np.dot(matrix, v1_1d[i,:])

        if present2: v2_out1d[i,:]= np.dot(matrix, v2_1d[i,:])
        if present3: v3_out1d[i,:]= np.dot(matrix, v3_1d[i,:])
    

    trans_return= coor1d_out.reshape(coor_shape)
    # print(coor1d_out.shape, coor_shape, trans_return.v_shape)
    if present1: trans_return=(trans_return,).__add__((v1_out1d.reshape(coor_shape),))
    if present2: trans_return=   trans_return.__add__((v2_out1d.reshape(coor_shape),))
    if present3: trans_return=   trans_return.__add__((v3_out1d.reshape(coor_shape),))

    # return coor1d_out.reshape(coor_shape)
    return trans_return
        



    
