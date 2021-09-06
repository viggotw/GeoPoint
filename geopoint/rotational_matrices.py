# The coordinate transformation from ENU to pan/tilt is based on the "kinematics of moving frames" from:
# https://ocw.mit.edu/courses/mechanical-engineering/2-017j-design-of-electromechanical-robotic-systems-fall-2009/course-text/MIT2_017JF09_ch09.pdf

import numpy as np
from functools import reduce

# Rotational matrices for rotating the coordinate system
def Rx(deg):
    sin_x = np.sin(np.deg2rad(deg))
    cos_x = np.cos(np.deg2rad(deg))
    return np.array([[1,      0,     0],
                     [0,  cos_x, sin_x],
                     [0, -sin_x, cos_x]])

def Ry(deg):
    sin_x = np.sin(np.deg2rad(deg))
    cos_x = np.cos(np.deg2rad(deg))
    return np.array([[cos_x, 0, -sin_x],
                     [    0, 1,      0],
                     [sin_x, 0,  cos_x]])

def Rz(deg):
    sin_x = np.sin(np.deg2rad(deg))
    cos_x = np.cos(np.deg2rad(deg))
    return np.array([[ cos_x, sin_x, 0],
                     [-sin_x, cos_x, 0],
                     [     0,     0, 1]])

def rotate_frame_zyx(x, y, z, z_offset, y_offset, x_offset):
    """
    Rotates the coordinate system along the Z-, Y-, and X-axis,
    and returns the new coordinates of (x,y,z) expressed in this
    new, shifted coordiante system
    """
    vector = np.array([[x],[y],[z]])
    vector_offset = reduce(np.dot, [Rx(x_offset), Ry(y_offset), Rz(z_offset+90), vector])
    x_new, y_new, z_new = (vector_offset[i,0] for i in range(3))
    return x_new, y_new, z_new

def offset_north_pitch_roll(x, y, z, north_offset, pitch_offset, roll_offset):
    return rotate_frame_zyx(x, y, z, north_offset, pitch_offset, roll_offset)