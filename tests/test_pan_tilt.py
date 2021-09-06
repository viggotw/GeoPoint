from numpy import tan, sqrt, pi
from geopoint.pan_tilt import xyz_to_pan_tilt

# Test decomposition
def test_rotation_matrices_decomposition_vector_1_0_0():
    assert xyz_to_pan_tilt(1,0,0) == (0, 0)

def test_rotation_matrices_decomposition_vector_0_1_0():
    assert xyz_to_pan_tilt(0,1,0) == (90, 0)

def test_rotation_matrices_decomposition_vector_0_0_1():
    assert xyz_to_pan_tilt(0,0,1) == (0, 90)

def test_rotation_matrices_decomposition_vector_1_1_0():
    assert xyz_to_pan_tilt(1,1,0) == (45, 0)

def test_rotation_matrices_decomposition_vector_1_0_1():
    assert xyz_to_pan_tilt(1,0,1) == (0, 45)

def test_rotation_matrices_decomposition_vector_0_1_1():
    assert xyz_to_pan_tilt(0,1,1) == (90, 45)

def test_rotation_matrices_decomposition_angle_45_45():
    assert xyz_to_pan_tilt(1,1,tan(pi/4)*sqrt(2)) == (45, 45)

def test_rotation_matrices_decomposition_vector_n1_0_0():
    assert xyz_to_pan_tilt(-1,0,0) == (180, 0)

def test_rotation_matrices_decomposition_vector_0_0_n1():
    assert xyz_to_pan_tilt(0,0,-1) == (0, -90)

def test_rotation_matrices_decomposition_vector_0_n1_n1():
    assert xyz_to_pan_tilt(0,-1,-1) == (-90, -45)