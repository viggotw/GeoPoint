from numpy import array, array_equal
from geopoint.rotational_matrices import Rx, Ry, Rz

# Test rotation matrices
def test_rotation_matrices_1():
    v = array([[0],[0],[1]])
    assert array_equal( Rz(90).dot(v).round(9),array([[0],[0],[1]]) )

def test_rotation_matrices_2():
    v = array([[0],[0],[1]])
    assert array_equal( Rx(90).dot(Rz(90).dot(v)).round(9), array([[0],[1],[0]]) )

def test_rotation_matrices_3():
    v = array([[0],[0],[1]])
    assert array_equal( Rz(90).dot(Rx(90).dot(Rz(90).dot(v))).round(9), array([[1],[0],[0]]) )

def test_rotation_matrices_4():
    v = array([[1],[0],[0]])
    assert array_equal( Rz(-90).dot(v).round(9), array([[0],[1],[0]]) )

def test_rotation_matrices_5():
    v = array([[1],[0],[0]])
    assert array_equal( Rx(-90).dot(Rz(-90).dot(v)).round(9), array([[0],[0],[1]]) )

def test_rotation_matrices_6():
    v = array([[1],[0],[0]])
    assert array_equal( Ry(-90).dot(Rx(-90).dot(Rz(-90).dot(v))).round(9), array([[1],[0],[0]]) )