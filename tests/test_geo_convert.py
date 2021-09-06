from geopoint.geo_convert import geodetic_to_ecef, geodetic_to_enu, ecef_to_geodetic, ecef_to_enu, enu_to_ecef

# Helper functions
def are_close(x0, x1):
    d = x1 - x0
    return (d * d) < 0.1

# Constants
latLA = 34.00000048
lonLA = -117.3335693
hLA = 251.702

################################################
# Testing against results in paper             #
################################################
def test_position_of_LA_ecef():
    x0, y0, z0 = geodetic_to_ecef(latLA, lonLA, hLA)
    assert are_close(x0, -2430601.8)
    assert are_close(y0, -4702442.7)
    assert are_close(z0,  3546587.4)

def enu_LA(dx, dy, dz):
    x, y, z = geodetic_to_ecef(latLA, lonLA, hLA)
    return ecef_to_enu(x+dx, y+dy, z+dz, latLA, lonLA, hLA)

def test_enu_LA():
    x_east, y_north, z_up = enu_LA(0, 0, 0)
    assert are_close(x_east,  0.0)
    assert are_close(y_north, 0.0)
    assert are_close(z_up,    0.0)

def test_enu_LA_dx():
    x_east, y_north, z_up = enu_LA(1, 0, 0)
    assert are_close(x_east,   0.88834836)
    assert are_close(y_north,  0.25676467)
    assert are_close(z_up,    -0.38066927)

def test_enu_LA_dy():
    x_east, y_north, z_up = enu_LA(0, 1, 0)
    assert are_close(x_east,  -0.45917011)
    assert are_close(y_north,  0.49675810)
    assert are_close(z_up,    -0.73647416)

def test_enu_LA_dz():
    x_east, y_north, z_up = enu_LA(0, 0, 1)
    assert are_close(x_east,  0.00000000)
    assert are_close(y_north, 0.82903757)
    assert are_close(z_up,    0.55919291)

def test_full_conversion_loop_with_LA():
    x, y, z = geodetic_to_ecef(latLA, lonLA, hLA)
    x_east, y_north, z_up = ecef_to_enu(x, y, z, latLA, lonLA, hLA)
    x_test, y_test, z_test = enu_to_ecef(x_east, y_north, z_up, latLA, lonLA, hLA)
    lat_test, lon_test, h_test = ecef_to_geodetic(x_test, y_test, z_test)

    assert are_close(latLA, lat_test)
    assert are_close(lonLA, lon_test)
    assert are_close(hLA, h_test)
    


################################################
# geodetic_to_enu                              #
################################################
# Target equals origin
def test_geodetic_to_enu_height_0m():
    assert geodetic_to_enu(0,0,0,0,0,0) == (0,0,0)

# Varying heights for origin
def test_geodetic_to_enu_origin_height_10m():
    assert geodetic_to_enu(0,0,0,0,0,10) == (0,0,-10)

def test_geodetic_to_enu_origin_height_n10m():
    assert geodetic_to_enu(0,0,0,0,0,-10) == (0,0,10)

# Varying heights for target
def test_geodetic_to_enu_target_height_10m():
    assert geodetic_to_enu(0,0,10,0,0,0) == (0,0,10)

def test_geodetic_to_enu_target_height_n10m():
    assert geodetic_to_enu(0,0,-10,0,0,0) == (0,0,-10)

# Varying latitudes for target
def test_geodetic_to_enu_varying_latitudes_90deg():
    enu_x, enu_y, enu_z = geodetic_to_enu(90,0,0,0,0,0)
    assert enu_x == 0 and enu_y > 0 and enu_z < 0

def test_geodetic_to_enu_varying_latitudes_n90deg():
    enu_x, enu_y, enu_z = geodetic_to_enu(-90,0,0,0,0,0)
    assert enu_x == 0 and enu_y < 0 and enu_z < 0

# Varying longitudes for target
def test_geodetic_to_enu_varying_longitudes_90deg():
    enu_x, enu_y, enu_z = geodetic_to_enu(0,90,0,0,0,0)
    assert enu_x > 0 and enu_y == 0 and enu_z < 0

def test_geodetic_to_enu_varying_longitudes_n90deg():
    enu_x, enu_y, enu_z = geodetic_to_enu(0,-90,0,0,0,0)
    assert enu_x < 0 and enu_y == 0 and enu_z < 0

def test_geodetic_to_enu_varying_longitudes_180deg():
    enu_x, enu_y, enu_z = geodetic_to_enu(180,0,0,0,0,0)
    assert enu_x == 0 and round(enu_y, 6) == 0 and enu_z < 0

# All four quadrants
def test_geodetic_to_enu_quadrant_1():
    enu_x, enu_y, enu_z = geodetic_to_enu(90,90,0,0,0,0)
    assert enu_x > 0 and enu_y > 0 and enu_z < 0

def test_geodetic_to_enu_quadrant_2():
    enu_x, enu_y, enu_z = geodetic_to_enu(-90,90,0,0,0,0)
    assert enu_x > 0 and enu_y < 0 and enu_z < 0

def test_geodetic_to_enu_quadrant_3():
    enu_x, enu_y, enu_z = geodetic_to_enu(-90,-90,0,0,0,0)
    assert enu_x < 0 and enu_y < 0 and enu_z < 0

def test_geodetic_to_enu_quadrant_4():
    enu_x, enu_y, enu_z = geodetic_to_enu(90,-90,0,0,0,0)
    assert enu_x < 0 and enu_y > 0 and enu_z < 0
