import numpy as np
from geopoint.geo_convert import geodetic_to_enu
from geopoint.rotational_matrices import offset_north_pitch_roll

# Decompose pan tilt angles from a cartesian coordinate system
def cartesian_to_spherical(x, y, z):
    # pan = tan^-1( y/x )
    azimuth = np.rad2deg(
        np.arctan2(y, x)
    )
    # tilt = tan^-1( sqrt(x^2 + y^2) )
    elevation = np.rad2deg(
        np.arctan2(z, np.sqrt(x**2 + y**2))
    )
    # radius = sqrt(x^2 + y^2 + z^2)
    r = np.sqrt(x**2 + y**2 + z**2)
    
    return azimuth, elevation, r

def spherical_to_cartesian(azimuth, elevation, r):
    azimuth = np.deg2rad(azimuth)
    elevation = np.deg2rad(elevation)
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z

def xyz_to_pan_tilt(x, y, z):
    pan, tilt, _ = cartesian_to_spherical(x, y, z)
    return pan, tilt

def geodetic_to_pan_tilt_with_offsets(target_lat, target_lon, target_h, origin_lat, origin_lon, origin_h, north_offset, roll_offset, pitch_offset):
    # Place the target in an ENU coordinate centered in the origin's position
    target_enu_x, target_enu_y, target_enu_z = geodetic_to_enu(target_lat, target_lon, target_h, origin_lat, origin_lon, origin_h)

    # Adjust the target's ENU coordinates based on offsets relative to north, roll and pitch (in that order).
    # This gives us the target's body coordinates.
    x_body, y_body, z_body = offset_north_pitch_roll(target_enu_x, target_enu_y, target_enu_z, north_offset, roll_offset, pitch_offset)

    # Decompose the XYZ coordinates to pan and tilt angles
    pan, tilt = xyz_to_pan_tilt(x_body, y_body, z_body)
    return pan, tilt

def geodetic_to_pan_tilt(target_lat, target_lon, target_h, origin_lat, origin_lon, origin_h, north_offset):
    return geodetic_to_pan_tilt_with_offsets(target_lat, target_lon, target_h, origin_lat, origin_lon, origin_h, north_offset, 0, 0)

def pan_tilt_to_target_geodetic(pan, tilt, dist, max_view_distance):
    """
    NOT IMPLEMENTED!
    The function should construct a dome around the PU with radius=max_view_distances.
    The target (lat,lon,h) is then calculated as the intersection between the ray given by
    the pan/tilt angle and on of the following:
    a. earth's surfce
    b. the dome's surface
    c. The dist parameter
    """
    pass
