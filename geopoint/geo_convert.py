# Some helpers for converting GPS readings from the WGS84 geodetic system to a local North-East-Up cartesian axis.
#
# The scripts that provide coordinate transformations between geodetic, ecef and enu in python are based on:
# https://gist.github.com/govert/1b373696c9a27ff4c72a
# #
# The implementation here is according to the paper:
# "Conversion of Geodetic coordinates to the Local Tangent Plane" Version 2.01.
# "The basic reference for this paper is J.Farrell & M.Barth 'The Global Positioning System & Inertial Navigation'"
# Also helpful is Wikipedia: http:#en.wikipedia.org/wiki/Geodetic_datum
# Also helpful are the guidance notes here: http:#www.epsg.org/Guidancenotes.aspx

from typing import Tuple
from numpy import sqrt, sin, cos, tan, arctan2, deg2rad, rad2deg, sign

'''
Alternative notation for the nested .dot()-notation
v_new = reduce(dot, [Ry(-90), Rx(-90), Rz(-90), v])
'''


# CONSTANTS
## WGS-84 geodetic constants
a = 6378137.0  # WGS-84 Earth semimajor axis (m)

b = 6356752.314245  # Derived Earth semiminor axis (m)
f = (a - b) / a  # Ellipsoid Flatness
f_inv = 1.0 / f  # Inverse flattening

a_sq = a * a
b_sq = b * b
e_sq = f * (2 - f)  # Square of Eccentricity

def geodetic_to_ecef(lat:float, lon:float, height:float) -> Tuple[float, float, float]:
    """
    Converts WGS-84 Geodetic point (lat, lon, h) to the Earth-Centered Earth-Fixed (ECEF) coordinates (x, y, z).

    :param lat: latitude in degrees
    :param lon: longitude in degrees
    :param height: meters
    :returns: (x, y, z) in [meters, meters, meters]
    """
    # (lat, lon) in WSG-84 degrees
    lambd = deg2rad(lat)
    phi = deg2rad(lon)
    s = sin(lambd)
    N = a / sqrt(1 - e_sq * s * s)

    sin_lambda = sin(lambd)
    cos_lambda = cos(lambd)
    cos_phi = cos(phi)
    sin_phi = sin(phi)

    x = (height + N) * cos_lambda * cos_phi
    y = (height + N) * cos_lambda * sin_phi
    z = (height + (1 - e_sq) * N) * sin_lambda

    return x, y, z


def ecef_to_geodetic(x:float, y:float, z:float) -> Tuple[float, float, float]:
    """
    Converts Earth-Centered Earth-Fixed (ECEF) coordinates (x, y, z) to the WGS-84 Geodetic point (lat, lon, height)
    
    :param x: meters
    :param y: meters
    :param z: meters
    :returns: (lat, lon, height) in [degrees, degrees, meters]
    """
    eps = e_sq / (1.0 - e_sq)
    p = sqrt(x * x + y * y)
    q = arctan2((z * a), (p * b))
    sin_q = sin(q)
    cos_q = cos(q)
    sin_q_3 = sin_q * sin_q * sin_q
    cos_q_3 = cos_q * cos_q * cos_q
    phi = arctan2((z + eps * b * sin_q_3), (p - e_sq * a * cos_q_3))
    lambd = arctan2(y, x)
    v = a / sqrt(1.0 - e_sq * sin(phi) * sin(phi))
    height = (p / cos(phi)) - v

    lat = rad2deg(phi)
    lon = rad2deg(lambd)

    return lat, lon, height


def ecef_to_enu(x:float, y:float, z:float, lat0:float, lon0:float, h0:float) -> Tuple[float, float, float]:
    """
    Converts the Earth-Centered Earth-Fixed (ECEF) coordinates (x, y, z) to 
    East-North-Up coordinates in a Local Tangent Plane that is centered at the 
    (WGS-84) Geodetic point (lat0, lon0, h0).
    
    :param x: meters
    :param y: meters
    :param z: meters
    :param lat0: latitude in degrees
    :param lon0: longitude in degrees
    :param height0: meters
    :returns: (xEast, yNorth, zUp) in [meters, meters, meters]
    """
    lamb = deg2rad(lat0)
    phi = deg2rad(lon0)
    s = sin(lamb)
    N = a / sqrt(1 - e_sq * s * s)

    sin_lambda = sin(lamb)
    cos_lambda = cos(lamb)
    sin_phi = sin(phi)
    cos_phi = cos(phi)

    x0 = (h0 + N) * cos_lambda * cos_phi
    y0 = (h0 + N) * cos_lambda * sin_phi
    z0 = (h0 + (1 - e_sq) * N) * sin_lambda

    xd = x - x0
    yd = y - y0
    zd = z - z0

    xEast = -sin_phi * xd + cos_phi * yd
    yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd
    zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd

    return xEast, yNorth, zUp


def enu_to_ecef(xEast:float, yNorth:float, zUp:float, lat0:float, lon0:float, h0:float) -> Tuple[float, float, float]:
    """
    Inverse of ecef_to_enu. Converts East-North-Up coordinates (xEast, yNorth, zUp) in a
    Local Tangent Plane that is centered at the (WGS-84) Geodetic point (lat0, lon0, h0)
    to the Earth-Centered Earth-Fixed (ECEF) coordinates (x, y, z).
    
    :param xEast: meters
    :param yNorth: meters
    :param zUp: meters
    :param lat0: latitude in degrees
    :param lon0: longitude in degrees
    :param height0: meters
    :returns: (x, y, z) in [meters, meters, meters]
    """
    lambd = deg2rad(lat0)
    phi = deg2rad(lon0)
    s = sin(lambd)
    N = a / sqrt(1 - e_sq * s * s)

    sin_lambda = sin(lambd)
    cos_lambda = cos(lambd)
    cos_phi = cos(phi)
    sin_phi = sin(phi)

    x0 = (h0 + N) * cos_lambda * cos_phi
    y0 = (h0 + N) * cos_lambda * sin_phi
    z0 = (h0 + (1 - e_sq) * N) * sin_lambda

    xd = -sin_phi * xEast - cos_phi * sin_lambda * yNorth + cos_lambda * cos_phi * zUp
    yd = cos_phi * xEast - sin_lambda * sin_phi * yNorth + cos_lambda * sin_phi * zUp
    zd = cos_lambda * yNorth + sin_lambda * zUp

    x = xd + x0
    y = yd + y0
    z = zd + z0

    return x, y, z

    
def geodetic_to_enu(lat:float, lon:float, h:float, lat0:float, lon0:float, h0:float) -> Tuple[float, float, float]:
    """
    Converts the geodetic WGS-84 coordinated (lat, lon, h) to 
    East-North-Up coordinates in a Local Tangent Plane that is centered at the 
    (WGS-84) Geodetic point (lat0, lon0, h0).
    
    :param lat: latitude of target in degrees
    :param lon: longitude of target in degrees
    :param h: height of target in meters
    :param lat0: latitude of origin in degrees
    :param lon0: longitude of origin in degrees
    :param h0: height of origin meters
    :returns: (xEast, yNorth, zUp) in [meters, meters, meters]
    """
    x, y, z = geodetic_to_ecef(lat, lon, h)
    xEast, yNorth, zUp = ecef_to_enu(x, y, z, lat0, lon0, h0)

    return xEast, yNorth, zUp


def sacrs_to_geodetic(loMeridian: int, yWesting: int, xSouthing: int) -> Tuple[float, float]:
    """
    South African Coordinate Reference System (Hartebeesthoek94) to Geodetic
    From "CDNGI Coordinate Conversion Utility v1 Sep 2009.xls".
    CM = central meridian

    :param loMeridian: central meridian (Lo.) degrees
    :param yWesting: coordinates measured in meters from the CM of the respective zone, increasing from the CM (where Y=0) in a westerly direction. Y is +ve west of the CM and -ve east of the CM.
    :param xSouthing: coordinates measured in meters southwards from the equator, increasing from the equator (where X = 0m) towards the South Pole.
    :returns: (lat, lon) in [degree, degree]
    """
    loMeridianRadians = deg2rad(loMeridian)
    ec_sq = (a_sq - b_sq) / a_sq   # e^2   in "CDNGI Coordinate Conversion Utility v1 Sep 2009.xls"
    ep_sq = (a_sq - b_sq) / b_sq   # e'^2
    n = (a - b) / (a + b)
    n_2 = n * n
    n_3 = n_2 * n
    n_4 = n_2 * n_2
    n_5 = n_2 * n_3
    p2 = 3.0 / 2.0 * n - 27.0 / 32.0 * n_3 + 269.0 / 512.0 * n_5
    p4 = 21.0 / 16.0 * n_2 - 55.0 / 32.0 * n_4
    p6 = 151.0 / 96.0 * n_3 - 417.0 / 128.0 * n_5
    p8 = 1097.0 / 512.0 * n_4
    p10 = 8011.0 / 2560.0 * n_5
    a0 = 1.0 / (n + 1.0) * (1.0 + 1.0 / 4.0 * n_2 + 1.0 / 64.0 * n_4)
    footBar = xSouthing / (a * a0)
    foot = footBar + p2 * sin(2.0 * footBar) + p4 * sin(4.0 * footBar) + p6 * sin(6.0 * footBar) + p8 * sin(8.0 * footBar) + p10 * sin(10.0 * footBar)
    Nf = a / sqrt(1.0 - ec_sq * sin(foot)**2)

    b1 = 1.0 / (Nf * cos(foot))
    b2 = tan(foot) / (2.0 * Nf * Nf * cos(foot))
    b3 = (1.0 + 2.0 * tan(foot)**2 + ep_sq * cos(foot)**2) / (6.0 * Nf**3 * cos(foot))
    b4 = (tan(foot) * (5.0 + 6.0 * tan(foot)**2 + ep_sq * cos(foot)**2)) / (24.0 * Nf**4 * cos(foot))
    b5 = (5.0 + 28.0 * tan(foot)**2 + 24.0 * tan(foot)**4) / (120.0 * Nf**5 * cos(foot))
    d1 = cos(foot) * (1.0 + ep_sq * cos(foot)**2)
    d2 = -1.0 / 2.0 * cos(foot)**2 * tan(foot) * (1.0 + 4.0 * ep_sq * cos(foot)**2)

    latRadians = -(foot - b2 * d1 * yWesting**2 + (b4 * d1 + b2 * b2 * d2) * yWesting**4)
    lonRadians = loMeridianRadians - (b1 * yWesting - b3 * yWesting**3 + b5 * yWesting**5)
    lat = rad2deg(latRadians)
    lon = rad2deg(lonRadians)

    return lat, lon


def geodetic_to_sacrs(lat: float, lon: float) -> Tuple[float, float, float]:
    """
    South African Coordinate Reference System (Hartebeesthoek94) to Geodetic
    From "CDNGI Coordinate Conversion Utility v1 Sep 2009.xls"

    :param lat: latitude in degrees
    :param lon: longitude in degrees
    :returns: (loMeridian, yWesting, xSouthing) in [degree, meters, meters]
    """
    # longitude of origin
    loMeridian = 2 * int(lon/2) + sign(lon)  # Take integer part of lon, rounded up to nearest odd number (or down if negative), so =ODD(TRUNC(lon))
    loMeridianRadians = deg2rad(loMeridian)
    ec_sq = (a_sq - b_sq) / a_sq   # e^2   in "CDNGI Coordinate Conversion Utility v1 Sep 2009.xls"
    ep_sq = (a_sq - b_sq) / b_sq   # e'^2
    n = (a - b) / (a + b)
    n_2 = n * n
    n_3 = n_2 * n
    n_4 = n_2 * n_2
    n_5 = n_2 * n_3
    A0 = 1.0 / (n + 1.0) * (1.0 + 1.0 / 4.0 * n_2 + 1.0 / 64.0 * n_4)
    A2 = -1.0 / (n + 1.0) * (3.0 / 2.0 * n - 3.0 / 16.0 * n_3 - 3.0 / 128.0 * n_5)
    A4 = 1.0 / (n + 1.0) * (15.0 / 16.0 * n_2 - 15.0 / 64.0 * n_4)
    A6 = -1.0 / (n + 1.0) * (35.0 / 48.0 * n_3 - 175.0 / 768.0 * n_5)
    A8 = 1.0 / (n + 1.0) * (315.0 / 512.0 * n_4)
    A10 = 1.0 / (n + 1.0) * (693.0 / 1280.0 * n_5)

    latRadians = deg2rad(-lat)
    lonRadians = deg2rad(lon)
    G = a * (A0 * latRadians + A2 * sin(2 * latRadians) + A4 * sin(4 * latRadians) + A6 * sin(6 * latRadians) + A8 * sin(8 * latRadians) + A10 * sin(10 * latRadians))
    N = a / sqrt(1.0 - ec_sq * sin(latRadians) * sin(latRadians))

    latCos = cos(latRadians)
    latCos_2 = latCos * latCos
    latCos_3 = latCos_2 * latCos
    latCos_4 = latCos_2 * latCos_2
    latCos_5 = latCos_4 * latCos
    latTan = tan(latRadians)
    latTan_2 = latTan * latTan
    latTan_4 = latTan_2 * latTan_2
    a1 = N * latCos
    a2 = -1.0 / 2.0 * N * latCos_2 * latTan
    a3 = -1.0 / 6.0 * N * latCos_3 * (1.0 - latTan_2 + ep_sq * latCos_2)
    a4 = 1.0 / 24.0 * N * latCos_4 * latTan * (5 - latTan_2 + 9.0 * ep_sq * latCos_2)
    a5 = 1.0 / 120.0 * N * latCos_5 * (5.0 - 18.0 * latTan_2 + latTan_4)
    l = lonRadians - loMeridianRadians
    l_2 = l * l
    l_3 = l * l_2
    l_4 = l_2 * l_2
    l_5 = l_2 * l_3
    xSouthing = G - a2 * l_2 + a4 * l_4
    yWesting = -(a1 * l - a3 * l_3 + a5 * l_5)

    return xSouthing, yWesting
