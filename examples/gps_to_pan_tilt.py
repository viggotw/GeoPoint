import numpy as np
import matplotlib.pyplot as plt
from geopoint.pan_tilt import geodetic_to_pan_tilt

# CONSTAINTS
TARGET_LAT, TARGET_H = 1e-4, 10  # 0.00001 deg ~= 1.11 m
ORIGIN_LAT, ORIGIN_LON, ORIGIN_H = 0, 0, 0  # Easy position at the equator
NORTH_OFFSET = 0

# Find the coordinates for the TARGET in a ENU coordinate system
# centered in a given origin, using the geodetic coordinates for
# both the TARGET and the ORIGIN
x = list(np.arange(-1e-3, 1e-3, 1e-5))
pans = []
tilts = []
for target_lon in x:
       pan, tilt = geodetic_to_pan_tilt(
              TARGET_LAT, target_lon, TARGET_H,
              ORIGIN_LAT, ORIGIN_LON, ORIGIN_H,
              NORTH_OFFSET
       )
       #print(f"TARGET lat: {TARGET_LAT:.2e}, lon: {target_lon:.2e}, h: {TARGET_H:.2e} ---> pan: {pan:.2f}, tilt: {tilt:.2f}")
       pans.append(pan)
       tilts.append(tilt)

plt.plot(x, pans)
plt.plot(x, tilts)
plt.legend(['pan', 'tilt'])
plt.xlabel('Target longitude')
plt.ylabel('Degrees')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.show()