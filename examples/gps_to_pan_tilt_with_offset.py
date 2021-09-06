from geopoint.pan_tilt import geodetic_to_pan_tilt_with_offsets

# CONSTAINTS
TARGET_LAT, TARGET_LON, TARGET_H = 1e-5, 1e-5, 10  # 0.00001 deg ~= 1.11 m
ORIGIN_LAT, ORIGIN_LON, ORIGIN_H = 0, 0, 0  # Easy position at the equator
NORTH_OFFSET, PITCH_OFFSET, ROLL_OFFSET = 90, 10, 10

# Find the coordinates for the TARGET in a ENU coordinate system
# centered in a given origin, using the geodetic coordinates for
# both the TARGET and the ORIGIN
pan, tilt = geodetic_to_pan_tilt_with_offsets(
       TARGET_LAT, TARGET_LON, TARGET_H,
       ORIGIN_LAT, ORIGIN_LON, ORIGIN_H,
       NORTH_OFFSET, PITCH_OFFSET, ROLL_OFFSET
)

print(f"pan: {pan:.2f}, tilt: {tilt:.2f}")