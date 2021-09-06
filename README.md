# GeoPoint: Essential transformations for PTZ units and conversion between global coordinate systems

## 1. Easily convert between different coordinate systems
There are many different coordinate systems out there, each tailored for different usecases. GeoPoint has the tools for converting between several of the common ones, like:
- WGS-84 Geodetic (latitude, longitude, height)
- Earth-Centered Earth-Fixed (ECEF)
- East-North-Up (ENU)
- South African Coordinate Reference System (SACRS)

## 2. Find pan/tilt angles necessary for looking at a specific lat/lon location
Given a pan-tilt-zoom (PTZ) or general pointing unit (PU), GeoPoint offers the tools for calculating the pan/tilt angles needed for a PU to look directly at a target with given lat/lon coordinates. This can also be used for tracking objects where either the PU, the target, or both are moving.

```
pan, tilt = geodetic_to_pan_tilt(
       TARGET_LAT, TARGET_LON, TARGET_H,
       ORIGIN_LAT, ORIGIN_LON, ORIGIN_H,
       NORTH_OFFSET
)
```

## 3. Compensate for offsets in pitch and roll
For PTZ units with powerful tele, GeoPoint offers tools for calculating the pan/tilt angles while also compensating for an offset in pitch and roll. This is usefull when the PTZ is not level.

```
pan, tilt = geodetic_to_pan_tilt_with_offsets(
       TARGET_LAT, TARGET_LON, TARGET_H,
       ORIGIN_LAT, ORIGIN_LON, ORIGIN_H,
       NORTH_OFFSET, PITCH_OFFSET, ROLL_OFFSET
)
```
