import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
!pip install skyfield
from skyfield.api import Topos, load
from skyfield import almanac
from datetime import datetime, timezone, timedelta
!pip install timezonefinder pytz
from timezonefinder import TimezoneFinder
import pytz

# Define observatory names
observatory_names = [
    'Palomar Observatory',
    'S.D.S.U. Observatory',
    'Table Mountain Observatory',
    'Gemini Observatory',
    'Keck Observatory',
    'Ishigakijima Astronomical Observatory',
    'Lulin Observatory',
    'Himalayan Chandra Observatory',
    'Girawali Observatory',
    'Wise Observatory',
    'Stella Robotic Observatory',
    'Nordic Optical Telescope',
    'Discovery Channel Telescope',
    'Fenton Hill Observatory',
    'The Liverpool Telescope',
    'GROWTH-India Telescope'
];

# Observatory coordinates (latitude and longitude in degrees)
observatory_coords = [
    (33.3564, -116.8650),  # Palomar Observatory
    (32.8422, -116.4276),  # S.D.S.U. Observatory
    (34.3820, -117.6818),  # Table Mountain Observatory
    (19.8238, -155.4689),  # Gemini Observatory
    (20.0243, -155.6655),  # Keck Observatory
    (24.3733, 124.1358),   # Ishigakijima Astronomical Observatory
    (23.4693, 120.8727),   # Lulin Observatory
    (32.7795, 78.9642),    # Himalayan Chandra Observatory
    (19.0725, 73.8452),    # Girawali Observatory
    (30.6102, 34.8019),    # Wise Observatory
    (28.3011, -16.5092),   # Stella Robotic Observatory
    (28.7572, -17.8851),   # Nordic Optical Telescope
    (34.7444, -111.4223),  # Discovery Channel Telescope
    (35.8813, -106.6734),  # Fenton Hill Observatory
    (28.7624, -17.8792),   # The Liverpool Telescope
    (32.7793, 78.9646)     # GROWTH-India Telescope
];

# Define the time of the GRB detection with astropy time object
utc_time = Time('2024-08-29T17:00:16', format='isot', scale='utc');

# Define the GRB coordinates
grb_coord = SkyCoord(ra=20.95*u.hour, dec=-35.7*u.deg, frame='icrs');

# Defining the observatory locations
observatory_locations = [
    EarthLocation(lat=33.3564*u.deg, lon=-116.8650*u.deg),  # Palomar Observatory
    EarthLocation(lat=32.8422*u.deg, lon=-116.4276*u.deg),  # S.D.S.U. Observatory
    EarthLocation(lat=34.3820*u.deg, lon=-117.6818*u.deg),  # Table Mountain Observatory
    EarthLocation(lat=19.8238*u.deg, lon=-155.4689*u.deg),  # Gemini Observatory
    EarthLocation(lat=20.0243*u.deg, lon=-155.6655*u.deg),  # Keck Observatory
    EarthLocation(lat=24.3733*u.deg, lon=124.1358*u.deg),   # Ishigakijima Astronomical Observatory
    EarthLocation(lat=23.4693*u.deg, lon=120.8727*u.deg),   # Lulin Observatory
    EarthLocation(lat=32.7795*u.deg, lon=78.9642*u.deg),    # Himalayan Chandra Observatory
    EarthLocation(lat=19.0725*u.deg, lon=73.8452*u.deg),    # Girawali Observatory
    EarthLocation(lat=30.6102*u.deg, lon=34.8019*u.deg),    # Wise Observatory
    EarthLocation(lat=28.3011*u.deg, lon=-16.5092*u.deg),   # Stella Robotic Observatory
    EarthLocation(lat=28.7572*u.deg, lon=-17.8851*u.deg),   # Nordic Optical Telescope
    EarthLocation(lat=34.7444*u.deg, lon=-111.4223*u.deg),  # Discovery Channel Telescope
    EarthLocation(lat=35.8813*u.deg, lon=-106.6734*u.deg),  # Fenton Hill Observatory
    EarthLocation(lat=28.7624*u.deg, lon=-17.8792*u.deg),   # The Liverpool Telescope
    EarthLocation(lat=32.7793*u.deg, lon=78.9646*u.deg)     # GROWTH-India Telescope
];

# Initialize a NumPy array to store Altitude and Azimuth for each observatory
altaz_array = np.zeros((len(observatory_locations), 2)) ;

# Loop through each observatory and calculate the Alt/Az for the GRB
for i, observatory in enumerate(observatory_locations):
    # Create AltAz frame for the observatory then transform the GRB coordinates
    altaz_coord = grb_coord.transform_to(AltAz(obstime=utc_time, location=observatory));

    # Store Altitude and Azimuth in the array
    altaz_array[i, 0] = altaz_coord.alt.deg;  # Altitude (in degrees)
    altaz_array[i, 1] = altaz_coord.az.deg;   # Azimuth (in degrees)

# Print the resulting Altitude for each observatory
for i, observatory_name in enumerate(observatory_names):
    print(f"The altitude of the GRB from the {observatory_name} is "
          f"{altaz_array[i, 0]:.2f}° at the time of detection.");

    # Check if the GRB is above the horizon
    if altaz_array[i, 0] > 0:
        print(f"The GRB may be visible from here! The GRB is above the "
              f"horizon from this location, but is the Sun out? {observatory_name} "
              "advances to the next step in our evaluation criteria.\n");
    else:
        print("This means the GRB is not visible from here, since the "
              "location of the GRB is below the horizon.\n");

# List to store whether the GRB is observable at each observatory
observability = [];

# Loading in ephemeris data on the Sun
eph = load('de421.bsp');
sun = eph['sun'];

# Loading in timescale with Skyfield time object
ts = load.timescale();
t = ts.utc(2024, 8, 29, 17, 0, 16)

# Loop through each observatory
for i, (name, (lat, lon)) in enumerate(zip(observatory_names, observatory_coords)):
    # Define the observatory's location on Earth
    observatory = Topos(latitude_degrees=lat, longitude_degrees=lon);

    # Calculate the Sun's position
    astrometric = (eph['earth'] + observatory).at(t).observe(sun)

    # Get the apparent position (apply light-time correction, etc.)
    apparent = astrometric.apparent();

    # Extract position information
    alt, az, distance = apparent.altaz();

    # Below -18 degrees is common for defining nighttime in astronomy
    if (alt.degrees < -18):
      print(f"It is nighttime at the {name} and the Sun is {alt.degrees:.2f}°"
      " below the horizon.");
    else:
      print(f"It is daytime at the {name}");

print();
# PART 2

# Initialize the TimezoneFinder object
tf = TimezoneFinder();

# Defining function that finds UTC offset for specific coordinates
def get_utc_offset(latitude, longitude, midnight_local, tf):

    # Find the time zone based on latitude and longitude
    timezone_str = tf.timezone_at(lat=latitude, lng=longitude);

    if timezone_str is None:
        raise ValueError(f"Could not find time zone for the location: "
        f"lat={latitude}, lon={longitude}");

    # Get the local time zone object using pytz
    local_timezone = pytz.timezone(timezone_str);

    # Makes midnight_local a time-zone aware object and assigns new value
    local_aware = local_timezone.localize(midnight_local);

    # Get the UTC offset in seconds (may include DST)
    utc_offset_seconds = local_aware.utcoffset().total_seconds();

    # Convert the offset to a timedelta object
    utc_offset_timedelta = timedelta(seconds=utc_offset_seconds);

    return utc_offset_timedelta;

# Function to get Moon's phase and position
def get_moon_phase_and_position(lat, lon, time):

    # Time of observation
    t = ts.utc(time.year, time.month, time.day, time.hour, time.minute, time.second);

    # Moon's position relative to Earth observer
    observer = eph['earth'] + Topos(latitude_degrees=lat, longitude_degrees=lon);
    astrometric = observer.at(t).observe(eph['moon']);
    apparent = astrometric.apparent();

    # Extracting information
    alt, az, distance = apparent.altaz();

    # Moon's phase angle
    phase = almanac.moon_phase(eph, t);
    phase_degrees = phase.degrees % 360;
    # Returns value where 0° is New Moon, 90° is First Quarter,
    # 180° is Full Moon, and 270° is Last Quarter

    phase_name = '';
    # Determine the phase name
    if 0 <= phase_degrees < 45 or 315 <= phase_degrees <= 360:
        phase_name = 'New Moon';
    elif 45 <= phase_degrees < 90:
        phase_name = 'Waxing Crescent';
    elif 90 <= phase_degrees < 135:
        phase_name = 'First Quarter';
    elif 135 <= phase_degrees < 180:
        phase_name = 'Waxing Gibbous';
    elif 180 <= phase_degrees < 225:
        phase_name = 'Full Moon';
    elif 225 <= phase_degrees < 270:
        phase_name = 'Waning Gibbous';
    elif 270 <= phase_degrees < 315:
        phase_name = 'Last Quarter';
    else:
        phase_name = 'Waning Crescent';

    return alt.degrees, az.degrees, phase_degrees, phase_name;

# Time resolution in days
time_step = timedelta(days=1);

# Set the start and end dates for the search period (Set to midnight)
start_date = datetime(2024, 8, 29);
end_date = datetime(2025, 8, 29);

# HUGE LOOP GOING THROUGH EACH OBSERVATORY
for i, observatory in enumerate(observatory_locations):

  # NESTED LOOP GOING THROUGH EVERY DAY OF THE YEAR
  current_date = start_date;
  while current_date <= end_date:
      # Define midnight local time for variable clarity
      midnight_local = current_date;

      # Get latitude and longitude
      lat = observatory_coords[i][0]
      lon = observatory_coords[i][1]

      # Convert local time to UTC depending on time zone
      utc_offset = get_utc_offset(lat, lon, midnight_local, tf);

      # Calculate the utc time for local midnight (consider DST)
      midnight_utc = midnight_local - utc_offset;

      # Calculate the Local Sidereal Time at midnight
      # Create Astropy Time object in UTC
      time = Time(midnight_utc);

      # Calculate Local Sidereal Time at the telescope's longitude
      lst = time.sidereal_time('mean', longitude=observatory.lon);

      # Giving lst units
      lst_hours = lst.to_value(u.hourangle);

      # Check if LST matches the target RA (e.g., 20.95 hours) within a tolerance
      target_ra_hours = 20.95;  # Target Right Ascension in hours
      tolerance = 0.1;  # Tolerance in hours

      # Calculate the difference
      lst_diff = np.abs((lst_hours - target_ra_hours));

      if lst_diff < tolerance:

        # Proceed with moon calculations for days of transits
        # Calculate the moon's position and phase at time of observation
        moon_altitude, moon_azimuth, phase_degrees, phase_name = get_moon_phase_and_position(lat, lon, midnight_utc);

        # Calculate the GRB's position to compare
        # AltAz frame for the time and location
        altaz_frame = AltAz(obstime=midnight_utc, location=observatory);

        # Transform GRB coordinates to AltAz frame
        grb_altaz = grb_coord.transform_to(altaz_frame)
        grb_altitude = grb_altaz.alt.deg
        grb_azimuth = grb_altaz.az.deg

        # Compare to moons azimuth and position
        # Create SkyCoord objects for the Moon and GRB positions in the AltAz frame
        moon_coord = SkyCoord(alt=moon_altitude*u.deg, az=moon_azimuth*u.deg, frame=altaz_frame)
        grb_coord_altaz = SkyCoord(alt=grb_altitude*u.deg, az=grb_azimuth*u.deg, frame=altaz_frame)

        # Compute the angular separation
        separation = grb_coord_altaz.separation(moon_coord).deg

        # Determine if the Moon is close to the GRB
        threshold = 70  # degrees
        if separation < threshold:
            # Print information to console
            print(f"On {current_date.strftime('%Y-%m-%d')} at {observatory_names[i]} the GRB"
            f" transits the meridian at local midnight. At time of transit, the Moon is {separation:.2f}°"
            f" from the GRB. Potential interference. The moon's phase degree is {phase_degrees:.2f}°,"
            f" which corresponds to a {phase_name}\n");
        else:
            print(f"On {current_date.strftime('%Y-%m-%d')} at {observatory_names[i]}, "
            f"the Moon is {separation:.2f}° from the GRB. Observation favorable."
            f" The moon's phase degree is {phase_degrees:.2f}°, which "
            f"corresponds to a {phase_name}\n");

      # Proceed to next day
      current_date += time_step;
