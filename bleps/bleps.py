# Binary ecLipsE Prediction Software (BLEPS)
# v 1.0 (April 2025)
# Nikki Miller - nikkimillerastro@gmail.com

from astropy.time import Time
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
from astroplan import Observer, moon_illumination
import yaml


class Target:
    def __init__(self, name, ra, de, t0, p, w1, ph2=None, w2=None):
        """
        Ingest the information about the eclipsing target(s)

        :param name: np.array(str) - name of system(s)
        :param ra: np.array(str) - right ascension of target in hour angle format 'XX XX XX.XX'
        :param de: np.array(str) - declination of target in degree format '+/-XX XX XX.X'
        :param t0: np.array(float) - time of primary eclipse minimum in HJD or BJD
        :param p: np.array(float) - orbital period of system in days
        :param w1: np.array(float) - width of primary eclipse in phase units
        :param ph2: np.array(float) - phase of secondary eclipse in phase units
        :param w2: np.array(float) - width of secondary eclipse in phase units
        """
        self.name = name
        self.t0 = t0
        self.p = p
        self.coord = SkyCoord(ra, de, unit=(u.hourangle, u.deg))

        t_pri_start = t0 - (0.5 * w1) * p
        if ph2 is not None and w2 is not None:
            t_sec_start = t0 + p * ph2 - (0.5 * w2) * p
            self.widths = [w1, w2]
            self.t_start = [t_pri_start, t_sec_start]
            self.ecl_type = [1, 2]
            self.labels = ["Primary", "Secondary"]
        else:
            self.widths = [w1]
            self.t_start = [t_pri_start]
            self.ecl_type = [1]
            self.labels = ["Primary"]


    @staticmethod
    def _sun_moon(observer_object, time):
        """
        Calculate moon % illumination, sunset and sunrise times given observer location and time.
        Prevents this intensive computation from being called unnecessarily.

        :param observer_object: astroplan.Observer object
        :param time: Time() object
        :return: moon illumination, sunset time, sunrise time
        """
        moon = f"{round(moon_illumination(Time(time)) * 100)}%"
        sunset_time = str(observer_object.sun_set_time(Time(time), which='nearest').iso)
        sunrise_time = str(observer_object.sun_rise_time(Time(time), which='nearest').iso)
        return moon, sunset_time, sunrise_time


    def predict_eclipses(self, observer_object, search_start, search_end, alt_limits=(30, 80), strategy='split'):
        """
        Predict visible eclipses for the provided target(s) during the specified search time period.

        :param observer_object: astroplan.Observer object containing information about the observatory / site.
        :param search_start: Time() object containing the start of the search window for eclipses
        :param search_end: Time() object containing the end of the search window for eclipses
        :param alt_limits: Altitude limits in degrees
        :param strategy: How to approach returning viable eclipses.
            - 'split' (default):
                Consider each half of an eclipse separately in feasibility checks (ingress / egress).
                This is useful for long-period and/or totally eclipsing systems.
            - 'full':
                The whole eclipse must be pass feasibility checks in order to be returned.
            - 'any':
                The observation is returned as viable if any the system is observable at any point during the eclipse.
        :return: astropy.Table containing details about the eclipses which meet observability criteria.
        """

        obs = observer_object
        min_alt, max_alt = alt_limits
        search_start, search_end = search_start.jd, search_end.jd

        # Create empty astropy.Table object to save viable eclipses to
        if strategy == 'split':
            colnames = ["Name", "EclType", "Half", "RA", "Dec", "Moon", "UTC_sunset", "UTC_sunrise",
                        "UTC_obs_start", "Alt_obs_start", "Az_obs_start", "UTC_obs_end", "Alt_obs_end", "Az_obs_end"]
            types = [str, int, str, str, str, str, str, str, str, float, float, str, float, float]
        elif strategy == 'full' or strategy == 'any':
            colnames = ["Name", "EclType", "RA", "Dec", "Moon", "UTC_sunset", "UTC_sunrise",
                        "UTC_start", "Alt_start", "Az_start", "UTC_mid", "Alt_mid", "Az_mid", "UTC_end", "Alt_end", "Az_end"]
            types = [str, int, str, str, str, str, str, str, float, float, str, float, float, str, float, float]
        else:
            raise ValueError(f"Invalid strategy specified: {strategy}. Should be ['split', 'full', 'any']")

        save = Table(names=colnames, dtype=types)

        # Loop through each system
        for j in range(len(self.labels)):
            for i, _ in enumerate(self.t0):
                # Find first eclipse time in window
                t_next = self.t_start[j][i]
                while t_next < search_start:
                    t_next += self.p[i]

                # Loop through all eclipses in window
                while t_next < search_end:
                    # Predict the start, mid and end times of eclipses
                    obstime = Time(t_next, format='jd')
                    midtime = Time((t_next + 0.5 * self.p[i] * self.widths[j][i]), format='jd')
                    endtime = Time((t_next + self.p[i] * self.widths[j][i]), format='jd')

                    # Construct alt-az frame for each of these times
                    aa_frame = obs.altaz(Time([obstime, midtime, endtime]))
                    azalt = self.coord[i].transform_to(aa_frame)

                    # Observability checks
                    ok_start = min_alt < azalt[0].alt.value < max_alt and obs.sun_altaz(Time(obstime)).alt <= -12 *u.deg
                    ok_mid = min_alt < azalt[1].alt.value < max_alt and obs.sun_altaz(Time(midtime)).alt <= -12 * u.deg
                    ok_end = min_alt < azalt[2].alt.value < max_alt and obs.sun_altaz(Time(endtime)).alt <= -12 * u.deg

                    if strategy == 'split':
                        # Split eclipses into ingress/egress halves – necessary for longer period eclipsing systems
                        ok_ingress = ok_start and ok_mid
                        ok_egress = ok_mid and ok_end

                        if ok_ingress:
                            moon, sunset, sunrise = self._sun_moon(obs, obstime)
                            save.add_row([self.name[i], self.ecl_type[j], 'ingress', str(self.coord.ra[i]),
                                          str(self.coord.dec[i]), moon, sunset, sunrise,
                                          obstime.iso, round(azalt[0].alt.value, 1), round(azalt[0].az.value, 1),
                                          midtime.iso, round(azalt[1].alt.value, 1), round(azalt[1].az.value, 1)])

                        elif ok_egress:
                            moon, sunset, sunrise = self._sun_moon(obs, midtime)
                            save.add_row([self.name[i], self.ecl_type[j], 'egress', str(self.coord.ra[i]),
                                          str(self.coord.dec[i]), moon, sunset, sunrise,
                                          midtime.iso, round(azalt[1].alt.value, 1), round(azalt[1].az.value, 1),
                                          endtime.iso, round(azalt[2].alt.value, 1), round(azalt[2].az.value, 1)])

                    elif strategy == 'full':
                        # The entire eclipse must be observable
                        if ok_start and ok_mid and ok_end:
                            moon, sunset, sunrise = self._sun_moon(obs, obstime)
                            save.add_row([self.name[i], self.ecl_type[j], str(self.coord.ra[i]), str(self.coord.dec[i]),
                                          moon, sunset, sunrise,
                                          obstime.iso, round(azalt[0].alt.value, 1), round(azalt[0].az.value, 1),
                                          midtime.iso, round(azalt[1].alt.value, 1), round(azalt[1].az.value, 1),
                                          endtime.iso, round(azalt[2].alt.value, 1), round(azalt[2].az.value, 1)])

                    elif strategy == 'any':
                        # Any part of the eclipse is observable
                        if ok_start or ok_mid or ok_end:
                            moon, sunset, sunrise = self._sun_moon(obs, obstime)
                            save.add_row([self.name[i], self.ecl_type[j], self.coord.ra[i], self.coord.dec[i],
                                         moon, sunset, sunrise,
                                         obstime.iso, round(azalt[0].alt.value, 1), round(azalt[0].az.value, 1), ok_start,
                                         midtime.iso, round(azalt[1].alt.value, 1), round(azalt[1].az.value, 1), ok_mid,
                                         endtime.iso, round(azalt[2].alt.value, 1), round(azalt[2].az.value, 1), ok_end])

                    t_next += self.p[i]
        return save


def gcal_format(table):
    """
    Reformat the eclipse timing info to be readable by Google Calendar.

    :param table: Astropy table generated by the Target.predict_eclipses method
    :type table: astropy.table.Table
    :return: Reformatted table ready to save as .csv
    """
    names = ['Subject', 'Start Date', 'Start Time', 'End Date', 'End Time', 'All Day Event', 'Description']
    gcal = Table(names=names, dtype=[str]*len(names))
    for row in table:
        subject = f"{row['Name']} primary eclipse" if row['EclType'] == 1 else f"{row['Name']} secondary eclipse"
        try:
            ecl_start = row['UTC_obs_start']
            ecl_end = row['UTC_obs_end']
            altaz_info = (f"Target altitude (start/end): {row["Alt_obs_start"]}, {row["Alt_obs_end"]}. "
                          f"Azimuth: {row["Az_obs_start"]}, {row["Az_obs_end"]}. ")
            io = row['Half']
        except KeyError:
            ecl_start = row['UTC_start']
            ecl_end = row['UTC_end']
            altaz_info = (f"Target altitude (start/mid/end): {row["Alt_start"]}, {row["Alt_mid"]}, {row["Alt_end"]}. "
                          f"Azimuth: {row["Az_start"]}, {row["Az_mid"]}, {row["Az_end"]}. ")
            io = 'entire eclipse'
        d_start = f"{ecl_start[8:10]}/{ecl_start[5:7]}/{ecl_start[0:4]}"
        t_start = f"{ecl_start[11:13]}:{ecl_start[14:16]}"
        d_end = f"{ecl_end[8:10]}/{ecl_end[5:7]}/{ecl_end[0:4]}"
        t_end = f"{ecl_end[11:13]}:{ecl_end[14:16]}"
        all_day = 'FALSE'
        description = (f"This is an {io} starting at {t_start} UTC - careful with timezones! "
                       f"Moon: {row['Moon']}, sunset: {row['UTC_sunset'][11:16]}, "
                       f"sunrise: {row['UTC_sunrise'][11:16]}. {altaz_info}")

        gcal.add_row([subject, d_start, t_start, d_end, t_end, all_day, description])

    return gcal


if __name__ == '__main__':
    print("====================================================")
    print("||   Binary ecLipsE Prediction Software (BLEPS)   ||")
    print("||              v.1.0 (April 2025)                ||")
    print("||   Nikki Miller - nikkimillerastro@gmail.com    ||")
    print("====================================================")

    # Load configuration file
    c = yaml.load(open('config.yaml', 'r'), Loader=yaml.SafeLoader)

    # Set up observer object and altitude limits
    o = c['observatory']
    print(f"\nSetting up {o['name']} site: {round(o['lat'], 3)} deg, {round(o['lon'], 3)} deg, +{o['elev']} m")
    observer = Observer(latitude=o['lat'] * u.deg, longitude=o['lon'] * u.deg, elevation=o['elev'] * u.m,
                        pressure=1 * u.bar, name=o['name'], timezone='UTC')
    limits = (c['altitude_limits'][0], c['altitude_limits'][1])  # Altitude limits in degrees

    # Search window definition
    t1 = Time.now()
    t2 = t1 + c['search'] * u.day
    start, end = sorted([t1, t2])

    # Load target file and set up Target object
    plan = Table.read(c['target_file'], format='csv')
    print(f"{len(plan)} eclipsing binary target(s) read from '{c['target_file']}'\n")

    target_name = np.array(plan['Name'])  # str - name of system(s)
    r_a = np.array(plan['RA'])            # str - right ascension of target in hour angle format 'XX XX XX.XX'
    dec = np.array(plan['Dec'])           # str - declination of target in degree format '+/-XX XX XX.X'
    t_zero = np.array(plan['T0'])         # float - time of primary eclipse minimum in HJD or BJD
    period = np.array(plan['Period'])     # float - orbital period of system in days
    width_1 = np.array(plan['Width1'])    # float - width of primary eclipse in phase units
    phase_2 = np.array(plan['Phase2'])    # float - phase of secondary eclipse in phase units
    width_2 = np.array(plan['Width2'])    # float - width of secondary eclipse in phase units

    t = Target(target_name, r_a, dec, t_zero, period, width_1, phase_2, width_2)

    # Predict upcoming eclipses
    print(f"Searching for eclipses between {start.iso[0:10]} and {end.iso[0:10]}, using '{c['strategy']}' strategy")
    s = t.predict_eclipses(observer, start, end, strategy=c['strategy'])
    print(f"Found {len(s)} eclipse(s).\n")

    # Save info on upcoming eclipses to file (or print to terminal)
    file_name = f"{o['name']}_eclipses_{Time(start, format='jd').iso[0:10]}_to_{Time(end, format='jd').iso[0:10]}.csv"
    if c['save']:
        print(f"Saving details to file: '{file_name}'.")
        s.write(file_name, format='csv', overwrite=True)
    else:
        print(s)

    if c['gcal']:
        g = gcal_format(s)
        file_name = 'gcal_' + file_name
        print(f"Saving details to Google Calendar-compatable format: '{file_name}'.")
        g.write(file_name, format='csv', overwrite=True)
