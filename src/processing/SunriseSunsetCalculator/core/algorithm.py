# Algorithm Source: https://edwilliams.org/sunrise_sunset_example.htm
# 	Almanac for Computers, 1990
# 	published by Nautical Almanac Office
# 	United States Naval Observatory

# day, month, year:      date of sunrise/sunset
# 	latitude, longitude:   location for sunrise/sunset
# 	zenith:                Sun's zenith for sunrise/sunset
# 	  offical      = 90 degrees 50'
# 	  civil        = 96 degrees
# 	  nautical     = 102 degrees
# 	  astronomical = 108 degrees

# 	NOTE: longitude is positive for East and negative for West

from math import floor, sin, atan, tan, cos, asin, acos, pi

zenith_angle = {'official': 90 + 50/60,
                'civil': 96,
                'nautical': 102,
                'astronomical': 108}


def get_day_of_Year(date_info):
    day = date_info.day
    month = date_info.month
    year = date_info.year
    # Explained here: https://astronomy.stackexchange.com/questions/2407/calculate-day-of-the-year-for-a-given-date
    N1 = floor(275 * month / 9)
    N2 = floor((month + 9) / 12)
    N3 = (1 + floor((year - 4 * floor(year / 4) + 2) / 3))
    N = N1 - (N2 * N3) + day - 30
    return N


def AlmanacSunrise(latitude, longitude, date_info, sunrise_sunset_type):
    # zenith angle- 90 degrees 50'
    zenith = pi * zenith_angle.get(sunrise_sunset_type) / 180
    # Step 1. calculate the day of the year
    N = get_day_of_Year(date_info)

    # Step 2. convert the longitude to hour value and calculate an approximate time

    lngHour = longitude / 15

    # if rising time is desired:
    t = N + ((6 - lngHour) / 24)

    # 3. calculate the Sun's mean anomaly
    M = (0.9856 * t) - 3.289

    # 4. calculate the Sun's true longitude
    L = M + (1.916 * sin(M * pi / 180)) + \
        (0.020 * sin(2 * M * pi / 180)) + 282.634

    if L > 360:
        L -= 360
    elif L < 0:
        L += 360

    # Step 5a. calculate the Sun's right ascension
    RA = (180 / pi) * atan(0.91764 * tan(L * pi / 180))

    # NOTE: RA potentially needs to be adjusted into the range [0,360) by adding/subtracting 360
    if RA > 360:
        RA -= 360
    elif RA < 0:
        RA += 360

    # right ascension value needs to be in the same quadrant as L

    Lquadrant = (floor(L / 90)) * 90
    RAquadrant = (floor(RA / 90)) * 90
    RA = RA + (Lquadrant - RAquadrant)

    # right ascension value needs to be converted into hours

    RA = RA / 15

    # calculate the Sun's declination

    sinDec = 0.39782 * sin(L * pi / 180)

    cosDec = cos(asin(sinDec))

    # calculate the Sun's local hour angle

    cosH = (cos(zenith) - (sinDec * sin(latitude * pi / 180))) / \
        (cosDec * cos(latitude * pi / 180))

    if cosH > 1:
        print("the sun never rises on this location (on the specified date)")
        return -2
    if cosH < -1:
        print("the sun never sets on this location (on the specified date)")
        return -1
    # finish calculating H and convert into hours

    # if if rising time is desired:
    H = 360 - acos(cosH) * 180 / pi
    # if setting time is desired:
    #   H = acos(cosH)

    H = H / 15

    # calculate local mean time of rising/setting

    T = H + RA - (0.06571 * t) - 6.622

    # adjust back to UTC

    UTC = T - lngHour
    # NOTE: UT potentially needs to be adjusted into the range [0,24) by adding/subtracting 24
    if UTC < 0:
        UTC += 24
    elif UTC > 24:
        UTC -= 24
    return UTC


def AlmanacSunset(latitude, longitude, date_info, sunrise_sunset_type):
    zenith = pi * zenith_angle.get(sunrise_sunset_type) / 180

    # Step 1. calculate the day of the year
    N = get_day_of_Year(date_info)
    # Step 2. convert the longitude to hour value and calculate an approximate time

    lngHour = longitude / 15

    # if rising time is desired:
    #     t = N + ((6 - lngHour) / 24)
    # if setting time is desired:
    t = N + ((18 - lngHour) / 24)

    # 3. calculate the Sun's mean anomaly
    M = (0.9856 * t) - 3.289

    # 4. calculate the Sun's true longitude
    L = M + (1.916 * sin(M * pi / 180)) + \
        (0.020 * sin(2 * M * pi / 180)) + 282.634

    if L > 360:
        L -= 360
    elif L < 0:
        L += 360

    # Step 5a. calculate the Sun's right ascension
    RA = (180 / pi) * atan(0.91764 * tan(L * pi / 180))

    # NOTE: RA potentially needs to be adjusted into the range [0,360) by adding/subtracting 360
    if RA > 360:
        RA -= 360
    elif RA < 0:
        RA += 360

    # right ascension value needs to be in the same quadrant as L

    Lquadrant = (floor(L / 90)) * 90
    RAquadrant = (floor(RA / 90)) * 90
    RA = RA + (Lquadrant - RAquadrant)

    # right ascension value needs to be converted into hours

    RA = RA / 15

    # calculate the Sun's declination

    sinDec = 0.39782 * sin(L * pi / 180)

    cosDec = cos(asin(sinDec))

    # calculate the Sun's local hour angle

    cosH = (cos(zenith) - (sinDec * sin(latitude * pi / 180))) / \
        (cosDec * cos(latitude * pi / 180))

    if cosH > 1:
        print("the sun never rises on this location (on the specified date)")
        return -1
    if cosH < -1:
        print("the sun never sets on this location (on the specified date)")
        return -2
    # finish calculating H and convert into hours

    # if if rising time is desired:
    #     H = 360 - acos(cosH) * 180 / pi
    # if setting time is desired:
    H = acos(cosH) * 180 / pi

    H = H / 15

    # calculate local mean time of rising/setting

    T = H + RA - (0.06571 * t) - 6.622

    # adjust back to UTC

    UTC = T - lngHour
    # NOTE: UT potentially needs to be adjusted into the range [0,24) by adding/subtracting 24

    if UTC < 0:
        UTC += 24
    elif UTC > 24:
        UTC -= 24
    return UTC
