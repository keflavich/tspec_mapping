import datetime

def dateobs_to_datetime(fitsstr):
    """
    Converte a DATE-OBS time string to a datetime instance

    string format should be like:
    '2012-11-05T08:47:31.549'

    (dateutil.parser.parse is better)
    """
    year = int(fitsstr[:4])
    month = int(fitsstr[5:7])
    day = int(fitsstr[8:10])
    hour = int(fitsstr[11:13])
    minute = int(fitsstr[14:16])
    second = int(fitsstr[17:19])
    microsecond = int(fitsstr[20:])
    return datetime.datetime(year,month,day,hour,minute,second,microsecond)


