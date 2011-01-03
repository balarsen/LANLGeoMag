"""
Lgm_CTrans module, this contains the necessary code for coordinate
transformations in Lgm

the ugly top part is
Generated by h2py from ../libLanlGeoMag/Lgm/Lgm_CTrans.h
@todo: go through and remove the bits we don't need

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""


###############################################################################
###############################################################################
###############################################################################
LGM_GOLD = 1.61803398874989484820
LGM_1O_GOLD = 0.61803398874989484820
LGM_1M_1O_GOLD = 0.38196601125010515180
LGM_ERROR = -1
LGM_FILL_VALUE = (-1e31)
LGM_JD_J2000 = 2451545.0
LGM_JD_GPS0 = 2444245.0
LGM_JD_TAI0 = 2436205.0
LGM_TIME_SYS_UTC = 0
LGM_TIME_SYS_TAI = 1
LGM_TIME_SYS_GPS = 2
LGM_TIME_SYS_TT = 3
LGM_TIME_SYS_TDB = 4
LGM_TIME_SYS_UT1 = 5
EME2000_COORDS = 1
ICRF2000_COORDS = 1
GEI2000_COORDS = 1
MOD_COORDS = 2
TOD_COORDS = 3
TEME_COORDS = 4
PEF_COORDS = 5
WGS84_COORDS = 6
IRTF_COORDS = 6
GEO_COORDS = 6
GSE_COORDS = 7
GSM_COORDS = 8
SM_COORDS = 9
EDMAG_COORDS = 10
CDMAG_COORDS = 11
EME2000_TO_EME2000 = 101
EME2000_TO_ICRF2000 = 101
EME2000_TO_GEI2000 = 101
EME2000_TO_MOD = 102
EME2000_TO_TOD = 103
EME2000_TO_TEME = 104
EME2000_TO_PEF = 105
EME2000_TO_WGS84 = 106
EME2000_TO_ITRF = 106
EME2000_TO_GEO = 106
EME2000_TO_GSE = 107
EME2000_TO_GSM = 108
EME2000_TO_SM = 109
EME2000_TO_EDMAG = 110
EME2000_TO_CDMAG = 111
ICRF2000_TO_EME2000 = 101
ICRF2000_TO_ICRF2000 = 101
ICRF2000_TO_GEI2000 = 101
ICRF2000_TO_MOD = 102
ICRF2000_TO_TOD = 103
ICRF2000_TO_TEME = 104
ICRF2000_TO_PEF = 105
ICRF2000_TO_WGS84 = 106
ICRF2000_TO_ITRF = 106
ICRF2000_TO_GEO = 106
ICRF2000_TO_GSE = 107
ICRF2000_TO_GSM = 108
ICRF2000_TO_SM = 109
ICRF2000_TO_EDMAG = 110
ICRF2000_TO_CDMAG = 111
GEI2000_TO_EME2000 = 101
GEI2000_TO_ICRF2000 = 101
GEI2000_TO_GEI2000 = 101
GEI2000_TO_MOD = 102
GEI2000_TO_TOD = 103
GEI2000_TO_TEME = 104
GEI2000_TO_PEF = 105
GEI2000_TO_WGS84 = 106
GEI2000_TO_ITRF = 106
GEI2000_TO_GEO = 106
GEI2000_TO_GSE = 107
GEI2000_TO_GSM = 108
GEI2000_TO_SM = 109
GEI2000_TO_EDMAG = 110
GEI2000_TO_CDMAG = 111
MOD_TO_EME2000 = 201
MOD_TO_ICRF2000 = 201
MOD_TO_GEI2000 = 201
MOD_TO_MOD = 202
MOD_TO_TOD = 203
MOD_TO_TEME = 204
MOD_TO_PEF = 205
MOD_TO_WGS84 = 206
MOD_TO_ITRF = 206
MOD_TO_GEO = 206
MOD_TO_GSE = 207
MOD_TO_GSM = 208
MOD_TO_SM = 209
MOD_TO_EDMAG = 210
MOD_TO_CDMAG = 211
TOD_TO_EME2000 = 301
TOD_TO_ICRF2000 = 301
TOD_TO_GEI2000 = 301
TOD_TO_MOD = 302
TOD_TO_TOD = 303
TOD_TO_TEME = 304
TOD_TO_PEF = 305
TOD_TO_WGS84 = 306
TOD_TO_ITRF = 306
TOD_TO_GEO = 306
TOD_TO_GSE = 307
TOD_TO_GSM = 308
TOD_TO_SM = 309
TOD_TO_EDMAG = 310
TOD_TO_CDMAG = 311
TEME_TO_EME2000 = 401
TEME_TO_ICRF2000 = 401
TEME_TO_GEI2000 = 401
TEME_TO_MOD = 402
TEME_TO_TOD = 403
TEME_TO_TEME = 404
TEME_TO_PEF = 405
TEME_TO_WGS84 = 406
TEME_TO_ITRF = 406
TEME_TO_GEO = 406
TEME_TO_GSE = 407
TEME_TO_GSM = 408
TEME_TO_SM = 409
TEME_TO_EDMAG = 410
TEME_TO_CDMAG = 411
PEF_TO_EME2000 = 501
PEF_TO_ICRF2000 = 501
PEF_TO_GEI2000 = 501
PEF_TO_MOD = 502
PEF_TO_TOD = 503
PEF_TO_TEME = 504
PEF_TO_PEF = 505
PEF_TO_WGS84 = 506
PEF_TO_ITRF = 506
PEF_TO_GEO = 506
PEF_TO_GSE = 507
PEF_TO_GSM = 508
PEF_TO_SM = 509
PEF_TO_EDMAG = 510
PEF_TO_CDMAG = 511
WGS84_TO_EME2000 = 601
WGS84_TO_ICRF2000 = 601
WGS84_TO_GEI2000 = 601
WGS84_TO_MOD = 602
WGS84_TO_TOD = 603
WGS84_TO_TEME = 604
WGS84_TO_PEF = 605
WGS84_TO_WGS84 = 606
WGS84_TO_ITRF = 606
WGS84_TO_GEO = 606
WGS84_TO_GSE = 607
WGS84_TO_GSM = 608
WGS84_TO_SM = 609
WGS84_TO_EDMAG = 610
WGS84_TO_CDMAG = 611
ITRF_TO_EME2000 = 601
ITRF_TO_ICRF2000 = 601
ITRF_TO_GEI2000 = 601
ITRF_TO_MOD = 602
ITRF_TO_TOD = 603
ITRF_TO_TEME = 604
ITRF_TO_PEF = 605
ITRF_TO_WGS84 = 606
ITRF_TO_ITRF = 606
ITRF_TO_GEO = 606
ITRF_TO_GSE = 607
ITRF_TO_GSM = 608
ITRF_TO_SM = 609
ITRF_TO_EDMAG = 610
ITRF_TO_CDMAG = 611
GEO_TO_EME2000 = 601
GEO_TO_ICRF2000 = 601
GEO_TO_GEI2000 = 601
GEO_TO_MOD = 602
GEO_TO_TOD = 603
GEO_TO_TEME = 604
GEO_TO_PEF = 605
GEO_TO_WGS84 = 606
GEO_TO_ITRF = 606
GEO_TO_GEO = 606
GEO_TO_GSE = 607
GEO_TO_GSM = 608
GEO_TO_SM = 609
GEO_TO_EDMAG = 610
GEO_TO_CDMAG = 611
GSE_TO_EME2000 = 701
GSE_TO_ICRF2000 = 701
GSE_TO_GEI2000 = 701
GSE_TO_MOD = 702
GSE_TO_TOD = 703
GSE_TO_TEME = 704
GSE_TO_PEF = 705
GSE_TO_WGS84 = 706
GSE_TO_ITRF = 706
GSE_TO_GEO = 706
GSE_TO_GSE = 707
GSE_TO_GSM = 708
GSE_TO_SM = 709
GSE_TO_EDMAG = 710
GSE_TO_CDMAG = 711
GSM_TO_EME2000 = 801
GSM_TO_ICRF2000 = 801
GSM_TO_GEI2000 = 801
GSM_TO_MOD = 802
GSM_TO_TOD = 803
GSM_TO_TEME = 804
GSM_TO_PEF = 805
GSM_TO_WGS84 = 806
GSM_TO_ITRF = 806
GSM_TO_GEO = 806
GSM_TO_GSE = 807
GSM_TO_GSM = 808
GSM_TO_SM = 809
GSM_TO_EDMAG = 810
GSM_TO_CDMAG = 811
SM_TO_EME2000 = 901
SM_TO_ICRF2000 = 901
SM_TO_GEI2000 = 901
SM_TO_MOD = 902
SM_TO_TOD = 903
SM_TO_TEME = 904
SM_TO_PEF = 905
SM_TO_WGS84 = 906
SM_TO_ITRF = 906
SM_TO_GEO = 906
SM_TO_GSE = 907
SM_TO_GSM = 908
SM_TO_SM = 909
SM_TO_EDMAG = 910
SM_TO_CDMAG = 911
EDMAG_TO_EME2000 = 1001
EDMAG_TO_ICRF2000 = 1001
EDMAG_TO_GEI2000 = 1001
EDMAG_TO_MOD = 1002
EDMAG_TO_TOD = 1003
EDMAG_TO_TEME = 1004
EDMAG_TO_PEF = 1005
EDMAG_TO_WGS84 = 1006
EDMAG_TO_ITRF = 1006
EDMAG_TO_GEO = 1006
EDMAG_TO_GSE = 1007
EDMAG_TO_GSM = 1008
EDMAG_TO_SM = 1009
EDMAG_TO_EDMAG = 1010
EDMAG_TO_CDMAG = 1011
CDMAG_TO_EME2000 = 1101
CDMAG_TO_ICRF2000 = 1101
CDMAG_TO_GEI2000 = 1101
CDMAG_TO_MOD = 1102
CDMAG_TO_TOD = 1103
CDMAG_TO_TEME = 1104
CDMAG_TO_PEF = 1105
CDMAG_TO_WGS84 = 1106
CDMAG_TO_ITRF = 1106
CDMAG_TO_GEO = 1106
CDMAG_TO_GSE = 1107
CDMAG_TO_GSM = 1108
CDMAG_TO_SM = 1109
CDMAG_TO_EDMAG = 1110
CDMAG_TO_CDMAG = 1111

###############################################################################
###############################################################################
###############################################################################

import ctypes
from Lgm_Types import LgmInt, LgmDouble, LgmLong, LgmChar
import Lgm_Vector
import Lgm_DateAndTime

class Lgm_DateTime(ctypes.Structure):
    def __init__(self, verbose = 0):
        self.nNutationTerms = 106
        self.Verbose = verbose
        self.DUT1 = 0.0
        self.xp = 0.0
        self.yp = 0.0
        self.ddPsi = 0.0
        self.ddEps = 0.0
    @classmethod
    def assign_fields(cls):
        cls._fields_ = [ ("Date", LgmLong ), # In basic ISO format (YYYYMMDD or YYYYDDD) Represented as a single long int
            ("Year", LgmInt), #  4-digit year
            ("Month", LgmInt), # [1-12]
            ("Day", LgmInt), # Day Of Month [1-31]
            ("Doy", LgmInt), # Day Of Year [1-31]
            ("Time", LgmDouble), # Decimal value of time in hours
            ("Hour", LgmInt), # Hours [0-23]
            ("Minute", LgmInt), # Minutes [0-59]
            ("Second", LgmDouble), # Seconds [0-60] (the 60 accommodates leap seconds)
            ("Week", LgmInt), # ISO Week number [1-53]
            ("wYear", LgmInt), #  ISO Year associated with the ISO Week Number (can be different from Year)
            ("Dow", LgmInt), # ISO Day Of Week number [1-7]
            ("DowStr", LgmChar*10), # ISO Day Of Week number [1-7]
            ("fYear", LgmDouble), # Decimal year (e.g. 2004.2345)
            ("JD", LgmDouble), # Julian Date
            ("T", LgmDouble), # Julian Centuries since J2000 for this time system
            ("DaySeconds", LgmDouble), # Number of seconds in the day.
            ("TZD_sgn", LgmInt), # Sign of Time zone offset
            ("TZD_hh", LgmInt), # Time zone offset hours
            ("TZD_mm", LgmInt), # Time zone offset minutes
            ("TimeSystem", LgmInt) ] # e.g. LGM_UTC, LGM_UT1, LGM_TAI, LGM_GPS, LGM_TT, LGM_TDB, LGM_TCG, etc..


class Lgm_CTrans(ctypes.Structure):
    @classmethod
    def assign_fields(cls):
        cls._fields_ = [ ("Verbose", LgmInt),
            ("l", Lgm_DateAndTime.Lgm_DateAndTime), # Structure containing Leap Second Info
            ("UT1", Lgm_DateTime), # UT is the mean solar time at Greenwich.
                                         # UT0 is a version of UT that uses data
                                         # from many different ground stations.
                                         # UT1 is a version of UT0 in which
                                         # corrections for polar motion have been
                                         # made so that time is independant of
                                         # observing location. There is also a UT2,
                                         # but we wont use UT0 or UT2 here.
                                         # Units: Decimal hours
            ("UTC", Lgm_DateTime), # Universal Time Coordinated.
                                         # Most commonly used time system. Derived
                                         # from atomic time. It is maintained to be
                                         # within +/- 0.9s of UT1 (via addition or
                                         # subtraction(?) of leap seconds).
                                         # Units: Decimal hours
            ("DUT1", LgmDouble), # Difference between UT1 and UTC.
                                         #      DUT1 = UT1 - UTC.
                                         # This is monitored and reported as part
                                         # of the Earth Orientation Parameters
                                         # (EOP). Can be predicted a short time
                                         # into the future, but definitive values
                                         # only available retrospectively.  We set
                                         # this value to 0.0 by default. Thus in
                                         # the absence of EOP data, we assume its
                                         # initial value.
                                         # Units: Decimal seconds
            ("LOD", LgmDouble), # Length Of Day (LOD). Its the amount of extra
                                         # time in seconds that the current day has. Not
                                         # predictable. Part of EOP values.
            ("TAI", Lgm_DateTime), # International Atomic Time.
                                          #     TAI = UTC + DAT
            ("GPS", Lgm_DateTime), # GPS time
                                          #     GPS = TAI - 19s
            ("DAT", LgmDouble), # Difference between UTC and TAI.
                                          #     DAT = TAI - UTC
                                          # DAT is essentially the number of leap seconds
                                          # and are an integral number of whole seconds.
                                          # Units: Decimal seconds.
            ("TT", Lgm_DateTime), # Terestrial Time (TT).
                                          # Essentially the same thing as
                                          # "Terrestrial Dynamical Time (TDT) or
                                          # Ephmeris Time (ET). Its defined to be,
                                          #      TT = TAI + 32.184s
                                          # Units: Decimal hours.
            ("TDB", Lgm_DateTime), # Barycentric Dynamical Time.
                                          # Not used here.
                                          # Units: Decimal hours
            ("TCG", Lgm_DateTime), # Geocentric Coordinate Time.
                                          # Not used here.
                                          # Units: Decimal hours
            ("gmst", LgmDouble), # Greenwich Mean Sidereal Time
                                          # units: in radians
            ("gast", LgmDouble), # Greenwich Apparent Sidereal Time
                                          # Units: in radians
            ("xp", LgmDouble), #  Pole wander parameters.
                                          # part of EOP data.
                                          # Units: radians
            ("yp", LgmDouble ), # Pole wander parameters.
                                          # part of EOP data.
                                          # Units: radians
            ("epsilon", LgmDouble ), # Mean Obliquity of the Ecliptic
                                         # (in radians)
            ("epsilon_true", LgmDouble ), # True Obliquity of the Ecliptic
                                         #  \f$\epsilon_{true} = \epsilon + dEps\f$
                                         # (in radians)
            ("eccentricity", LgmDouble ), # Eccentricity of Earth-Sun orbit
            ("lambda_sun", LgmDouble ), #  Ecliptic Long. of Sun (in radians)
            ("earth_sun_dist", LgmDouble ), #  Earth-Sun distance (in units of earth radii)
            ("RA_sun", LgmDouble ), #  Right Ascention of Sun (in degrees)
            ("DEC_sun", LgmDouble ), # Declination of Sun (in degrees)
            ("lambda_sun_ha", LgmDouble ), # high accuracy eccliptic coords of sun
            ("r_sun_ha", LgmDouble ), # high accuracy eccliptic coords of sun
            ("beta_sun_ha", LgmDouble ), # high accuracy eccliptic coords of sun
            ("RA_sun_ha", LgmDouble ), # high accuracy Right Ascention of Sun (in degrees)
            ("DEC_sun_ha", LgmDouble ), # high accuracy Declination of Sun (in degrees)
            ("Sun", Lgm_Vector.Lgm_Vector), # direction of Sun in GEI system (unit vector)
            ("EcPole", Lgm_Vector.Lgm_Vector), # direction of Ecliptic Pole in GEI system (unit vector)
            ("psi", LgmDouble), # Geodipole tilt angle, \f$\psi\f$ (in radians)
            ("sin_psi", LgmDouble), # \f$\sin(\psi)\f$
            ("cos_psi", LgmDouble), # \f$\cos(\psi)\f$
            ("tan_psi", LgmDouble), # \f$\tan(\psi)\f$
            ("RA_moon", LgmDouble), # Right Ascention of Moon (in degrees)
            ("DEC_moon", LgmDouble), # Declination of Moon (in degrees)
            ("MoonPhase", LgmDouble), # The Phase of the Moon (in days)
            ("EarthMoonDistance", LgmDouble), # Distance between the Earth and Moon (in earth-radii)
         #  The following are various important parameters derived from
         #  the IGRF field. Note that these are the basis for defining
         #  Mag coord systems. That's why they are here and not somewhere else...
            ("M_cd", LgmDouble), # centered  dipole Magnetic moment. (nT Re^3)
            ("M_cd_McIllwain", LgmDouble), # magnetic dipole moment used by McIllwain to compute L. Sometimes want to use this for consistency?
            ("CD_gcolat", LgmDouble), #  Geographic colat of centered dipole axis (deg.)
            ("CD_glon", LgmDouble), # Geographic long. of centered dipole axis (deg.)
            ("ED_x0", LgmDouble), # x-comp of dipole displacement from center. Used in eccentric dipole field.
            ("ED_y0", LgmDouble), # y-comp of dipole displacement from center. Used in eccentric dipole field.
            ("ED_z0", LgmDouble), # z-comp of dipole displacement from center. Used in eccentric dipole field.
            ("Zeta", LgmDouble), # Precession angle, \f$\zeta\f$
            ("Theta", LgmDouble), # Precession angle, \f$\theta\f$
            ("Zee", LgmDouble), # Precession angle, \f$z\f$
            ("nNutationTerms", LgmInt), # number of terms to usek in the dPsi/dEps Nutation series.
            ("dPsi", LgmDouble), #
            ("dEps", LgmDouble), #
            ("dPsiCosEps", LgmDouble), #
            ("dPsiSinEps", LgmDouble), #
            ("ddPsi", LgmDouble), # radians additional corrections to dPsi -- part of EOP data
            ("ddEps", LgmDouble), # radians additional corrections to dEps -- part of EOP data
            ("EQ_Eq", LgmDouble), # Equation of the equinoxes.
            ("OmegaMoon", LgmDouble), # Ascending node of Moon.
            ("dX", LgmDouble), #  for IUA-2000A reduction (not used yet)
            ("dY", LgmDouble), # for IUA-2000A reduction (not used yet)
            # Transformation matrices between various ccord systems
            ("Agei_to_mod", LgmDouble * 3 * 3), #
            ("Amod_to_gei", LgmDouble * 3 * 3), #
            ("Amod_to_tod", LgmDouble * 3 * 3), #
            ("Atod_to_mod", LgmDouble * 3 * 3), #
            ("Ateme_to_pef", LgmDouble * 3 * 3), #
            ("Apef_to_teme", LgmDouble * 3 * 3), #
            ("Apef_to_tod", LgmDouble * 3 * 3), #
            ("Atod_to_pef", LgmDouble * 3 * 3), #
            ("Awgs84_to_pef", LgmDouble * 3 * 3), #
            ("Apef_to_wgs84", LgmDouble * 3 * 3), #
            ("Agse_to_mod", LgmDouble * 3 * 3), #
            ("Amod_to_gse", LgmDouble * 3 * 3), #
            ("Asm_to_gsm", LgmDouble * 3 * 3), #
            ("Agsm_to_sm", LgmDouble * 3 * 3), #
            ("Agsm_to_mod", LgmDouble * 3 * 3), #
            ("Amod_to_gsm", LgmDouble * 3 * 3), #
            ("Agsm_to_gse", LgmDouble * 3 * 3), #
            ("Agse_to_gsm", LgmDouble * 3 * 3), #
            ("Awgs84_to_mod", LgmDouble * 3 * 3), #
            ("Amod_to_wgs84", LgmDouble * 3 * 3), #
            ("Awgs84_to_gei", LgmDouble * 3 * 3), #
            ("Agei_to_wgs84", LgmDouble * 3 * 3), #
            ("Agsm_to_wgs84", LgmDouble * 3 * 3), #
            ("Awgs84_to_gsm", LgmDouble * 3 * 3), #
            ("Awgs84_to_cdmag", LgmDouble * 3 * 3), #
            ("Acdmag_to_wgs84", LgmDouble * 3 * 3), #
            ("Agei_to_mod", LgmDouble * 3 * 3), #
            ("Amod_to_gei", LgmDouble * 3 * 3), #
            ("Amod_to_tod", LgmDouble * 3 * 3), #
            ("Atod_to_mod", LgmDouble * 3 * 3), #
            ("Ateme_to_pef", LgmDouble * 3 * 3), #
            ("Apef_to_teme", LgmDouble * 3 * 3), #
            ("Apef_to_tod", LgmDouble * 3 * 3), #
            ("Atod_to_pef", LgmDouble * 3 * 3), #
            ("Awgs84_to_pef", LgmDouble * 3 * 3), #
            ("Apef_to_wgs84", LgmDouble * 3 * 3), #
            ("Agse_to_mod", LgmDouble * 3 * 3), #
            ("Amod_to_gse", LgmDouble * 3 * 3), #
            ("Asm_to_gsm", LgmDouble * 3 * 3), #
            ("Agsm_to_sm", LgmDouble * 3 * 3), #
            ("Agsm_to_mod", LgmDouble * 3 * 3), #
            ("Amod_to_gsm", LgmDouble * 3 * 3), #
            ("Agsm_to_gse", LgmDouble * 3 * 3), #
            ("Agse_to_gsm", LgmDouble * 3 * 3), #
            ("Awgs84_to_mod", LgmDouble * 3 * 3), #
            ("Amod_to_wgs84", LgmDouble * 3 * 3), #
            ("Awgs84_to_gei", LgmDouble * 3 * 3), #
            ("Agei_to_wgs84", LgmDouble * 3 * 3), #
            ("Agsm_to_wgs84", LgmDouble * 3 * 3), #
            ("Awgs84_to_gsm", LgmDouble * 3 * 3), #
            ("Awgs84_to_cdmag", LgmDouble * 3 * 3), #
            ("Acdmag_to_wgs84", LgmDouble * 3 * 3), #
            # These variables are needed to make IGRF Calls reentrant/thread-safe.
            ("Lgm_IGRF_FirstCall", LgmInt), #
            ("Lgm_IGRF_OldYear", LgmDouble), #
            ("Lgm_IGRF_g", LgmDouble * 13 * 13), #
            ("Lgm_IGRF_h", LgmDouble * 13 * 13), #
            ("Lgm_IGRF_R", LgmDouble * 13 * 13), #
            ("Lgm_IGRF_K", LgmDouble * 13 * 13), #
            ("Lgm_IGRF_S", LgmDouble * 13 * 13), #
            ("Lgm_IGRF_TwoNm1_Over_NmM", LgmDouble * 13 * 13), #
            ("Lgm_IGRF_NpMm1_Over_NmM", LgmDouble * 13 * 13), #
            ("Lgm_IGRF_SqrtNM1", LgmDouble * 13 * 13), #
            ("Lgm_IGRF_SqrtNM2", LgmDouble * 13 * 13) ]  #
        # Mike has Lgm_init_ctrans that sets a few vars, set them here







#class Lgm_DateTime(ctypes.Structure):
#    _fields_ = [ ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ), #
#        ("", ) ] #
