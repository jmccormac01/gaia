"""
Script to query the GAIA catalogue for parallaxes

Modified version of queryGaia to use output paramfit
file from Barry Smalley
"""
import os
import re
import argparse as ap
from collections import (
    OrderedDict,
    defaultdict
    )
import pymysql
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astroquery.vizier import Vizier

# pylint: disable = invalid-name
# pylint: disable = redefined-outer-name
# pylint: disable = no-member
# pylint: disable = undefined-loop-variable

GAIA_CATALOGUE_ID = 'I/337/tgasptyc'

def argParse():
    """
    Parse the command line arguments

    Parameters
    ----------
    None

    Returns
    -------
    args : array-like
        Parsed ArgumentParser object containing the
        command line arguments

    Raises
    ------
    None
    """
    parser = ap.ArgumentParser()
    parser.add_argument('--radius',
                        help='search radius in arcsec',
                        type=int,
                        default=5)
    return parser.parse_args()

def swaspIdToSkyCoord(swasp_id):
    """
    Convert a swasp_id to SkyCoord object

    Parameters
    ----------
    swasp_id : str
        ID of the swasp object

    Returns
    -------
    coords : astropy.coordinates.SkyCoord
        SkyCoord object for the given swasp_id

    Raises
    ------
    None
    """
    # set up a regex string for the swasp_id format
    p = re.compile(r'1SWASPJ(?P<ra1>\d\d)(?P<ra2>\d\d)(?P<ra3>\d\d.\d\d)'
                   r'(?P<dec1>.\d\d)(?P<dec2>\d\d)(?P<dec3>\d\d.\d)')
    match = re.findall(p, swasp_id)[0]
    if len(match) == 6:
        ra = ":".join(match[:3])
        dec = ":".join(match[3:])
        coords = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.degree), frame='icrs')
    else:
        coords = None
    return coords

def queryGaiaAroundSwaspId(swasp_id, radius):
    """
    Query GAIA around the coordinates of a given
    SuperWASP object

    Parameters
    ----------
    swasp_id : str
        SuperWASP object ID, e.g. 1SWASPJ102030.65+222839.0
    radius : int
        Search radius in arcsec

    Returns
    -------
    objects : array-like
        Dictionary of matching objects. Returns empty OrderedDict
        if nothing is found

    Raises
    ------
    None
    """
    coordinates = swaspIdToSkyCoord(swasp_id)
    objects = OrderedDict()
    if coordinates:
        vizier = Vizier(columns=['TYC', 'HIP', '_RAJ2000', '_DEJ2000', 'Plx', 'e_Plx',
                                 'pmRA', 'pmDE', 'Source'],
                        keywords=['optical'])
        vizier.ROW_LIMIT = -1
        results = vizier.query_region(coordinates,
                                      radius=Angle(radius*u.arcsec),
                                      catalog=GAIA_CATALOGUE_ID)
        try:
            for result in results[0]:
                objects[result['Source']] = {'tyc': result['TYC'],
                                             'hip': result['HIP'],
                                             'ra': round(float(result['_RAJ2000']), 8),
                                             'dec': round(float(result['_DEJ2000']), 8),
                                             'plx': round(float(result['Plx']), 4),
                                             'eplx': round(float(result['e_Plx']), 4),
                                             'pmra': round(float(result['pmRA']), 4),
                                             'pmdec': round(float(result['pmDE']), 4)}
        except (TypeError, IndexError):
            print('No match in GAIA for {}...'.format(swasp_id))
    return objects

def rStar(theta, parallax):
    """
    Calculate r_star from IRFM theta and parallax

    Parameters
    ----------
    theta : float
        IRFM theta value from SWASP paramfit (mas, from Barry Smalley)
    parallax : float
        GAIA parallax (mas)

    Returns
    -------
    r_star : float
        Stellar radius in R_sun

    Raises
    ------
    None
    """
    r_star = 214.9*(theta/2.)/parallax
    return r_star

def getParamfitResults():
    """
    Grab the swasp_id and theta values from the database

    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    None
    """
    swasp_ids, thetas, chi2s, residuals = [], [], [], []
    qry = """
        SELECT swasp_id, theta, chi2, residuals
        FROM paramfit
         """
    with pymysql.connect(host='localhost', db='eblm', password='mysqlpassword') as cur:
        cur.execute(qry)
        for row in cur:
            swasp_ids.append(row[0])
            thetas.append(row[1])
            chi2s.append(row[2])
            residuals.append(row[3])
    return swasp_ids, thetas, chi2s, residuals

if __name__ == "__main__":
    # parse the command line args
    args = argParse()
    # open database connection
    db = pymysql.connect(host='localhost', db='eblm', password='mysqlpassword')
    # hold the output data
    matches = OrderedDict()
    n_gaia_matches = {}
    swasp_ids, thetas, chi2s, residuals = getParamfitResults()

    # loop over the objects and get the matches
    # and r_stars for any matches
    for swasp_id, theta, chi2, resid in zip(swasp_ids, thetas, chi2s, residuals):
        print('\nQuerying GAIA for {}:'.format(swasp_id))
        print('ID                   TYC         RA            DEC           '
              'PMRA      PMDEC     PLX       ePLX      THETA     CHI2      RESID     RSTAR')

        matches[swasp_id] = queryGaiaAroundSwaspId(swasp_id, args.radius)
        n_gaia_matches[swasp_id] = len(matches[swasp_id])

        # put a warning when theta might be bad
        if chi2 > 3 or resid > 0.2:
            token = '?'
        else:
            token = ''

        for match in matches[swasp_id]:
            parallax = matches[swasp_id][match]['plx']
            r_star = rStar(theta, parallax)
            if r_star >= 2.0:
                giant_flag = 'giant'
            else:
                giant_flag = 'dwarf'
            # add the ? token if chi2 or residual too large to be sure
            giant_flag = giant_flag + token
            # check for multiple matches, if more than 1 just put the number of matches in the db table
            # those with >1 will need looked at separately
            if n_gaia_matches[swasp_id] == 1:
                qry = """
                    UPDATE eblm_parameters
                    SET n_gaia_matches = {},
                    r_star = {},
                    giant_flag = '{}'
                    WHERE swasp_id = '{}'
                    """.format(n_gaia_matches[swasp_id],
                               round(float(r_star), 3),
                               giant_flag,
                               swasp_id)
            else:
                qry = """
                    UPDATE eblm_parameters
                    SET n_gaia_matches = {}
                    WHERE swasp_id = '{}'
                    """.format(n_gaia_matches[swasp_id],
                               swasp_id)
            with db.cursor() as cur:
                cur.execute(qry)
                db.commit()

            matches[swasp_id][match]['rstar'] = round(float(r_star), 4)
            matches[swasp_id][match]['theta'] = round(float(theta), 4)
            matches[swasp_id][match]['chi2'] = round(float(chi2), 4)
            matches[swasp_id][match]['resid'] = round(float(resid), 4)
            print("{:<20} {:<11} {:<13} {:<13} {:<9} {:<9} {:<9} {:<9} {:<9} {:<9} {:<9} {:<9}{}".format(match,
                                matches[swasp_id][match]['tyc'].decode('ascii'),
                                matches[swasp_id][match]['ra'],
                                matches[swasp_id][match]['dec'],
                                matches[swasp_id][match]['pmra'],
                                matches[swasp_id][match]['pmdec'],
                                matches[swasp_id][match]['plx'],
                                matches[swasp_id][match]['eplx'],
                                matches[swasp_id][match]['theta'],
                                matches[swasp_id][match]['chi2'],
                                matches[swasp_id][match]['resid'],
                                matches[swasp_id][match]['rstar'],
                                token))
            print(qry)
    db.close()

