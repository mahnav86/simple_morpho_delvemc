"""Drop-in replacement for query_dl.query: DELVE-MC DR3 cone search with
morphocloud star/galaxy separation (no SQL star/galaxy cuts).

Author: Mahdieh Navabi

Same signature as query_dl.query. `profile` is accepted for compatibility
but not used (the public Data Lab TAP service needs no profile).
"""
from functools import lru_cache

import numpy as np
from astropy.table import Table

from simple_adl_morpho.morpho_separation import fetch, classify


@lru_cache(maxsize=8)
def _labelled(ra, dec, radius, gmax):
    """Fetch + classify once per cone; the stars and galaxies calls share it."""
    dl = fetch(ra, dec, radius, rmax=gmax)
    dl = dl[dl["brickuniq"].map({"t": True, "f": False, True: True, False: False}).astype(bool)]
    if len(dl) == 0:
        return None
    return classify(dl)


def query(profile, ra, dec, radius, gmax, stars=True, galaxies=False):
    labelled = _labelled(float(ra), float(dec), float(radius), float(gmax))
    if labelled is None:
        empty = np.array([])
        return (empty, empty) if (stars and galaxies) else empty

    star_tab = labelled[labelled["STAR_LEAK1"]].reset_index(drop=True).to_records(index=False)
    gal_tab = labelled[~labelled["STAR_LEAK1"]].reset_index(drop=True).to_records(index=False)

    # return numpy structured arrays, as query_dl did (plot code relies on this)
    if stars and galaxies:
        return star_tab, gal_tab
    if stars:
        return star_tab
    if galaxies:
        return gal_tab
    return None
