"""Cone query on DELVE-MC DR3 + morphocloud star/galaxy separation.

query(ra, dec, radius, stars=True, galaxies=False, rmax=None)
    -> astropy Table of stars, galaxies, or (stars, galaxies)
    rmax=None uses catalog.mag_max from config.yaml in this package directory.

Author: Mahdieh Navabi
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyvo
import yaml
from astropy.table import Table

from morphocloud.features import engineer_features
from morphocloud.infer import StarGalaxyClassifier

R_g, R_r, R_i, R_z = 3.1863, 2.1401, 1.5690, 1.196   # DECam A_band / E(B-V)
TAP_URL = "https://datalab.noirlab.edu/tap"
MAG_SENT, ERR_SENT, FWHM_SENT = 99.0, 9.9, 999990.0  # DELVE-MC missing-value conventions
CONFIG_FILE = Path(__file__).resolve().parent / "config.yaml"

_clf = None


def _model():
    global _clf
    if _clf is None:
        _clf = StarGalaxyClassifier.load()
    return _clf


def default_rmax():
    """r-band magnitude limit from catalog.mag_max in config.yaml."""
    with open(CONFIG_FILE, "r") as f:
        cfg = yaml.safe_load(f)
    return float(cfg["catalog"]["mag_max"])


def fetch(ra, dec, radius, rmax=None):
    if rmax is None:
        rmax = default_rmax()
    # RA box is a cheap pre-filter; the ACOS distance cut below stays for exactness
    dra = radius / np.cos(np.radians(dec))
    adql = f"""
    SELECT objid, brick, brickuniq,
           ra, dec,
           gmag, gerr, grms, ndetg,
           rmag, rerr, rrms, ndetr,
           imag, ierr, irms, ndeti,
           zmag, zerr, zrms, ndetz,
           ebv, ndet, nalfdetiter,
           chi, sharp, prob, mag_auto, fwhm, ellipticity, asemi, bsemi
    FROM delve_dr3.mc_object
    WHERE dec BETWEEN {dec - radius} AND {dec + radius}
      AND ra BETWEEN {ra - dra} AND {ra + dra}
      AND DEGREES(ACOS(SIN(RADIANS(dec))*SIN(RADIANS({dec})) + COS(RADIANS(dec))*COS(RADIANS({dec}))*COS(RADIANS(ra - {ra})))) <= {radius}
      AND rmag < {rmax}
    """
    svc = pyvo.dal.TAPService(TAP_URL)
    # The Data Lab server sometimes stalls; retry a few times before giving up
    for attempt in range(1, 4):
        try:
            df = svc.run_async(adql, maxrec=10_000_000).to_table().to_pandas()
            break
        except Exception as exc:  # timeouts show up as requests/pyvo errors
            if attempt == 3:
                raise
            print(f"TAP query failed ({type(exc).__name__}), retrying ({attempt}/3)...")
            time.sleep(10 * attempt)
    # Computed SQL columns (AS ...) come back from the TAP service as 1-character
    # strings, so derive the dereddened magnitudes here instead.
    df["ebv"] = df["ebv"].astype(float)
    for lo, R in (("g", R_g), ("r", R_r), ("i", R_i), ("z", R_z)):
        mag = df[f"{lo}mag"].astype(float)
        df[f"{lo}mag0"] = (mag - R * df["ebv"]).where(mag < MAG_SENT)
    return df


def _clean(x, sent):
    v = np.asarray(x, dtype=np.float64).copy()
    v[v >= sent] = np.nan
    return v


def to_raw(dl):
    """Map Data Lab columns to the RAW_INPUT_COLUMNS used by morphocloud."""
    raw = pd.DataFrame({"OBJID": dl["objid"].astype(str), "RA": dl["ra"], "DEC": dl["dec"]})
    for b, lo in [("G", "g"), ("R", "r"), ("I", "i"), ("Z", "z")]:
        raw[f"{b}MAG"] = _clean(dl[f"{lo}mag"], MAG_SENT)
        raw[f"{b}ERR"] = _clean(dl[f"{lo}err"], ERR_SENT)
        raw[f"{b}SCATTER"] = _clean(dl[f"{lo}rms"], MAG_SENT)
        raw[f"NDET{b}"] = dl[f"ndet{lo}"].astype(np.int32)
    for col, src in [("CHI", "chi"), ("SHARP", "sharp"), ("PROB", "prob"),
                     ("EBV", "ebv"), ("MAG_AUTO", "mag_auto")]:
        raw[col] = _clean(dl[src], MAG_SENT)
    for col, src in [("FWHM", "fwhm"), ("ELLIPTICITY", "ellipticity"),
                     ("ASEMI", "asemi"), ("BSEMI", "bsemi")]:
        raw[col] = _clean(dl[src], FWHM_SENT)
    raw["BRICKUNIQ"] = dl["brickuniq"].astype(bool).values
    return raw


def classify(dl):
    """Add P_STAR and STAR_LEAK1 (leak1 operating point) to the fetched table."""
    raw = to_raw(dl)
    # seeing: median FWHM (arcsec) of point-like detections in this cone
    ok = raw["FWHM"].notna() & (raw["PROB"] >= 0.8)
    seeing = float(raw.loc[ok, "FWHM"].median())

    feats = engineer_features(raw, seeing=seeing)
    clf = _model()
    p = clf.predict(feats)
    thr = clf.smooth_threshold(operating_point="leak1")
    star = p >= thr(feats["RMAG0"].to_numpy(), strictness=0.0, flat=0.5)

    out = dl.copy()
    out["P_STAR"] = p
    out["STAR_LEAK1"] = star
    return out


def query(ra, dec, radius, stars=True, galaxies=False, rmax=None):
    if rmax is None:
        rmax = default_rmax()
    dl = fetch(ra, dec, radius, rmax=rmax)
    dl = dl[dl["brickuniq"].map({"t": True, "f": False, True: True, False: False}).astype(bool)]
    if len(dl) == 0:
        empty = Table()
        if stars and galaxies:
            return empty, empty
        return empty

    labelled = classify(dl)
    star_tab = Table.from_pandas(labelled[labelled["STAR_LEAK1"]].reset_index(drop=True))
    gal_tab = Table.from_pandas(labelled[~labelled["STAR_LEAK1"]].reset_index(drop=True))

    if stars and galaxies:
        return star_tab, gal_tab
    if stars:
        return star_tab
    if galaxies:
        return gal_tab
    return None
