import sys

from scipy.stats import binom

from keys import Keys


def _h(cov: int, vr: int, err: float, k: int, l: int):
    """
    helper function
    """
    return binom.pmf(l, k, err) * binom.pmf(vr - k + l, cov - k, err / 3)


def fpProb(cov: int, vr: int, err: float, mvr: int = 1):
    fp = 0
    for k in range(0, mvr):
        for l in range(0, k + 1):
            fp += _h(cov, vr, err, k, l)
    return fp


def getTPFP(cov: int, vr: int, err: float, mvr: int = 1):
    if cov < 0:
        raise ValueError("Coverage must be non-negative (COV >= 0)")
    if vr < 0:
        raise ValueError("Variant reads must be non-negative (VR >= 0)")
    if vr > cov:
        raise ValueError("Variant reads must be less than coverage (VR <= COV)")
    vaf = vr / cov
    if err < 0:
        raise ValueError("Sequencing error must be non-negative (ERR >= 0)")
    if err > vaf:
        raise ValueError("Sequencing error must be less than VAF (ERR <= VAF)")
    if mvr < 0:
        raise ValueError("Minimum variant reads must be non-negative (MVR >= 0)")
    if mvr > vr:
        raise ValueError(
            "Minimum variant reads must be less than variant reads (MVR <= VR)"
        )
    # NOTE: useless
    if not 0 < mvr <= vr <= cov:
        raise ValueError("Condition not satisfied: 0 < MVR <= VR <= COV")

    fp = fpProb(cov, vr, err, mvr)
    fp = min(fp, 1)
    tp = 1 - fp

    return {
        Keys.FP: fp,
        Keys.TP: tp,
    }
