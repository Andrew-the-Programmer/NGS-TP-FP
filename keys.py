class Keys:
    VAF = "vaf"
    COV = "cov"
    VR = "variant reads"
    ERR = "error"
    FP = "false positive"
    TP = "true positive"
    MVR = "minimum variant reads"
    ERR_MSG = "error message"


def _format_percent(p, decimals=1):
    p = p * 100
    if p < 1:
        return f"{p:,.{decimals}g}"
    if p > 99:
        return 100 - float(_format_percent(1 - p / 100, decimals))
    return f"{p:,.{decimals}f}"


def format_percent(p, decimals=1):
    return f"{_format_percent(p, decimals)}%"
