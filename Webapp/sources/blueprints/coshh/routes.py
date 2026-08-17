import time

from . import coshh_bp


@coshh_bp.route("/coshh_setup/time")
def coshh_setup():
    print(time.time())
    return {"time": time.time()}
