from . import coshh_api_bp


@coshh_api_bp.route("/test_response", methods=["GET", "POST"])
def test_response():
    return {
        "message": "Connection Successful!",
        "user_status": "Authenticated",
        # "workgroup": workgroup
    }
