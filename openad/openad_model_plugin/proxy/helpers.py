"""
Proxy helper functions for interacting with guardian
"""

import jwt, time
from typing import Dict, Any

from openad.openad_model_plugin.utils import get_logger


logger = get_logger(__name__)


def jwt_decode(bearer_token) -> Dict[str, Any]:
    """Returns a dictionary object of the jwt information

    Args:
        bearer_token (str): bearer token for proxy header

    Returns:
        Dict[str, Any]: JWT decoded information including expiry, and available services, etc.
    """
    try:
        decoded_token: dict = jwt.decode(
            bearer_token, options={"verify_at_hash": False, "verify_signature": False}, verify=False
        )
        # Convert expiry time to a human-readable format and update dict
        decoded_token["exp_formatted"] = time.strftime("%a %b %e, %G  at %R", time.localtime(decoded_token["exp"]))
    except Exception as e:
        logger.debug("Could not decode time from proxy jwt")
        logger.debug(str(e))
        decoded_token = {}
    return decoded_token
