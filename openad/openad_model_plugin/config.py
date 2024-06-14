import os


##########################################################################
# region - Gloabls
##########################################################################

DISPATCHER_SERVICE_PATH = os.path.expanduser("~/.servicing/")

SERVICES_PATH = "/definitions/services/"

SERVICE_MODEL_PATH = os.path.expanduser("~/.openad_model_services/")
if not os.path.exists(SERVICE_MODEL_PATH):
    # create dir on runtime
    os.makedirs(SERVICE_MODEL_PATH)
