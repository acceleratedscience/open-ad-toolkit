# Example command:
# list rxn models

import pandas as pd
from openad.helpers.output import output_table, output_error


def list_models(inputs: dict, cmd_pointer):
    """
    List available RXN models
    """

    rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][cmd_pointer.login_settings["toolkits"].index("RXN")]

    # Load models
    try:
        x = rxn4chemistry_wrapper.list_models()
        # raise Exception('This is a test error')
    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(["Unable to load models", err], return_val=False)
        return

    # Print result
    results = []
    results2 = []
    for i in x:
        outstr = []
        for ii in x[i]:
            outstr.append(ii["name"])

        results.append(i)
        results2.append(outstr)
    res_dict = {"Models": results, "versions": results2}
    df = pd.DataFrame.from_dict(res_dict)

    return output_table(df, headers=["Models", "Versions"])
