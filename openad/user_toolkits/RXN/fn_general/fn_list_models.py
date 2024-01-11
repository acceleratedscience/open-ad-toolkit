# Example command:
# list rxn models

from openad.helpers.output import output_table
from openad.app.global_var_lib import GLOBAL_SETTINGS


def list_models(inputs: dict, cmd_pointer):
    """
    List available RXN models
    """

    rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][cmd_pointer.login_settings["toolkits"].index("RXN")]
    # Prepare the data query
    try:
        x = rxn4chemistry_wrapper.list_models()
    except Exception as e:  # pylint: disable=broad-exception-caught
        raise Exception("unable to load models :" + str(e)) from e  # pylint: disable=broad-exception-raised
    results = []
    results2 = []
    for i in x:
        outstr = []
        for ii in x[i]:
            outstr.append(ii["name"])

        results.append(i)
        results2.append(outstr)
    res_dict = {"Models": results, "versions": results2}
    import pandas as pd

    df = pd.DataFrame.from_dict(res_dict)
    df.style.hide(axis="index")
    if GLOBAL_SETTINGS["display"] == "notebook":
        from IPython.display import HTML

        return HTML(df.to_html(index=False))

    output_table(df, headers=["Models", "Versions"])
