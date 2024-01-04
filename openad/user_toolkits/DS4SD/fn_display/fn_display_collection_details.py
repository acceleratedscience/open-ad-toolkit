from openad.helpers.output import output_text, output_error
from openad.app.global_var_lib import GLOBAL_SETTINGS

_tableformat = "simple"


def display_collection_details(inputs: dict, cmd_pointer):
    """Displays the details for a given collection
    inputs: parser inputs from pyparsing
       cmd_pointer: pointer to runtime"""
    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    try:
        collections = api.elastic.list()
    except Exception as e:  # pylint: disable=broad-exception-caught
        output_error("Error in calling deepsearch:" + str(e), return_val=False)
        return False

    collection = None
    for c in collections:
        if inputs["collection"] == c.name:
            collection = c
            break
        if inputs["collection"] == c.source.index_key:
            collection = c
            break

    if collection == None:
        output_error(" No Collection Found", return_val=False)
        return False

    output = [
        " ",
        "<h1>" + collection.name + " </h1> ",
        "<success>Collection Name:</success> " + collection.name,
        " <success>Domains: </success>" + str(" / ".join(collection.metadata.domain)),
    ]
    output.append("<success>Type: </success>" + collection.metadata.type)
    output.append("<success>Collection Key: </success>" + collection.source.index_key)
    output.append("<success>Documents: </success>" + str(collection.documents))
    output.append("<success>Created: </success>" + collection.metadata.created.strftime("%Y-%m-%d"))
    output.append("<success>URL: </success>" + collection.metadata.source)
    output.append("<success>Description:</success> " + collection.metadata.description)

    if GLOBAL_SETTINGS["display"] == "notebook":  # Have to do individual rows because of Markdown Behaviour
        for x in output:
            output_text(x, return_val=False)
        return True
    else:
        return output_text("\n".join(output) + "\n", return_val=True)
