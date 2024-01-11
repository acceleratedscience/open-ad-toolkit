# Example command:
# display collection details 'Patents from USPTO'

from openad.helpers.output import output_text, output_error
from openad.helpers.output_msgs import msg


def display_collection_details(inputs: dict, cmd_pointer):
    """
    Displays the details for a given collection.

    Parameters
    ----------
    inputs:
        Parser inputs from pyparsing.
    cmd_pointer:
        Pointer to runtime.
    """

    api = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("DS4SD")]
    try:
        collections = api.elastic.list()
        # raise Exception('This is a test error')
    except Exception as err:  # pylint: disable=broad-exception-caught
        output_error(msg("err_deepsearch", err), return_val=False)
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
        output_error("No collection found", return_val=False)
        return False

    output = [
        "<bold>" + collection.name + "</bold> ",
        "<green>Collection Name: </green>" + collection.name,
        "<green>Domains: </green>" + str(" / ".join(collection.metadata.domain)),
        "<green>Type: </green>" + collection.metadata.type,
        "<green>Collection Key: </green>" + collection.source.index_key,
        "<green>Documents: </green>" + str(collection.documents),
        "<green>Created: </green>" + collection.metadata.created.strftime("%Y-%m-%d"),
        "<green>URL: </green><link>" + collection.metadata.source + "</link>",
        "<green>Description: </green>" + collection.metadata.description,
    ]

    return output_text("\n".join(output), pad=1)
