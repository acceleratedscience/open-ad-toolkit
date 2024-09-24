import pyparsing as py
from openad.openad_model_plugin.utils import bcolors, get_logger
from openad.core.help import help_dict_create

# fmt: off
# Repeated description sections
MODEL_SERVICE_QUOTES = "<soft>Note: '<service_name>' can optionally be in quotes to support spaces.</soft>"
MODEL_SERVICE_AUTH_GROUP_QUOTES = "<soft>Note: '<service_name>' and '<auth_group>' can optionally be in quotes to support spaces.</soft>"
AUTH_GROUP_QUOTES = "<soft>Note: '<auth_group>' can optionally be in quotes to support spaces.</soft>"
# fmt: on

# Logger
logger = get_logger(__name__, color=bcolors.OKCYAN + bcolors.UNDERLINE)


def model_service_grammar_add(statements: list, help: list):
    """Create the required grammar for managing cataloging services and model up or down."""

    logger.debug("catalog model service grammer")
    catalog = py.CaselessKeyword("catalog")
    uncatalog = py.CaselessKeyword("uncatalog")
    model = py.CaselessKeyword("model")
    up = py.CaselessKeyword("up")
    local = py.CaselessKeyword("local")
    down = py.CaselessKeyword("down")
    service = py.CaselessKeyword("service")
    status = py.CaselessKeyword("status")
    fr_om = py.CaselessKeyword("from")
    _list = py.CaselessKeyword("list")
    quoted_string = py.QuotedString("'", escQuote="\\")

    auth_group = quoted_string | py.Word(py.alphanums + "_")
    service_name = quoted_string | py.Word(py.alphanums + "_")

    a_s = py.CaselessKeyword("as")
    describe = py.CaselessKeyword("describe")
    remote = py.CaselessKeyword("remote")
    auth = py.CaselessKeyword("auth")
    group = py.CaselessKeyword("group")
    _with = py.CaselessKeyword("with")
    add = py.CaselessKeyword("add")
    remove = py.CaselessKeyword("remove")
    to = py.CaselessKeyword("to")

    # catalog service
    using_keyword = py.CaselessKeyword("USING").suppress()
    quoted_identifier = py.QuotedString("'", escChar="\\", unquoteResults=True)
    parameter = py.Word(py.alphas, py.alphanums + "-_") | quoted_identifier
    value = py.Word(py.alphanums + "-_") | quoted_identifier
    param_value_pair = py.Group(parameter + py.Suppress("=") + value)
    using_clause = py.Optional(
        using_keyword + py.Suppress("(") + py.Optional(py.OneOrMore(param_value_pair))("params") + py.Suppress(")")
    )

    statements.append(py.Forward(model + auth + _list)("list_auth_services"))
    help.append(
        help_dict_create(
            name="model auth list",
            category="Model Service",
            command="model auth list",
            description="Show authentication group mapping",
        )
    )

    # fmt: off
    statements.append(py.Forward(model + auth + add + group + auth_group("auth_group") + _with + quoted_string("api_key"))("add_service_auth_group"))
    # fmt: on
    help.append(
        help_dict_create(
            name="model auth add group",
            category="Model Service",
            command="model auth add group <auth_group> with '<api_key>'",
            description=f"""
Create an authentication group for model services to use.

{AUTH_GROUP_QUOTES}
""",
        )
    )

    statements.append(py.Forward(model + auth + remove + group + auth_group("auth_group"))("remove_service_auth_group"))
    help.append(
        help_dict_create(
            name="model auth remove group",
            category="Model Service",
            command="model auth remove group <auth_group>",
            description=f"""
Remove an authentication group.

{AUTH_GROUP_QUOTES}
""",
        )
    )

    # fmt: off
    statements.append(py.Forward(model + auth + add + service + service_name("service_name") + to + group + auth_group("auth_group"))("attach_service_auth_group"))
    # fmt: on
    help.append(
        help_dict_create(
            name="model auth add service",
            category="Model Service",
            command="model auth add service <service_name> to group <auth_group>",
            description=f"""
Attach an authentication group to a model service.

{MODEL_SERVICE_AUTH_GROUP_QUOTES}
""",
        )
    )

    # fmt: off
    statements.append(py.Forward(model + auth + remove + service + service_name("service_name"))("detach_service_auth_group"))
    # fmt: on
    help.append(
        help_dict_create(
            name="model auth remove service",
            category="Model Service",
            command="model auth remove service <service_name>",
            description=f"""
Detach an authentication group from a model service.

{MODEL_SERVICE_QUOTES}
""",
        )
    )

    statements.append(py.Forward(model + service + status)("model_service_status"))
    help.append(
        help_dict_create(
            name="model service status",
            category="Model Service",
            command="model service status",
            description="Get the status of all currently cataloged services.",
        )
    )

    statements.append(py.Forward(model + service + describe + (service_name)("service_name"))("model_service_config"))
    help.append(
        help_dict_create(
            name="model service describe",
            category="Model Service",
            command="model service describe <service_name>",
            description=f"""
Get the configuration of a service.

{MODEL_SERVICE_QUOTES}
""",
        )
    )

    statements.append(py.Forward(model + catalog + _list)("get_catalog_namespaces"))
    help.append(
        help_dict_create(
            name="model catalog list",
            category="Model Service",
            command="model catalog list",
            description="Get the list of currently cataloged services.",
        )
    )

    statements.append(py.Forward(model + uncatalog + service + service_name("service_name"))("uncatalog_model_service"))
    statements.append(
        py.Forward(uncatalog + model + service + service_name("service_name"))("uncatalog_model_service")
    )  # Backward compatibility for inconsistent verb+noun
    help.append(
        help_dict_create(
            name="model uncatalog service",
            category="Model Service",
            command="model uncatalog service <service_name>",
            description=f"""
Uncatalog a model service.

{MODEL_SERVICE_QUOTES}

Example:
<cmd>uncatalog model service 'gen'</cmd>
""",
        )
    )

    # fmt: off
    statements.append(py.Forward(model + catalog + service + fr_om + py.Optional(remote("remote")) + quoted_string("path") + a_s + (quoted_string | py.Word(py.alphanums + "_"))("service_name") + using_clause)("catalog_add_model_service"))
    statements.append(py.Forward(catalog + model + service + fr_om + py.Optional(remote("remote")) + quoted_string("path") + a_s + (quoted_string | py.Word(py.alphanums + "_"))("service_name") + using_clause)("catalog_add_model_service")) # Backward compatibility for inconsistent verb+noun
    # fmt: on
    help.append(
        help_dict_create(
            name="model catalog service",
            category="Model Service",
            command="model catalog service from '<path>' | '<github>' | remote '<service_url>' as <service_name> using (<parameter>=<value> <parameter>=<value>)",
            description=f"""
Catalog a model service from a directory path or remotely from a GitHub repository or an existing OpenAD service url.

Use the <cmd>remote</cmd> clause when cataloging a service from an OpenAD service url (ip address and port).

The <cmd>service_name</cmd> is how you will refer to the service once it's up, e.g. <cmd>prop</cmd> for the property inference service.

The <cmd>using</cmd> clause can contain optional header parameters for communication with the service backend. When using a remote service, the following parameters are required:
-<cmd>Inference-Service</cmd>: The name of the inference service that is hosted.
-<cmd>auth_group</cmd> OR <cmd>Authorization</cmd>: The name of an authorization group which contains the <cmd>api_key</cmd> linked to the service access OR the <cmd>bearer_token</cmd>/<cmd>api_key</cmd> provided directly.

{MODEL_SERVICE_QUOTES}

Examples:
- <cmd>catalog model service from '/services/property_service' as prop</cmd>
- <cmd>catalog model service from 'git@github.com:acceleratedscience/generation_inference_service.git' as 'gen'</cmd>
- <cmd>catalog model service from remote 'http://54.235.3.243:30001' as molf using (Inference-Service=molformer auth_group=group_1)</cmd>
- <cmd>catalog model service from remote 'http://54.235.3.243:30001' as 'gen' using (Inference-Service=generation Authorization='KW6BSP2M9...')</cmd>
""",
        )
    )
    # fmt: off
    statements.append(py.Forward(model + service + up + service_name("service_name") + py.Optional(py.CaselessKeyword("NO_GPU")("no_gpu")))("service_up"))
    # fmt: on
    help.append(
        help_dict_create(
            name="model up",
            category="Model Service",
            command="model service up <service_name> [ no_gpu ]}",
            description=f"""
Launch a model service after it was cataloged using the <cmd>catalog model service from ...</cmd> command.
You can launch your service with GPU disabled by appending <cmd>no_gpu</cmd>.

{MODEL_SERVICE_QUOTES}

Examples:
- <cmd>model service up gen</cmd>
- <cmd>model service up 'gen'</cmd>
- <cmd>model service up gen no_gpu</cmd>""",
        )
    )

    statements.append(
        py.Forward(
            model
            + service
            + local
            + up
            + service_name("service_name")
            + py.Optional(py.CaselessKeyword("NO_GPU")("no_gpu"))
        )("local_service_up")
    )
    help.append(
        help_dict_create(
            name="model local up",
            category="Model Service",
            command="model service local up <service_name> ",
            description=f"""
Launch a model service locally.

{MODEL_SERVICE_QUOTES}

Example:
<cmd>model service local up gen</cmd>
""",
        )
    )

    statements.append(py.Forward(model + service + down + service_name("service_name"))("service_down"))
    help.append(
        help_dict_create(
            name="model down",
            category="Model Service",
            command="model service down <service_name>",
            description=f"""
Shut down a model service.

{MODEL_SERVICE_QUOTES}

Examples:
- <cmd>model service down gen</cmd>
- <cmd>model service down 'gen'</cmd>
""",
        )
    )
