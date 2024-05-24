"""Implements Support for both Property and Generation Service Integration into OpenAD"""

# import openad_model_property_service.service_defs as new_prop_services

# from openad.core.help import help_dict_create
import requests
import os
import pandas as pd
import glob
import json
from openad.helpers.output import (
    output_error,
    output_text,
    output_success,
    output_warning,
)
from openad.helpers.spinner import spinner
from openad.openad_model_plugin.catalog_model_services import (
    help_dict_create,
    get_service_endpoint,
)


# from openad.molecules.mol_functions import MOL_PROPERTIES as m_props
# from openad.helpers.general import is_notebook_mode

from pyparsing import (
    Word,
    delimitedList,
    alphas,
    alphanums,
    OneOrMore,
    ZeroOrMore,
    CharsNotIn,
    Forward,
    CaselessKeyword,
    MatchFirst,
    Keyword,
    QuotedString,
    ParserElement,
    Suppress,
    Optional,
    Group,
    Combine,
    nums,
    oneOf,
    Literal,
    # replaceWith,
    # Combine,
    # pyparsing_test,
    # ParseException,
)

(
    get,
    lister,
    description,
    using,
    create,
    s_et,
    unset,
    workspace,
    workspaces,
    context,
    jobs,
    e_xec,
    a_s,
    optimize,
    w_ith,
    toolkits,
    toolkit,
    GPU,
    experiment,
    add,
    run,
    save,
    runs,
    show,
    file,
    d_isplay,
    history,
    data,
    remove,
    result,
    f_rom,
    inchi,
    inchikey,
    smiles,
    formula,
    name,
    l_ast,
    load,
    results,
    export,
    create,
    rename,
    merge,
    pubchem,
    sources,
    basic,
    force,
    append,
    only,
    upsert,
) = map(
    CaselessKeyword,
    "get list description using create set unset workspace workspaces context jobs exec\
          as optimize with toolkits toolkit gpu experiment add run save runs show \
              file display history data remove result from inchi inchikey smiles formula name last load results export create rename merge pubchem sources basic force append only upsert".split(),
)
name_expr = Word(alphanums + "_" + ".")
key_val_expr = Word(alphanums + "_" + ".")
key_val_expr_num = Word(nums)
key_val_expr_alpha = Word(alphanums + "_" + ".")
number_type = Combine((Word(nums) + "." + Word(nums))) | Word(nums)
key_val_line = Group(name_expr("key") + Suppress("=") + key_val_expr("val"))
boolean_var = Keyword("True") | Keyword("False")

mol = ["molecule", "mol"]
mols = ["molecules", "mols"]
molset = ["molecule-set", "molset"]
molsets = ["molecule-sets", "molsets"]
clear = CaselessKeyword("clear")
cache = CaselessKeyword("cache")
analysis = CaselessKeyword("analysis")
enrich = CaselessKeyword("enrich")
mol_properties = ["synonyms"]
# mol_properties.extend(m_props)
mol_properties = MatchFirst(map(CaselessKeyword, mol_properties))
molecules = MatchFirst(map(CaselessKeyword, mols))
molecule = MatchFirst(map(CaselessKeyword, mol))
molecule_set = MatchFirst(map(CaselessKeyword, molset))
molecule_sets = MatchFirst(map(CaselessKeyword, molsets))
numbers = (
    Word(nums + "." + nums)
    | Word(nums)
    | Word("-" + nums + "." + nums)
    | Word("-" + nums)
    | Word("+" + nums + "." + nums)
    | Word("+" + nums)
)
molecule_identifier = (
    Word(
        alphanums + "_" + "[" + "]" + "(" + ")" + "=" + "," + "-" + "+" + "/" + "\\" + "#" + "@" + "." + "*" + ";",
    )
    | Word(
        alphas,
        alphanums + "_" + "[" + "]" + "(" + ")" + "=" + "," + "-" + "+" + "/" + "\\" + "#" + "@" + "." + "*" + ";",
    )
    | Word(nums)
) | Suppress(Word("'")) + (
    Word(
        alphanums + "_" + "[" + "]" + "(" + ")" + "=" + "," + "-" + "+" + "/" + "\\" + "#" + "@" + "." + "*" + ";",
    )
    | Word(
        alphas,
        alphanums + "_" + "[" + "]" + "(" + ")" + "=" + "," + "-" + "+" + "/" + "\\" + "#" + "@" + "." + "*" + ";",
    )
    | Word(nums)
) + Suppress(
    Word("'")
)


desc = QuotedString("'", escQuote="\\")
input_object = QuotedString('"', end_quote_char='"', escQuote="\\")

####################################################################################
# here we are setting out the help and pyparsing grammar string components for creating grammar from memtadata
# service_command_start  is used to contruct the start of a string and is broken up by type of service
#       # get_molecule_property is used for properties of molecules defined by smiles strings, this is restrictied to valid smiles characters.
#       # get_protein_property is for protein properties, this uses a different dtring format and currently we just take any alphanum
#       # get_crystal_property is for crystaline properties and takes descriptors from files, yet to determine if this is strategic or not
#       # generate_data is for any of the different types of data set generation options

save_as_clause = "+ Optional(CaselessKeyword('save_as')('save_as')+desc('results_file'))"
save_as_clause_help = " (save_as '<filename.csv>')"

service_command_start = {}
service_command_subject = {}
service_command_help = {}

service_command_start["get_molecule_property"] = 'get + CaselessKeyword("molecule") + CaselessKeyword("property")'
service_command_start["get_crystal_property"] = 'get + CaselessKeyword("crystal") + CaselessKeyword("property")'
service_command_start["get_protein_property"] = 'get + CaselessKeyword("protein") + CaselessKeyword("property")'
service_command_start["generate_data"] = 'CaselessKeyword("generate") + CaselessKeyword("with")'

service_command_subject[
    "get_molecule_property"
] = '+CaselessKeyword("for")+((Word("[")+delimitedList(molecule_identifier,delim=",")("molecules")+Word("]")|molecule_identifier("molecule")))'
service_command_subject[
    "get_protein_property"
] = '+CaselessKeyword("for")+((Word("[")+ delimitedList(molecule_identifier,delim=",")("proteins")+Word("]")|molecule_identifier("protein")))'
service_command_subject[
    "get_crystal_property"
] = '+CaselessKeyword("for")+((Word("[")+ delimitedList(desc,delim=",")("crystal_files")+Word("]")|desc("crystal_file")("crystal_PATH")))'
service_command_subject[
    "generate_data"
] = '+CaselessKeyword("data")+<TARGET>Optional(CaselessKeyword("Sample")+Word(nums)("sample_size"))'

###################################################################
# targets for generate Data
generation_targets = {
    "conditional_generation": {
        "object": "CaselessKeyword('for')+input_object('target@object')+",
        "string": "CaselessKeyword('for')+(molecule_identifier|desc)('target@string')+",
        "list": "CaselessKeyword('for')+Word('[')+delimitedList(desc|numbers)('target@list')+Word(']')+",
        "array": "CaselessKeyword('for')+Word('[')+delimitedList(desc|numbers)('target@list')+Word(']')+",
    },
    "prediction": {
        # "object": "CaselessKeyword('for')+input_object('target@object')+",
        "object": "CaselessKeyword('for')+(molecule_identifier|desc)('target@string')+",
        "string": "CaselessKeyword('for')+(molecule_identifier|desc)('target@string')+",
        "list": "CaselessKeyword('for')+Word('[')+delimitedList(desc)('target@list')+Word(']')+",
        "array": "CaselessKeyword('for')+Word('[')+delimitedList(desc)('target@list')+Word(']')+",
        "number": "CaselessKeyword('for')+number_type('target@snumber')+",
    },
    "generation": {
        "object": "Optional(CaselessKeyword('for')+input_object('target@object'))+",
        "string": "Optional(CaselessKeyword('for')+(molecule_identifier|desc)('target@string'))+",
        "list": "Optional(CaselessKeyword('for')+Word('[')+delimitedList(desc|numbers)('target@list')+Word(']'))+",
        "array": "Optional(CaselessKeyword('for')+Word('[')+delimitedList(desc|numbers)('target@list')+Word(']'))+",
    },
    "controlled_sampling": {
        "object": "Optional(CaselessKeyword('for')+input_object('target@object'))+",
        "string": "Optional(CaselessKeyword('for')+ (molecule_identifier|desc)('target@string')+Word(']'))+",
        "list": "Optional(CaselessKeyword('for')+Word('[')+delimitedList(desc)('target@list')+Word(']'))+",
        "number": "CaselessKeyword('for')+number_type('target@number')+",
        "array": "CaselessKeyword('for')+delimitedList(desc)('target@list')+",
    },
}
# """ algorithm_version="solubility", search="sample", temperature=1.5, tolerance=60,
#         sampling_wrapper={'fraction_to_mask': mask, 'property_goal': {'<esol>': 0.234}}"""


service_command_help[
    "get_molecule_property"
] = "get molecule property <property> for [<list of SMILES>] | <SMILES>   USING (<parameter>=<value> <parameter>=<value>)"
service_command_help[
    "get_crystal_property"
] = "get crystal property <property> for <directory>   USING (<parameter>=<value> <parameter>=<value>)"
service_command_help[
    "get_protein_property"
] = "get protein property <property> for [<list of Proteins>] | <Protein>   USING (<parameter>=<value> <parameter>=<value>)"
service_command_help[
    "generate_data"
] = "generate with <property> data <TARGET> (sample <sample_size>)  USING (<parameter>=<value> <parameter>=<value>) "


def service_grammar_add(statements: list, help: list, service_catalog: dict):
    """defines the grammar available for managing molecules"""
    for service in service_catalog.keys():
        service_list = service_catalog[service]
        for schema in service_list:
            command = "CaselessKeyword(service)('service')+" + service_command_start[schema["service_type"]]
            valid_types = None  # noqa: F841
            valid_type = None

            # in properties this refers to generateable properties for a Target in Generators this is simple a single name of a generarion Algorithm
            # for some propertyy statements there can be multiple properties in a single statement

            if len(list(schema["valid_types"])) > 1:
                valid_types = list(schema["valid_types"])  # noqa: F841  # Used in eval statement
            else:
                valid_type = f'CaselessKeyword("{list(schema["valid_types"])[0]}")("type")'
                help_type = list(schema["valid_types"])[0]

            if valid_type is None:
                valid_type = '( (Word("[")+delimitedList(oneOf(valid_types)|Suppress(Word("\'"))+oneOf(valid_types)+Suppress(Word("\'")),delim=",")("types")+Word("]")) | ( oneOf(valid_types)("type")) ) '
                help_type = "[ " + ", ".join(list(schema["valid_types"])) + " ] | <valid_type>  "
            expression = ""

            # if parameters  exist for command build parameter grammar
            if len(list(schema["parameters"])) > 0:
                if len(list(schema["required_parameters"])) > 0:
                    expression = "+"
                else:
                    expression = "+ Optional"

                expression = (
                    expression
                    + '( (CaselessKeyword ("USING")+ Suppress("(") +'
                    + optional_parameter_list(schema, "parameters")
                    + '+Suppress(")") )("USING"))'
                )
            # prepare command for generator type
            if "generator_type" in schema.keys():
                if schema["target"]:
                    target_type = schema["target"]["type"]
                    try:
                        cmd_subject = str(service_command_subject[schema["service_type"]]).replace(
                            "<TARGET>",
                            generation_targets[schema["generator_type"]["algorithm_type"]][target_type],
                        )
                    except Exception as e:
                        print(schema)
                        output_error(e)
                        continue
                else:
                    cmd_subject = str(service_command_subject[schema["service_type"]]).replace("<TARGET>", "")
            else:
                cmd_subject = service_command_subject[schema["service_type"]]

            # below is simply for when debugging is required
            """"
            print(
                "Forward( "
                + command
                + "+"
                + valid_type
                + cmd_subject
                + expression
                + ")"
                + f'("{schema["service_name"]}@{schema["service_type"]}")'
            )
            print("--------------------------------------------------------------------")
            """

            # Compile the pyparsing grammar for the statement
            try:
                stmt = eval(
                    "Forward( "
                    + command
                    + "+"
                    + valid_type
                    + cmd_subject
                    + expression
                    + save_as_clause
                    + ")"
                    + f'("{schema["service_name"]}@{schema["service_type"]}")'
                )
            except:
                print(2)
                output_error("error for schema")
                output_error(schema)
                continue

            # Add pyparsing to statament grammar array for pyparsing.
            statements.append(stmt)

            ## the following generates the Help statements for a given command
            try:
                target_description = ""
                function_description = ""
                if "target" in schema:
                    if schema["target"] is not None:
                        target_description = "<h2>Target:</h2>\r"
                        for key, value in schema["target"].items():
                            target_description = target_description + f"- <cmd>{key}</cmd> : {value}\n  "

                if schema["description"] is not None:
                    function_description = "\r<h2>Function Description:</h2>\r" + schema["description"]
                    while "  " in function_description:
                        function_description = function_description.replace("  ", " ")
            except Exception as e:
                output_error(e)

            parameter_help = "<h2>Parameters:</h2>"
            num_params = 0
            for parameter, description in dict(schema["parameters"]).items():
                num_params += 1
                print_description = ""
                for key, value in description.items():
                    print_description = print_description + f"- <cmd>{key}</cmd> : {value}\n  "

                parameter_help = parameter_help + f"\n<cmd>{parameter}</cmd> \r {print_description}\n  "

            required_parameters = ""
            for i in schema["required_parameters"]:
                if required_parameters == "":
                    required_parameters = "\n<h2>Required Parameters:</h2> \n"
                required_parameters = required_parameters + f"\n - <cmd>{i}</cmd>"
            algo_versions = ""
            if "algorithm_versions" in schema:
                algo_versions = " \n <h2> Algorithm Versions </h2> \n"
                for i in schema["algorithm_versions"]:
                    algo_versions = algo_versions + f"\n - <cmd>{i}</cmd>"
            try:
                if "generator_type" in schema.keys():
                    if not schema["target"]:
                        command_str = (
                            str(service + " " + service_command_help["generate_data"])
                            .replace("<TARGET>", "")
                            .replace("<property>", help_type)
                        )
                    else:
                        if schema["generator_type"]["algorithm_type"] in [
                            "conditional_generation",
                            "controlled_sampling",
                        ]:
                            command_str = (
                                service
                                + " "
                                + str(service_command_help["generate_data"])
                                .replace(
                                    "<TARGET>",
                                    " for  <" + schema["target"]["type"] + ">",
                                )
                                .replace("<property>", help_type)
                            )
                        else:
                            command_str = (
                                service
                                + " "
                                + str(service_command_help["generate_data"])
                                .replace(
                                    "<TARGET>",
                                    " for (<" + schema["target"]["type"] + ">)",
                                )
                                .replace("<property>", help_type)
                            )
                else:
                    command_str = str(service + " " + service_command_help[schema["service_type"]]).replace(
                        "<property>", help_type
                    )
            except Exception as e:
                output_error("-------")
                output_error(e)
                output_error("-------")
            if "generator_type" in schema:
                category = schema["generator_type"]["algorithm_type"]
            else:
                category = "Model-" + schema["sub_category"]

            # if command has no parameters simply remove the USING Clause
            if num_params == 0:
                command_str = command_str.replace(" USING (<parameter>=<value> <parameter>=<value>)", "")

            # add help statement to help array psed through into function
            help.append(
                help_dict_create(
                    name=schema["service_type"],
                    category=service + "->" + category,
                    parent=None,
                    command=command_str + save_as_clause_help,
                    description=target_description
                    + parameter_help
                    + algo_versions
                    + required_parameters
                    + function_description,
                )
            )

    return statements


def optional_parameter_list(inp_statement: dict, clause: str):
    """Create an optional parameter list for a clause"""
    ii = 0
    expression = " "
    for i in inp_statement[clause]:
        if "allOf" in inp_statement[clause][i] and "type" not in inp_statement[clause][i]:
            type_str = "allOf"
        elif "anyOf" in inp_statement[clause][i] and "type" not in inp_statement[clause][i]:
            type_str = "anyOf"
        elif "type" in inp_statement[clause][i]:
            type_str = "type"
        elif "allOf" in inp_statement[clause][i]:
            type_str = "type"
        if isinstance(inp_statement[clause][i][type_str], list):
            if isinstance(i, int):
                parameter = "param_" + i + "@" + "integer"
            elif isinstance(i, float):
                parameter = "param_" + i + "@" + "float"
            else:
                parameter = "param_" + i + "@" + "other"
        else:
            parameter = "param_" + i + "@" + inp_statement[clause][i][type_str]

        if i in inp_statement["required_parameters"]:
            status = ""

        else:
            status = "ZeroOrMore"
        if ii == 0:
            expression = expression + " "
            if inp_statement[clause][i] == "allOf":
                expression = (
                    expression
                    + f" {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i][type_str] == "string":
                expression = (
                    expression
                    + f" {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i][type_str] == "desc":
                expression = (
                    expression
                    + f" {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+desc('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i][type_str] == "object":
                expression = (
                    expression
                    + f" {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+input_object('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i][type_str] == "boolean":
                expression = (
                    expression
                    + f" {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+boolean_var('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            else:
                expression = (
                    expression
                    + f" {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+number_type('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
        else:
            if type_str == "allOf":
                expression = (
                    expression
                    + f" & {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )

            elif inp_statement[clause][i][type_str] == "string":
                expression = (
                    expression
                    + f" & {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i][type_str] == "desc":
                expression = (
                    expression
                    + f" & {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+desc('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i][type_str] == "object":
                expression = (
                    expression
                    + f" & {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+input_object('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i][type_str] == "boolean":
                expression = (
                    expression
                    + f" & {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+boolean_var('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            else:
                expression = (
                    expression
                    + f" & {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+number_type('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
        ii = 1

    return expression


def subject_files_repository(file_directory, suffix):
    file_directory = os.path.expanduser(file_directory)
    parameter_files = glob.glob(file_directory + "/*." + suffix)
    parameter_list_files = []
    for file in parameter_files:
        with open(file, "r", encoding="utf-8") as file_handle:
            name = os.path.basename(file)
            payload = file_handle.read()
            parameter_list_files.append([name, payload])

    return parameter_list_files


def request_generate(request_input):
    """This function constructs the request to be passed to the remote Server"""

    name = request_input.getName()
    Sample_Size = None
    subjects = []
    if name.split("@")[1] == "generate_data":
        if "sample_size" in request_input.as_dict():
            Sample_Size = request_input.as_dict()["sample_size"]
        if "target@object" in request_input.as_dict():
            subjects = [json.loads(str(request_input.as_dict()["target@object"]).replace("'", '"'))]
        if "target@string" in request_input.as_dict():
            subjects = request_input.as_dict()["target@string"]
        if "target@number" in request_input.as_dict():
            subjects = request_input.as_dict()["target@number"]

        if "target@list" in request_input.as_dict():
            subjects = request_input.as_dict()["target@list"]

    if name.split("@")[1] == "get_molecule_property":
        if "molecules" in request_input.as_dict():
            subjects = request_input.as_dict()["molecules"]

        if "molecule" in request_input.as_dict():
            if isinstance(request_input.as_dict()["molecule"], list):
                subjects = request_input.as_dict()["molecule"]
            else:
                subjects = [request_input.as_dict()["molecule"]]
    if name.split("@")[1] == "get_protein_property":
        if "proteins" in request_input.as_dict():
            subjects = request_input.as_dict()["proteins"]
        if "protein" in request_input.as_dict():
            if isinstance(request_input.as_dict()["protein"], list):
                subjects = request_input.as_dict()["protein"]
            else:
                subjects = [request_input.as_dict()["protein"]]

    if "types" in request_input.as_dict():
        property_types = request_input.as_dict()["types"]
    if "type" in request_input.as_dict():
        property_types = [request_input.as_dict()["type"]]
    if name.split("@")[1] == "get_crystal_property":
        subjects = subject_files_repository(request_input.as_dict()["crystal_PATH"], "cif")
        subjects.extend(subject_files_repository(request_input.as_dict()["crystal_PATH"], "csv"))
    template = {
        "service_name": request_input.getName().split("@")[0],
        "service_type": request_input.getName().split("@")[1],
        "parameters": {
            "property_type": property_types,
            "subjects": subjects,
        },
        "api_key": "api-dthgwrhrthrtrth",
    }
    if Sample_Size is not None:
        template["sample_size"] = Sample_Size

    for param in request_input.as_dict().keys():
        if str(param).startswith("param_"):
            actual_param = str(param)[6:]
            if actual_param.split("@")[1] == "number" or actual_param.split("@")[1] == "float":
                template["parameters"][actual_param.split("@")[0]] = float(
                    "".join(list(request_input.as_dict()[param]["val"]))
                )
            elif actual_param.split("@")[1] == "integer":
                template["parameters"][actual_param.split("@")[0]] = int(
                    "".join(list(request_input.as_dict()[param]["val"]))
                )
            elif actual_param.split("@")[1] == "object":
                # x = dict(json.loads(request_input.as_dict()[param]["val"].replace("'", '"')))
                template["parameters"][actual_param.split("@")[0]] = json.loads(
                    request_input.as_dict()[param]["val"].replace("'", '"')
                )
            elif actual_param.split("@")[1] == "integer":
                template["parameters"][actual_param.split("@")[0]] = bool(request_input.as_dict()[param]["val"])

            else:
                template["parameters"][actual_param.split("@")[0]] = request_input.as_dict()[param]["val"]
    return template


def convert(lst):
    """Used for for converting lists to strings."""
    return str(lst).translate("[],'")


def openad_model_requestor(cmd_pointer, parser):
    """The Procedure handles communication with external services"""
    if "service" in parser.as_dict():
        service = parser.as_dict()["service"]
    else:
        service = None

    a_request = request_generate(parser)
    Endpoint = get_service_endpoint(service)

    if Endpoint is not None and len((Endpoint.strip())) > 0:
        if "http" not in Endpoint:
            Endpoint = "http://" + Endpoint
    # Endpoint = "http://34.205.69.8:8080"
    else:
        Endpoint = None

    if Endpoint is None:
        return output_error(
            "No Service Cataloged or service not up. \n Check Service Status <cmd>model service status</cmd> "
        )

    spinner.start("Executing Request Against Server")

    try:
        response = requests.post(Endpoint + "/service", json=a_request)
    except Exception as e:
        spinner.fail("Request Failed")
        spinner.stop()
        output_error(str(e))
        return output_error("Error: \n Server not reachable at " + str(Endpoint))

    spinner.succeed("Request Returned")
    spinner.stop()
    try:
        response_result = response.json()
        try:
            if isinstance(response_result, dict):
                if "error" in response_result:
                    run_error = "Request Error:\n"

                    for key, value in response_result["error"].items():
                        value = str(value).replace("<", "`<")
                        value = str(value).replace(">", ">`")
                        run_error = run_error + f"- <cmd>{key}</cmd> : {value}\n  "
                    return output_error(run_error)
                if "detail" in response_result:
                    return output_warning(response_result["detail"])

            result = pd.DataFrame(response_result)
            if "save_as" in parser:
                results_file = str(parser["results_file"])
                if not results_file.endswith(".csv"):
                    results_file = results_file + ".csv"
                result.to_csv(
                    cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + results_file,
                    index=False,
                )
            return result

        except:
            result = response_result

        if isinstance(result, dict):
            if "error" in result:
                run_error = "Request Error:\n"
                for key, value in result["error"].items():
                    run_error = run_error + f"- <cmd>{key}</cmd> : {value}\n  "
                return output_text(run_error)

    except Exception as e:
        run_error = "HTTP Request Error:\n"
        spinner.fail("Request Failed")
        spinner.stop()
        return output_error(run_error + "\n" + str(e))

    return result
