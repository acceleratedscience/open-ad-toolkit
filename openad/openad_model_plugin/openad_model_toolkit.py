""" Implements Support for both Property and Generation Service Integration into OpenAD"""

# import openad_model_property_service.service_defs as new_prop_services

# from openad.core.help import help_dict_create
import requests
import os
import pandas as pd
import glob
import json
from openad.helpers.output import output_error, output_text, output_success
from openad.helpers.spinner import spinner
from openad.openad_model_plugin.catalog_model_services import help_dict_create, get_service_endpoint


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
    # Literal,
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
molecule_identifier = Word(
    alphas, alphanums + "_" + "[" + "]" + "(" + ")" + "=" + "," + "-" + "+" + "/" + "#" + "@" + "." + "*" + ";"
) | Word(nums)


desc = QuotedString("'", escQuote="\\")
input_object = QuotedString('"', end_quote_char='"', escQuote="\\")

service_command_start = {}
service_command_subject = {}
service_command_help = {}
service_command_start["get_molecule_property"] = 'get + CaselessKeyword("molecule") + CaselessKeyword("property")'
service_command_start["get_crystal_property"] = 'get + CaselessKeyword("crystal") + CaselessKeyword("property")'
service_command_start["get_protein_property"] = 'get + CaselessKeyword("protein") + CaselessKeyword("property")'
service_command_start["generate_data"] = 'CaselessKeyword("generate") + CaselessKeyword("with")'

service_command_subject["get_molecule_property"] = (
    '+CaselessKeyword("for")+((Word("[")+ delimitedList(molecule_identifier,delim=",")("molecules")+Word("]")|molecule_identifier("molecule")))'
)
service_command_subject["get_protein_property"] = (
    '+CaselessKeyword("for")+((Word("[")+ delimitedList(molecule_identifier,delim=",")("proteins")+Word("]")|molecule_identifier("protein")))'
)
service_command_subject["get_crystal_property"] = (
    '+CaselessKeyword("for")+((Word("[")+ delimitedList(desc,delim=",")("crystal_files")+Word("]")|desc("crystal_file")("crystal_PATH")))'
)
service_command_subject["generate_data"] = (
    '+CaselessKeyword("data")+<TARGET>Optional(CaselessKeyword("Sample")+Word(nums)("sample_size"))'
)

###################################################################
# targets for generate Data
generation_targets = {
    "conditional_generation": {
        "object": "CaselessKeyword('for')+input_object('target@object')+",
        "string": "CaselessKeyword('for')+(molecule_identifier|desc)('target@string')+",
        "list": "CaselessKeyword('for')+delimitedList(desc)('target@list')+",
    },
    "generation": {
        "object": "Optional(input_object('target@object'))+",
        "string": "Optional(CaselessKeyword('for')+(molecule_identifier|desc)('target@string'))+",
        "list": "Optional(CaselessKeyword('for')+delimitedList(desc))('target@list')+",
    },
    "controlled_generation": {
        "object": "Optional(CaselessKeyword('for')+input_object('target@object'))+",
        "string": "Optional(CaselessKeyword('for')+ (molecule_identifier|desc)('target@string'))+",
        "list": "Optional(CaselessKeyword('for')+delimitedList(desc))('target@list')+",
    },
}
""" algorithm_version="solubility", search="sample", temperature=1.5, tolerance=60,
        sampling_wrapper={'fraction_to_mask': mask, 'property_goal': {'<esol>': 0.234}}"""


service_command_help["get_molecule_property"] = (
    "get molecule property <property> for [<list of SMILES>] | <SMILES>   USING (<parameter>=<value> <parameter>=<value>)"
)
service_command_help["get_crystal_property"] = (
    "get crystal property <property> for <directory>   USING (<parameter>=<value> <parameter>=<value>)"
)
service_command_help["get_protein_property"] = (
    "get protein property <property> for [<list of Proteins>] | <Protein>   USING (<parameter>=<value> <parameter>=<value>)"
)
service_command_help["generate_data"] = (
    "generate with <property> data <TARGET> (sample <sample_size>)  USING (<parameter>=<value> <parameter>=<value>) "
)

def service_grammar_add(statements: list, help: list, service_catalog: dict):
    """defines the grammar available for managing molecules"""
    for service in service_catalog.keys():
        service_list = service_catalog[service]
        for schema in service_list:
            command = f"CaselessKeyword(service)('service')+" + service_command_start[schema["service_type"]]
            valid_types = None
            valid_type = None
            first = True

            if len(list(schema["valid_types"])) > 1:
                valid_types = list(schema["valid_types"])  # Useed in eval statement
            else:
                valid_type = f'CaselessKeyword("{list(schema["valid_types"])[0]}")("type")'
                help_type = list(schema["valid_types"])[0]

            if valid_type is None:
                valid_type = f'( (Word("[")+delimitedList(oneOf(valid_types),delim=",")("types")+Word("]")) | ( oneOf(valid_types)("type")) ) '
                help_type = "[ " + ", ".join(list(schema["valid_types"])) + " ] | <valid_type>  "
            expression = ""

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

            if "generator_type" in schema.keys():
                if schema["target"]:
                    target_type = schema["target"]["type"]
                    try:
                        cmd_subject = str(service_command_subject[schema["service_type"]]).replace(
                            "<TARGET>", generation_targets[schema["generator_type"]["algorithm_type"]][target_type]
                        )
                    except Exception as e:
                        print(e)
                        continue
                else:
                    cmd_subject = str(service_command_subject[schema["service_type"]]).replace("<TARGET>", "")
            else:
                cmd_subject = service_command_subject[schema["service_type"]]

            try:
                stmt = eval(
                    "Forward( "
                    + command
                    + "+"
                    + valid_type
                    + cmd_subject
                    + expression
                    + ")"
                    + f'("{schema["service_name"]}@{schema["service_type"]}")'
                )
            except:
                print("error for schema")
                print(schema)
                continue

            statements.append(stmt)
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
                print(e)

            parameter_help = "<h2>Parameters:</h2>"
            for parameter, description in dict(schema["parameters"]).items():
                print_description = ""
                for key, value in description.items():
                    print_description = print_description + f"- <cmd>{key}</cmd> : {value}\n  "

                parameter_help = parameter_help + f"\n<cmd>{parameter}</cmd> \r {print_description}\n  "

            required_parameters = ""
            for i in schema["required_parameters"]:
                if required_parameters == "":
                    required_parameters = "\n<h2>Required Parameters:<\h2> \n"
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
                        if schema["generator_type"]["algorithm_type"] == "conditional_generation":
                            command_str = (
                                service
                                + " "
                                + str(service_command_help["generate_data"])
                                .replace("<TARGET>", " for  <" + schema["target"]["type"] + ">")
                                .replace("<property>", help_type)
                            )
                        else:
                            command_str = (
                                service
                                + " "
                                + str(service_command_help["generate_data"])
                                .replace("<TARGET>", " for (<" + schema["target"]["type"] + ">)")
                                .replace("<property>", help_type)
                            )
                else:
                    command_str = str(service + " " + service_command_help[schema["service_type"]]).replace(
                        "<property>", help_type
                    )
            except Exception as e:
                print("-------")
                print(e)
                print("-------")
            if "generator_type" in schema:
                category = schema["generator_type"]["algorithm_type"]
            else:
                category = "Model-" + schema["sub_category"]
            help.append(
                help_dict_create(
                    name=schema["service_type"],
                    category=service + "->" + category,
                    parent=None,
                    command=command_str,
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
            subjects = [json.loads(request_input.as_dict()["target@object"].replace("'", '"'))]
        if "target@string" in request_input.as_dict():
            subjects = request_input.as_dict()["target@string"]
        if "target@list" in request_input.as_dict():
            subjects = request_input.as_dict()["target@list"]

    if name.split("@")[1] == "get_molecule_property":
        if "molecules" in request_input.as_dict():
            subjects = request_input.as_dict()["molecules"]

        if "molecule" in request_input.as_dict():
            subjects = [request_input.as_dict()["molecule"]]
    if name.split("@")[1] == "get_protein_property":
        if "proteins" in request_input.as_dict():
            subjects = request_input.as_dict()["proteins"]
        if "protein" in request_input.as_dict():
            subjects = [request_input.as_dict()["protein"]]

    if "types" in request_input.as_dict():
        property_types = request_input.as_dict()["types"]
    if "type" in request_input.as_dict():
        property_types = [request_input.as_dict()["type"]]
    if name.split("@")[1] == "get_crystal_property":
        subjects = subject_files_repository(request_input.as_dict()["crystal_PATH"], "cif")

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

    if Endpoint is None:
        output_error("No Service Cataloged")
        return None

    spinner.start("Executing Request Against Server")

    try:
        response = requests.post(Endpoint + "/service", json=a_request)
    except:
        spinner.fail("Request Failed")
        spinner.stop()
        return output_error("Error Server not reachable")

    spinner.succeed("Request Returned")
    spinner.stop()
    try:

        response_result = response.json()
        try:

            if isinstance(response_result, dict):
                if "error" in response_result:

                    run_error = "Request Error:\n"
                    for key, value in response_result["error"].items():
                        run_error = run_error + f"- <cmd>{key}</cmd> : {value}\n  "
                    return output_error(run_error)

            result = pd.DataFrame(response_result)
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
        spinner.fail("Request Failed")
        spinner.stop()
        return output_error(run_error + "\n" + str(e))

    return result
