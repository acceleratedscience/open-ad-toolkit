"""Implements Support for both Property and Generation Service Integration into OpenAD"""

# import openad_model_property_service.service_defs as new_prop_services

import glob
import json
import os
import shutil, re
import openad.helpers.general as helpers_general
import pandas as pd

# from openad.core.help import help_dict_create
import requests

from openad.helpers.output import output_error, output_success, output_text, output_warning
from openad.helpers.spinner import Spinner
from openad.openad_model_plugin.catalog_model_services import get_service_requester, help_dict_create
from openad.openad_model_plugin.auth_services import get_service_api_key
from openad.openad_model_plugin.catalog_model_services import Dispatcher
from openad.app.global_var_lib import GLOBAL_SETTINGS
from openad.smols.smol_batch_files import merge_molecule_property_data
from pyparsing import (  # replaceWith,; Combine,; pyparsing_test,; ParseException,
    CaselessKeyword,
    CharsNotIn,
    Combine,
    Forward,
    Group,
    Keyword,
    Literal,
    MatchFirst,
    OneOrMore,
    Optional,
    ParserElement,
    QuotedString,
    Suppress,
    Word,
    ZeroOrMore,
    alphanums,
    alphas,
    delimitedList,
    nums,
    oneOf,
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
desc = QuotedString("'", escQuote="\\")
name_expr = Word(alphanums + "_" + ".")
key_val_expr = Word(alphanums + "_" + ".") | desc
key_val_expr_num = Word(nums)
key_val_expr_alpha = Word(alphanums + "_" + ".")
number_type = Combine(Optional("-") + Word(nums) + Word(".") + Word(nums)) | Word(nums)
array_var = Suppress(Word("[")) + delimitedList(OneOrMore(desc | number_type)) + Suppress(Word("]"))
key_val_line = Group(name_expr("key") + Suppress("=") + key_val_expr("val"))
boolean_var = Keyword("True") | Keyword("False")

mol = ["molecule", "mol"]
mols = ["molecules", "mols"]
molset = ["molecule-set", "molset"]
molsets = ["molecule-sets", "molsets"]
clear = CaselessKeyword("clear")
cache = CaselessKeyword("cache")
molecules_list = ["@molecules", "@mols"]
analysis = CaselessKeyword("analysis")
enrich = CaselessKeyword("enrich")
mol_properties = ["synonyms"]
# mol_properties.extend(m_props)
mol_list = MatchFirst(map(CaselessKeyword, molecules_list))
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

CLI_WIDTH = helpers_general.get_print_width(full=True)

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
async_clause = "+ Optional(CaselessKeyword('async')('async')) "
async_clause_help = " (async)"
service_command_start = {}
service_command_subject = {}
service_command_help = {}
service_command_description = {}
service_command_merge = {}

service_command_start["get_molecule_property"] = 'get + CaselessKeyword("molecule") + CaselessKeyword("property")'
service_command_start["get_crystal_property"] = 'get + CaselessKeyword("crystal") + CaselessKeyword("property")'
service_command_start["get_protein_property"] = 'get + CaselessKeyword("protein") + CaselessKeyword("property")'
service_command_start["generate_data"] = 'CaselessKeyword("generate") + CaselessKeyword("with")'

service_command_merge[
    "get_molecule_property"
] = '+ Optional((CaselessKeyword("merge with mols")|CaselessKeyword("merge with molecules"))("merge_with_mws"))'
service_command_merge["get_crystal_property"] = ""
service_command_merge["get_protein_property"] = ""
service_command_merge["generate_data"] = ""

service_command_subject[
    "get_molecule_property"
] = '+CaselessKeyword("for")+(mol_list("mol_list")|(Word("[")+delimitedList(molecule_identifier,delim=",")("molecules")+Word("]")|molecule_identifier("molecule")))'
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
] = "get molecule property <property> FOR @mols | [<list of SMILES>] | <SMILES>   USING (<parameter>=<value> <parameter>=<value>) (merge with mols|molecules)"
service_command_help[
    "get_crystal_property"
] = "get crystal property <property> FOR <directory> USING (<parameter>=<value> <parameter>=<value>)"
service_command_help[
    "get_protein_property"
] = "get protein property <property> FOR [<list of Proteins>] | <Protein> USING (<parameter>=<value> <parameter>=<value>)"
service_command_help[
    "generate_data"
] = "generate with <property> data <TARGET> (sample <sample_size>) USING (<parameter>=<value> <parameter>=<value>) "

service_command_description[
    "get_molecule_property"
] = """
This command gets (generate/predict) a molecules property for one or molecules specified with a SMILES string in the <cmd>FOR</cmd> clause. SMILES can be provided as a single SMILES string or multiple smiles in a comma seperated list in square brackets e.g. <cmd> FOR [CCO, CC(C)CC1=CC=C(C=C1)C(C)C(=O)O ] </cmd>.
SMILES strings can be specified with or without single quotes, but when in a list smiles with square brackets should be enclosed in single quotes e.g <cmd>[ 'C([H])([H])([H])[H]' ,CCO ]</cmd>

This command gets (generate/predict) the following properties:\n<cmd><property_list></cmd>

The clause <cmd>merge with mols </cmd> will merge the resulting molecule properties with the memory molecule working set.

Note: <cmd> @mols </cmd>  specifies the list in the current molecules working set.
    Example: 
        The following will generate the specified properties for the molecules in the molecule working set.
        <cmd>get molecule property <property> for @mols </cmd>
        The following will generate the specified properties for the molecules in the molecule working set and will merge the resulting molecule properties with the memory molecule working set.
        <cmd>get molecule property <property> for @mols merge with mols</cmd>

"""
service_command_description[
    "get_crystal_property"
] = """
This command gets (generate/predict) crystal properties
"""
service_command_description[
    "get_protein_property"
] = """
This command gets (generate/predict) a proteins property for one or protiens specified with a FASTA string in the <cmd>FOR</cmd> clause.
FASTA strings can be provided as a single  string or multiple FASTA strings in a comma seperated list in square brackets 
e.g. <cmd> FOR ['NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC','NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC' ]</cmd>.
FASTA strings must be provided in single quotes.
This command gets (generate/predict) the following properties:\n<cmd><property_list></cmd>\n
"""
service_command_description[
    "generate_data"
] = """
    This function generates a data set based on the following parameters 
 """

async_help_clause = "\n \n Note: If <cmd> async clause </cmd> is defined the user will be returned an id for the given job and will use the <cmd> `model service <service name> result '<job_id>' </cmd> command to retrieve it when it is ready. use this command to test for readiness."


def service_grammar_add(statements: list, help: list, service_catalog: dict):
    """defines the grammar available for managing molecules"""
    for service in service_catalog.keys():
        service_list = service_catalog[service]
        for schema in service_list:
            # Allow Async for command if supported
            if "async_allow" in schema and schema["async_allow"]:
                async_allow = True
            else:
                async_allow = False

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
                help_type = "[ " + ", ".join(list(schema["valid_types"])) + " ] | <valid_property>  "
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
                    + service_command_merge[schema["service_type"]]
                    + save_as_clause
                    + (async_clause if async_allow else "")
                    + ")"
                    + f'("{schema["service_name"]}@{schema["service_type"]}")'
                )

            except:
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
                        target_description = "<h2>Target:</h2>\n"
                        for key, value in schema["target"].items():
                            target_description = target_description + f"- <cmd>{key}</cmd> : {value}\n  "

                if schema["description"] is not None:
                    function_description = "\n<h2>Function Description:</h2>\n" + schema["description"]
                    while "  " in function_description:
                        function_description = function_description.replace("  ", " ")
            except Exception as e:
                output_error(e)
            parameter_help = ""
            num_params = 0
            for parameter, description in dict(schema["parameters"]).items():
                if parameter in ["selected_property", "property_type", "domain", "algorithm_type"]:
                    continue
                num_params += 1
                print_description = ""
                for key, value in description.items():
                    print_description = print_description + f"- <cmd>{key}</cmd> : {value}\n  "

                parameter_help = parameter_help + f"<cmd>{parameter}</cmd> \n {print_description}\n  "
            if "generator_type" in schema.keys():
                key = "generate_data"
            else:
                key = schema["service_type"]

            if num_params != 0:
                try:
                    parameter_help = (
                        service_command_description[key]
                        #    .replace("<property_list>", help_type.split("|")[0])
                        .replace(
                            "<property_list>",
                            format_properties(
                                list(
                                    help_type.split("|")[0]
                                    .replace(" ", "")
                                    .replace("[", "")
                                    .replace("]", "")
                                    .split(",")
                                )
                            ),
                        ).replace("<property>", help_type.split("|")[0].split(",")[0].replace("[", "").lstrip())
                        + "<h2>Parameters:</h2>\n   <warning>--Note: Parameters should be entered for <cmd> USING Clause </cmd> in the order they are below. </warning>\n"
                        + parameter_help
                    )
                except Exception as e:
                    print(e)
            else:
                parameter_help = (
                    service_command_description[key]
                    # .replace("<property_list>", help_type.split("|")[0])
                    .replace(
                        "<property_list>",
                        format_properties(
                            list(help_type.split("|")[0].replace(" ", "").replace("[", "").replace("]", "").split(","))
                        ),
                    ).replace("<property>", help_type.split("|")[0].split(",")[0].replace("[", "").lstrip())
                    + "\n"
                    + " <h2>No Parameters</h2>\n"
                    + parameter_help
                )

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
                    command=command_str + save_as_clause_help + (async_clause_help if async_allow else ""),
                    description=target_description
                    + parameter_help
                    + algo_versions
                    + required_parameters
                    + function_description
                    + (async_help_clause if async_allow else ""),
                )
            )

    return statements


def optional_parameter_list(inp_statement: dict, clause: str):
    """Create an optional parameter list for a clause"""
    ii = 0
    expression = " "
    type_dict = {
        "allOf": "key_val_expr",
        "anyOf": "key_val_expr",
        "string": "key_val_expr",
        "desc": "desc",
        "object": "input_object",
        "boolean": "boolean_var",
        "array": "array_var",
        "number": "number_type",
        "integer": "number_type",
    }
    for i in inp_statement[clause]:
        if i in ["selected_property", "property_type", "domain", "algorithm_type"]:
            continue
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
            if type_str in ["anyOf", "allOf"]:
                expression = (
                    expression
                    + f" {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i][type_str] in type_dict:
                expression = (
                    expression
                    + f" {status}(Group( CaselessKeyword ('"
                    + i
                    + f"') +Suppress('=')+{type_dict[inp_statement[clause][i][type_str]]}('val'))('"
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
            if type_str in ["anyOf", "allOf"]:
                expression = (
                    expression
                    + f" & {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
            elif inp_statement[clause][i][type_str] in type_dict:
                expression = (
                    expression
                    + f" & {status}(Group( CaselessKeyword ('"
                    + i
                    + f"') +Suppress('=')+{type_dict[inp_statement[clause][i][type_str]]}('val'))('"
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


def mol_list_gen(cmd_pointer):
    mol_list = []
    for molecule in cmd_pointer.molecule_list:
        mol_list.append(molecule["identifiers"]["canonical_smiles"])
    return mol_list


def request_generate(cmd_pointer, request_input):
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
        if "mol_list" in request_input.as_dict():
            subjects = mol_list_gen(cmd_pointer)

        elif "molecules" in request_input.as_dict():
            subjects = request_input.as_dict()["molecules"]

        elif "molecule" in request_input.as_dict():
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
        "api_key": "reserved for federated APIs",
    }
    if "async" in request_input.as_dict():
        template["async"] = True
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
            elif actual_param.split("@")[1] == "string":
                template["parameters"][actual_param.split("@")[0]] = str(request_input.as_dict()[param]["val"])
            else:
                template["parameters"][actual_param.split("@")[0]] = request_input.as_dict()[param]["val"]
    return template


def convert(lst):
    """Used for for converting lists to strings."""
    return str(lst).translate("[],'")


def openad_model_requestor(cmd_pointer, parser):
    """The Procedure handles communication with external services"""
    if "service" in parser.as_dict():
        service_name = parser.as_dict()["service"]
    else:
        service_name = None

    a_request = request_generate(cmd_pointer, parser)

    spinner = Spinner(GLOBAL_SETTINGS["VERBOSE"])
    spinner.start("Executing Request Against Server")

    with Dispatcher as servicer:
        service_status = servicer.get_short_status(service_name)
    try:
        # response = Dispatcher.service_request(
        #     name=service_name, method="POST", timeout=None, verify=not service_status.get("is_remote"), _json=a_request
        # )
        response = Dispatcher.service_request(
            name=service_name, method="POST", timeout=None, verify=False, _json=a_request
        )
        # response = requests.post(Endpoint + "/service", json=a_request, headers=headers, verify=False)
    except Exception as e:
        spinner.fail("Request Failed")
        spinner.stop()
        output_error(str(e))
        return output_error("Error: \n Server not reachable  " + str(service_status.get("url")))

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
            if "merge_with_mws" in parser.as_dict():
                merge_molecule_property_data(cmd_pointer=cmd_pointer, dataframe=result)
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


def format_properties(props: list):
    """formats synonyms for display"""

    props_string = ""

    prop_length = 10
    for i in props:
        if len(i) > prop_length:
            prop_length = len(i)

    props_string = props_string + single_value_columns(props, helpers_general.get_print_width(full=True), 40)

    props_string = re.sub(r"<(.*?:)> ", r"<success>\1</success>", props_string)
    return "\n" + props_string


def single_value_columns(values, sys_cli_width, designated_display_width):
    """displays columns of single value"""
    return_string = ""
    i = 0
    if GLOBAL_SETTINGS["display"] == "notebook":
        cli_width = 150
    else:
        cli_width = min(sys_cli_width, 150)

    for value in values:
        if len(str(value).strip()) == 0:
            continue
        display_width = designated_display_width
        spacing = 0
        while len(value) > display_width:
            display_width += designated_display_width
        if (len(f"{value:<{display_width}}") + i + spacing < cli_width) or return_string == "":
            if spacing == 0:
                return_string = return_string + f"{value:<{display_width}}"
            else:
                return_string = return_string + " " + f"{value:<{display_width}}"

            i = i + len(f"{value:<{display_width}}") + spacing
            spacing = 1
        else:
            return_string = return_string + f"\n{value:<{display_width}}"
            i = len(f"{value:<{display_width}}")
    return return_string
