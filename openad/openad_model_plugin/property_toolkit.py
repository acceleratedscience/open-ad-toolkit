# import openad_model_property_service.service_defs as new_prop_services

# from openad.core.help import help_dict_create
import requests
import pandas as pd

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

service_command_start = {}
service_command_subject = {}
service_command_help = {}
service_command_start["get_molecule_property"] = 'get + CaselessKeyword("molecule") + CaselessKeyword("property")'
service_command_start["get_crystal_property"] = 'get + CaselessKeyword("crystal") + CaselessKeyword("property")'
service_command_start["get_protein_property"] = 'get + CaselessKeyword("protein") + CaselessKeyword("property")'
service_command_subject["get_molecule_property"] = (
    '+CaselessKeyword("for")+((Word("[")+ OneOrMore(molecule_identifier)("molecules")+Word("]")|molecule_identifier("molecule")))'
)
service_command_subject["get_protein_property"] = (
    '+CaselessKeyword("for")+((Word("[")+ OneOrMore(molecule_identifier)("proteins")+Word("]")|molecule_identifier("protein")))'
)
service_command_subject["get_crystal_property"] = (
    '+CaselessKeyword("for")+((Word("[")+ OneOrMore(desc)("crystal_files")+Word("]")|desc("crystal_file")("crystal_PATH")))'
)

service_command_help["get_molecule_property"] = (
    "get molecule property <property> for [<list of SMILES>] | <SMILES>   USING (<parameter>=<value>)"
)
service_command_help["get_crystal_property"] = (
    "get crystal property <property> for <directory>   USING (<parameter>=<value>)"
)
service_command_help["get_protein_property"] = (
    "get protein property <property> for [<list of Proteins>] | <Protein>   USING (<parameter>=<value>)"
)


import glob
import json


def help_dict_create(
    name: str,  # Name of the comand - used for ...?
    command: str,  # Command structure, used for help, docs, training
    description: str,  # Description of the command, used for help, docs, training
    note: str = None,  # Additional note to the command, only used in help (eg. To learn more about runs, run `run ?`)
    url: str = None,  # Currently not used - URL to the documentation of the command?
    category: str = "Uncategorized",  # Category used to organize the commands in help & docs
    parent: str = None,  # Parent command, only relevant for follow-up commands like `result open`
):
    """Create a help dictionary"""
    return {
        "category": category,
        "name": name,
        "command": command,
        "description": description,
        "note": note,
        "url": url,
        "parent": parent,
    }


def get_services(reference, type) -> list:
    """pulls the list of available services for"""

    service_list = []
    service_files = glob.glob(reference + "/*.json")

    for file in service_files:

        with open(file, "r") as file_handle:
            try:
                jdoc = json.load(file_handle)

                service_list.append(jdoc)
            except Exception as e:
                print(e)
                print("invalid service json definition  " + file)
    return service_list


def service_grammar_add(statements: list, help: list, service_list: list):
    """defines the grammar available for managing molecules"""
    for schema in service_list:
        command = service_command_start[schema["service_type"]]
        valid_types = None
        valid_type = None
        first = True
        if len(list(schema["valid_types"])) > 1:

            for prop_type in list(schema["valid_types"]):
                if first == True:
                    valid_types = f"CaselessKeyword('{prop_type}')"
                    first = False
                else:
                    valid_types = valid_types + f" | CaselessKeyword('{prop_type}')"
            print(valid_types)
        else:
            valid_type = f'CaselessKeyword("{list(schema["valid_types"])[0]}")("type")'
            help_type = list(schema["valid_types"])[0]
        if valid_type is None:
            valid_type = f'Word("[")+OneOrMore({valid_types})("types")+Word("]")'
            help_type = "[ " + ", ".join(list(schema["valid_types"])) + " ] "
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

        stmt = eval(
            "Forward( "
            + command
            + "+"
            + valid_type
            + service_command_subject[schema["service_type"]]
            + expression
            + ")"
            + f'("{schema["service_name"]}@{schema["service_type"]}")'
        )

        statements.append(stmt)
        parameter_help = "\nParameters:\n"
        for parameter, description in dict(schema["parameters"]).items():
            print_description = ""
            for key, value in description.items():
                print_description = print_description + f"\n          {key} : {value}"

            parameter_help = parameter_help + f"      {parameter} : {print_description} \n"

        required_parameters = ""
        for i in schema["required_parameters"]:
            if required_parameters == "":
                required_parameters = "\n Required Parameters: \n"
            required_parameters = required_parameters + f"\n      - {i}"

        help.append(
            help_dict_create(
                name=schema["service_type"],
                category="Model",
                command=str(service_command_help[schema["service_type"]]).replace("<property>", help_type),
                description=parameter_help + required_parameters,
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
        elif "type" in inp_statement[clause][i]:
            type_str = "type"
        elif "allOf" in inp_statement[clause][i]:
            type_str = "type"
        parameter = "param_" + i
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
            else:
                expression = (
                    expression
                    + f" {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr_num('val'))('"
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
            else:
                expression = (
                    expression
                    + f" & {status}(Group( CaselessKeyword ('"
                    + i
                    + "') +Suppress('=')+key_val_expr_num('val'))('"
                    + parameter
                    + "'))"
                    + " "
                )
        ii = 1
    return expression


import os


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
    name = request_input.getName()

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
    for param in request_input.as_dict().keys():
        if str(param).startswith("param_"):
            actual_param = str(param)[6:]
            template["parameters"][actual_param] = request_input.as_dict()[param]["val"]

    return template


def convert(lst):
    """Used for for converting lists to strings."""
    return str(lst).translate("[],'")


from openad.helpers.output import output_error, output_table, output_success


def openad_model_requestor(cmd_pointer, parser):

    a_request = request_generate(parser)

    Endpoint = "http:// <insert here>"

    response = requests.post(Endpoint + "/service", json=a_request)
    try:
        result = pd.DataFrame(response.json())
    except:
        return "Error returned from service"

    return output_table(result)


###############################################################################
# getting code
if __name__ == "__main__":
    defs = get_services(
        "/Users/phildowney/services-build/Open-AD-Model-Service/openad-services/properties_service/openad_model_property_service/service_defs"
    )

    stmt_list = []
    help_list = []
    stmt_defs = Forward()
    service_grammar_add(stmt_list, help_list, defs)
    for x in help_list:
        print(x["command"])
        print(x["description"])
    for i in stmt_list:
        stmt_defs |= i
    request = "get crystal property absolute_energy for '/Users/phildowney/services-build/Open-AD-Model-Service/openad-model-inference/gt4sd_common/gt4sd_common/properties/tests/' using(algorithm_version=v0)"

    temp = stmt_defs.parseString(convert(request), parse_all=True)

    test = {
        "service_type": "get_crystal_property",
        "service_name": "get crystal absolute_energy",
        "parameters": {
            "property_type": [
                "absolute_energy",
            ],
            "algorithm_version": "v0",
            "subjects": [
                [
                    "1000041.cif",
                    "#------------------------------------------------------------------------------\n#$Date: 2015-01-27 21:58:39 +0200 (Tue, 27 Jan 2015) $\n#$Revision: 130149 $\n#$URL: svn://www.crystallography.net/cod/cif/1/00/00/1000041.cif $\n#------------------------------------------------------------------------------\n#\n# This file is available in the Crystallography Open Database (COD),\n# http://www.crystallography.net/\n#\n# All data on this site have been placed in the public domain by the\n# contributors.\n#\ndata_1000041\nloop_\n_publ_author_name\n'Abrahams, S C'\n'Bernstein, J L'\n_publ_section_title\n;\nAccuracy of an automatic diffractometer. measurement of the sodium\nchloride structure factors\n;\n_journal_coden_ASTM              ACCRA9\n_journal_name_full               'Acta Crystallographica (1,1948-23,1967)'\n_journal_page_first              926\n_journal_page_last               932\n_journal_paper_doi               10.1107/S0365110X65002244\n_journal_volume                  18\n_journal_year                    1965\n_chemical_formula_structural     'Na Cl'\n_chemical_formula_sum            'Cl Na'\n_chemical_name_systematic        'Sodium chloride'\n_space_group_IT_number           225\n_symmetry_cell_setting           cubic\n_symmetry_Int_Tables_number      225\n_symmetry_space_group_name_Hall  '-F 4 2 3'\n_symmetry_space_group_name_H-M   'F m -3 m'\n_cell_angle_alpha                90\n_cell_angle_beta                 90\n_cell_angle_gamma                90\n_cell_formula_units_Z            4\n_cell_length_a                   5.62\n_cell_length_b                   5.62\n_cell_length_c                   5.62\n_cell_volume                     177.5\n_refine_ls_R_factor_all          0.022\n_cod_database_code               1000041\nloop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\ny,z,x\nz,x,y\nx,z,y\ny,x,z\nz,y,x\nx,-y,-z\ny,-z,-x\nz,-x,-y\nx,-z,-y\ny,-x,-z\nz,-y,-x\n-x,y,-z\n-y,z,-x\n-z,x,-y\n-x,z,-y\n-y,x,-z\n-z,y,-x\n-x,-y,z\n-y,-z,x\n-z,-x,y\n-x,-z,y\n-y,-x,z\n-z,-y,x\n-x,-y,-z\n-y,-z,-x\n-z,-x,-y\n-x,-z,-y\n-y,-x,-z\n-z,-y,-x\n-x,y,z\n-y,z,x\n-z,x,y\n-x,z,y\n-y,x,z\n-z,y,x\nx,-y,z\ny,-z,x\nz,-x,y\nx,-z,y\ny,-x,z\nz,-y,x\nx,y,-z\ny,z,-x\nz,x,-y\nx,z,-y\ny,x,-z\nz,y,-x\nx,1/2+y,1/2+z\n1/2+x,y,1/2+z\n1/2+x,1/2+y,z\ny,1/2+z,1/2+x\n1/2+y,z,1/2+x\n1/2+y,1/2+z,x\nz,1/2+x,1/2+y\n1/2+z,x,1/2+y\n1/2+z,1/2+x,y\nx,1/2+z,1/2+y\n1/2+x,z,1/2+y\n1/2+x,1/2+z,y\ny,1/2+x,1/2+z\n1/2+y,x,1/2+z\n1/2+y,1/2+x,z\nz,1/2+y,1/2+x\n1/2+z,y,1/2+x\n1/2+z,1/2+y,x\nx,1/2-y,1/2-z\n1/2+x,-y,1/2-z\n1/2+x,1/2-y,-z\ny,1/2-z,1/2-x\n1/2+y,-z,1/2-x\n1/2+y,1/2-z,-x\nz,1/2-x,1/2-y\n1/2+z,-x,1/2-y\n1/2+z,1/2-x,-y\nx,1/2-z,1/2-y\n1/2+x,-z,1/2-y\n1/2+x,1/2-z,-y\ny,1/2-x,1/2-z\n1/2+y,-x,1/2-z\n1/2+y,1/2-x,-z\nz,1/2-y,1/2-x\n1/2+z,-y,1/2-x\n1/2+z,1/2-y,-x\n-x,1/2+y,1/2-z\n1/2-x,y,1/2-z\n1/2-x,1/2+y,-z\n-y,1/2+z,1/2-x\n1/2-y,z,1/2-x\n1/2-y,1/2+z,-x\n-z,1/2+x,1/2-y\n1/2-z,x,1/2-y\n1/2-z,1/2+x,-y\n-x,1/2+z,1/2-y\n1/2-x,z,1/2-y\n1/2-x,1/2+z,-y\n-y,1/2+x,1/2-z\n1/2-y,x,1/2-z\n1/2-y,1/2+x,-z\n-z,1/2+y,1/2-x\n1/2-z,y,1/2-x\n1/2-z,1/2+y,-x\n-x,1/2-y,1/2+z\n1/2-x,-y,1/2+z\n1/2-x,1/2-y,z\n-y,1/2-z,1/2+x\n1/2-y,-z,1/2+x\n1/2-y,1/2-z,x\n-z,1/2-x,1/2+y\n1/2-z,-x,1/2+y\n1/2-z,1/2-x,y\n-x,1/2-z,1/2+y\n1/2-x,-z,1/2+y\n1/2-x,1/2-z,y\n-y,1/2-x,1/2+z\n1/2-y,-x,1/2+z\n1/2-y,1/2-x,z\n-z,1/2-y,1/2+x\n1/2-z,-y,1/2+x\n1/2-z,1/2-y,x\n-x,1/2-y,1/2-z\n1/2-x,-y,1/2-z\n1/2-x,1/2-y,-z\n-y,1/2-z,1/2-x\n1/2-y,-z,1/2-x\n1/2-y,1/2-z,-x\n-z,1/2-x,1/2-y\n1/2-z,-x,1/2-y\n1/2-z,1/2-x,-y\n-x,1/2-z,1/2-y\n1/2-x,-z,1/2-y\n1/2-x,1/2-z,-y\n-y,1/2-x,1/2-z\n1/2-y,-x,1/2-z\n1/2-y,1/2-x,-z\n-z,1/2-y,1/2-x\n1/2-z,-y,1/2-x\n1/2-z,1/2-y,-x\n-x,1/2+y,1/2+z\n1/2-x,y,1/2+z\n1/2-x,1/2+y,z\n-y,1/2+z,1/2+x\n1/2-y,z,1/2+x\n1/2-y,1/2+z,x\n-z,1/2+x,1/2+y\n1/2-z,x,1/2+y\n1/2-z,1/2+x,y\n-x,1/2+z,1/2+y\n1/2-x,z,1/2+y\n1/2-x,1/2+z,y\n-y,1/2+x,1/2+z\n1/2-y,x,1/2+z\n1/2-y,1/2+x,z\n-z,1/2+y,1/2+x\n1/2-z,y,1/2+x\n1/2-z,1/2+y,x\nx,1/2-y,1/2+z\n1/2+x,-y,1/2+z\n1/2+x,1/2-y,z\ny,1/2-z,1/2+x\n1/2+y,-z,1/2+x\n1/2+y,1/2-z,x\nz,1/2-x,1/2+y\n1/2+z,-x,1/2+y\n1/2+z,1/2-x,y\nx,1/2-z,1/2+y\n1/2+x,-z,1/2+y\n1/2+x,1/2-z,y\ny,1/2-x,1/2+z\n1/2+y,-x,1/2+z\n1/2+y,1/2-x,z\nz,1/2-y,1/2+x\n1/2+z,-y,1/2+x\n1/2+z,1/2-y,x\nx,1/2+y,1/2-z\n1/2+x,y,1/2-z\n1/2+x,1/2+y,-z\ny,1/2+z,1/2-x\n1/2+y,z,1/2-x\n1/2+y,1/2+z,-x\nz,1/2+x,1/2-y\n1/2+z,x,1/2-y\n1/2+z,1/2+x,-y\nx,1/2+z,1/2-y\n1/2+x,z,1/2-y\n1/2+x,1/2+z,-y\ny,1/2+x,1/2-z\n1/2+y,x,1/2-z\n1/2+y,1/2+x,-z\nz,1/2+y,1/2-x\n1/2+z,y,1/2-x\n1/2+z,1/2+y,-x\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_symmetry_multiplicity\n_atom_site_Wyckoff_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n_atom_site_attached_hydrogens\n_atom_site_calc_flag\nNa1 Na1+ 4 a 0. 0. 0. 1. 0 d\nCl1 Cl1- 4 b 0.5 0.5 0.5 1. 0 d\nloop_\n_atom_type_symbol\n_atom_type_oxidation_number\nNa1+ 1.000\nCl1- -1.000\n",
                ]
            ],
            "subject_type": "CIF",
        },
        "api_key": "api-dthgwrhrthrtrth",
    }
    test2 = {
        "service_name": "get molecule activity_against_target",
        "service_type": "get_molecule_property",
        "parameters": {
            "property_type": [
                "activity_against_target",
            ],
            "subjects": [
                "CCO",
                "C1C2C(COC2O)C(O1)C3=CC4=C(C=C3)OCO4",
                "C1=C(N=NN1CC(C(=O)O)N)P(=O)(O)O",
                "C1=CC=C(C=C1)C2(C(=O)[N-]C(=O)N2)C3=CC=CC=C3.[Na+]",
                "CC1(CC(C2(C1C(OC=C2)OC3C(C(C(C(O3)CO)O)O)O)O)O)O ",
                "CC(C)(C)OC(=O)NC(CC1=CC=CC=C1)C(C[N+](C)(CC2=CC=CC=C2)NC(=O)C3=CC=CC=C3)O",
                "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            ],
            "target": "drd2",
        },
        "api_key": "api-dthgwrhrthrtrth",
    }
    scma2 = {
        "service_type": "get_molecule_property",
        "service_name": "get molecule activity_against_target",
        "subject": "smiles_string",
        "service_description": "",
        "valid_types": ["activity_against_target"],
        "type_description": {},
        "parameters": {
            "target": {"title": "Target", "description": "name of the target.", "example": "drd2", "type": "string"}
        },
        "required_parameters": ["target"],
        "category": "properties",
        "sub_category": "molecules",
        "wheel_package": "",
        "GPU": False,
        "persistent": False,
        "help": "",
    }
    scma = {
        "service_type": "get_crystal_property",
        "service_name": "get crystal absolute_energy",
        "subject": "cif_files",
        "service_description": "",
        "valid_types": ["absolute_energy"],
        "type_description": {},
        "parameters": {
            "algorithm_type": {"title": "Algorithm Type", "default": "prediction", "type": "string"},
            "domain": {"default": "crystals", "allOf": [{"$ref": "#/definitions/DomainSubmodule"}]},
            "algorithm_name": {"title": "Algorithm Name", "default": "cgcnn", "type": "string"},
            "algorithm_version": {
                "title": "Algorithm Version",
                "description": "Version of the algorithm",
                "example": "v0",
                "type": "string",
            },
            "algorithm_application": {"title": "Algorithm Application", "default": "AbsoluteEnergy", "type": "string"},
            "batch_size": {
                "title": "Batch Size",
                "description": "Prediction batch size",
                "default": 256,
                "type": "integer",
            },
            "workers": {
                "title": "Workers",
                "description": "Number of data loading workers",
                "default": 0,
                "type": "integer",
            },
        },
        "required_parameters": ["algorithm_version"],
        "category": "properties",
        "sub_category": "crystal",
        "wheel_package": "",
        "GPU": False,
        "persistent": True,
        "help": "",
    }
    # astate = service_grammar_add([], [dict(scma)])
    # print(astate.parse_string("get crystal property absolute_energy for 'CC0' using(algorithm_version=v0)").as_dict())
    # print(astate.runTests("get crystal property absolute_energy for 'CC0' using(algorithm_version=v0)"))
    # print(astate.runTests("get crystal property absolute_energy for 'CC0' "))

    # astate2 = service_grammar_add([], [dict(scma2)])

    """get molecule property activity_against_target for  [ CCO,C1C2C(COC2O)C(O1)C3=CC4=C(C=C3)OCO4,C1=C(N=NN1CC(C(=O)O)N)P(=O)(O)O,C1=CC=C(C=C1)C2(C(=O)[N-]C(=O)N2)C3=CC=CC=C3.[Na+],CC1(CC(C2(C1C(OC=C2)OC3C(C(C(C(O3)CO)O)O)O)O)O)O,CC(C)(C)OC(=O)NC(CC1=CC=CC=C1)C(C[N+](C)(CC2=CC=CC=C2)NC(=O)C3=CC=CC=C3)O,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O ] using (target=ddr2)"""

    """get molecule property activity_against_target for  [ CCO,C1C2C(COC2O)C(O1)C3=CC4=C(C=C3)OCO4,C1=C(N=NN1CC(C(=O)O)N)P(=O)(O)O,C1=CC=C(C=C1)C2(C(=O)[N-]C(=O)N2)C3=CC=CC=C3.[Na+],CC1(CC(C2(C1C(OC=C2)OC3C(C(C(C(O3)CO)O)O)O)O)O)O,CC(C)(C)OC(=O)NC(CC1=CC=CC=C1)C(C[N+](C)(CC2=CC=CC=C2)NC(=O)C3=CC=CC=C3)O,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O ] using (target=ddr2)"""

    """get molecule property activity_against_target for  [ CCO,C1C2C(COC2O)C(O1)C3=CC4=C(C=C3)OCO4,C1=C(N=NN1CC(C(=O)O)N)P(=O)(O)O,C1=CC=C(C=C1)C2(C(=O)[N-]C(=O)N2)C3=CC=CC=C3.[Na+],CC1(CC(C2(C1C(OC=C2)OC3C(C(C(C(O3)CO)O)O)O)O)O)O,CC(C)(C)OC(=O)NC(CC1=CC=CC=C1)C(C[N+](C)(CC2=CC=CC=C2)NC(=O)C3=CC=CC=C3)O,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O ] using (target=ddr2)"""

    """
    a_request = request_generate(astate_parse)

    import requests
    import pandas as pd

    Endpoint = "http://1"

    response = requests.post(Endpoint + "/service", json=a_request)
    print()
    print("-------------------------------------------------")
    print("Testing " + test["service_name"])
    print("---------------------------------------------")

    try:
        if response.json() != {}:
            print(pd.DataFrame(response.json()))
    except:
        print("fail")
    """
    # import openad_model_property_service.service_defs as new_prop_services
