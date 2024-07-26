import os
import json
import pickle
import base64
import string
import hashlib
from openad.app.global_var_lib import GLOBAL_SETTINGS
from time import sleep


def generate_smiles_hash(smiles_string):
    # Convert the SMILES string to bytes (required by hashlib)
    bytes_smiles = smiles_string.encode("utf-8")

    # Calculate the SHA-256 hash
    sha256_hash = hashlib.sha256(bytes_smiles).hexdigest()

    return sha256_hash


def create_valid_filename_from_smiles(smiles_string):
    # Generate the SHA-256 hash from the SMILES string
    hash_value = generate_smiles_hash(smiles_string)

    # Create the filename by adding the hash value and the .smiles extension
    filename = hash_value + ".smiles"

    return filename


def format_filename(s):
    """Take a string and return a valid filename constructed from the string.
    Uses a whitelist approach: any characters not present in valid_chars are
    removed. Also spaces are replaced with underscores.

    Note: this method may produce invalid filenames such as ``, `.` or `..`
    When I use this method I prepend a date string like '2009_01_15_19_46_32_'
    and append a file extension like '.txt', so I avoid the potential of using
    an invalid filename.
    """
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    filename = "".join(c for c in s if c in valid_chars)
    filename = filename.replace(" ", "_")  # I don't like spaces in filenames.
    filename = filename.replace("CON", "C0N")  # I don't like spaces in filenames.
    return filename


class rxn_helper:
    import pandas as pd

    _RXN_VARS_TEMPLATE = {"current_project": None, "current_project_id": None}

    def __init__(self) -> None:
        pass

    # Creates_the Cache Key
    def get_dataframe_from_file(self, cmd_pointer, filename):
        import pandas

        if not os.path.isfile(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + filename):
            raise BaseException("file does not exist")
        else:
            df = pandas.read_csv(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + filename)
        return df

    def get_column_as_list_from_dataframe(self, a_dataframe, column_name) -> list:
        import pandas

        if column_name in a_dataframe:
            return a_dataframe[column_name].tolist()
        else:
            return []

    def append_project(self, cmd_pointer, project_name, project_id):
        if not os.path.isdir(cmd_pointer.home_dir + "/RXN_Projects/"):
            os.mkdir(cmd_pointer.home_dir + "/RXN_Projects/")

        try:
            with open(cmd_pointer.home_dir + "/RXN_Projects/rxn_projects.pkl", "r") as handle:
                projects = json.loads(handle.read())
                handle.close()

        except Exception as e:
            projects = {}
        projects[project_name] = project_id

        try:
            with open(cmd_pointer.home_dir + "/RXN_Projects/rxn_projects.pkl", "w") as handle:
                json.dump(dict(projects), handle)
                handle.close()
        except Exception as e:
            # print(e)
            return False

        # print("Appending 3")
        import shutil
        from datetime import datetime

        shutil.copyfile(
            cmd_pointer.home_dir + "/RXN_Projects/rxn_projects.pkl",
            cmd_pointer.home_dir + "/RXN_Projects/rxn_projects_" + datetime.now().strftime("%Y-%m-%d_%H%M%S") + ".bup",
        )
        return True

    def get_all_projects(self, cmd_pointer):
        try:
            with open(cmd_pointer.home_dir + "/RXN_Projects/rxn_projects.pkl", "r") as handle:
                projects = json.loads(handle.read())
                handle.close()
        except:
            projects = {}

        return projects

    def get_project_id(self, cmd_pointer, project_name):
        try:
            with open(cmd_pointer.home_dir + "/RXN_Projects/rxn_projects.pkl", "r") as handle:
                projects = json.loads(handle.read())
                handle.close()
                result = projects[project_name]
        except Exception as e:
            # print(e)

            result = False

        return result

    def gen_cache_key(self, chem_list: list) -> str:
        chem_list.sort()
        key = ""
        delimiter = ""
        for i in chem_list:
            key = key + delimiter + i
            delimiter = "_"

        key = key

        return key

    # saves cache_handle
    def save_to_results_cache(self, cmd_pointer, chem_list: str, payload, call_type: str) -> bool:
        try:
            self.create_cache(cmd_pointer)
            key = self.gen_cache_key(chem_list)

            payload_file = cmd_pointer.toolkit_dir + "/RXN/cache/" + key + "_" + call_type + ".result"

            # payload_file = json.loads(payload_file)
            with open(payload_file, "wb") as handle:
                settings = pickle.dump(dict(payload), handle)
        except BaseException as e:
            # print('failed to save result')
            # print(e)
            return False
        return True

    # creates the Cache directory for results
    def create_cache(self, cmd_pointer):
        if not os.path.isdir(cmd_pointer.toolkit_dir + "/RXN/cache"):
            os.mkdir(cmd_pointer.toolkit_dir + "/RXN/cache")

    # retrieves a matching result from the cache
    def retrieve_cache(self, cmd_pointer, chem_list: str, call_type: str):
        import sys

        try:
            self.create_cache(cmd_pointer)
            key = self.gen_cache_key(chem_list)

            payload_file = cmd_pointer.toolkit_dir + "/RXN/cache/" + key + "_" + call_type + ".result"

            with open(payload_file, "rb") as handle:
                result = pickle.loads(handle.read())
            return result
        except BaseException as e:
            return False

    # checks to see if a smiles string is valid
    def valid_smiles(self, input_molecule) -> bool:
        from rdkit import Chem

        m = Chem.MolFromSmiles(input_molecule, sanitize=False)
        if m is None:
            return False
        else:
            try:
                Chem.SanitizeMol(m)
            except:
                False
        return True

    # sets the current project
    def set_current_project(self, cmd_pointer, project_name: str) -> bool:
        # projects = self.get_all_projects(cmd_pointer)
        # project_id=1
        # project_description=1
        # result = self.validate_project(cmd_pointer,project_name)
        result = self.get_project_id(cmd_pointer, project_name)
        # print("setting project")
        # print(result)
        if result != False:
            rxn_position = cmd_pointer.login_settings["toolkits"].index("RXN") - 1
            # cmd_pointer.login_settings['session_vars'][rxn_position]['env_vars']= self._RXN_VARS_TEMPLATE
            cmd_pointer.login_settings["session_vars"][rxn_position]["current_project"] = project_name

            cmd_pointer.login_settings["session_vars"][rxn_position]["current_project_id"] = result
            # cmd_pointer.login_settings['session_vars'][rxn_position]['current_project_description']=result['description']
            rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][
                cmd_pointer.login_settings["toolkits"].index("RXN")
            ]
            # print('prj_id')
            # print(result)
            rxn4chemistry_wrapper.set_project(result)
            return result
        else:
            return False

    ### only for function checks not for login.py
    def sync_up_workspace_name(self, cmd_pointer, reset=False):
        # print("syncing)")
        name, id = self.get_current_project(cmd_pointer)
        # print("current_project: "+str(name)+" "+str(id))
        if name == cmd_pointer.settings["workspace"] and reset != True:
            return True

        try:
            result = self.set_current_project(cmd_pointer, cmd_pointer.settings["workspace"].upper())
        except BaseException as e:
            raise BaseException(
                "Unable to set current project due to API issue, check server connections Try Again: " + str(e)
            )

        if result == False:
            retries = 0
            while retries < 5 and result == False:
                if retries > 1:
                    sleep(3)
                retries = retries + 1

                # sys.stdout = open(os.devnull, "w")
                # sys.stderr = open(os.devnull, "w")

                try:
                    rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][
                        cmd_pointer.login_settings["toolkits"].index("RXN")
                    ]
                    sleep(3)

                    x = rxn4chemistry_wrapper.create_project(cmd_pointer.settings["workspace"])

                    if len(x) == 0:
                        continue
                    else:
                        self.append_project(
                            cmd_pointer, cmd_pointer.settings["workspace"].upper(), x["response"]["payload"]["id"]
                        )
                    # print('here')
                except BaseException as e:
                    print(e)
                    raise BaseException("Unable to create project :" + str(e))
                try:
                    sleep(3)
                    result = self.set_current_project(cmd_pointer, cmd_pointer.settings["workspace"])
                    # print(cmd_pointer.settings['workspace'])
                    # print(result)
                except BaseException as e:
                    # sys.stdout = sys.__stdout__
                    # sys.stderr = sys.__stderr__
                    raise BaseException(
                        "Unable to set current project due to API issue, check server connections Try Again: " + str(e)
                    )

            # sys.stdout = sys.__stdout__
            # sys.stderr = sys.__stderr__
            if result == False:
                if GLOBAL_SETTINGS["display"] == "notebook":
                    from IPython.display import display, Markdown

                    display(Markdown("Failed to setup RXN project for this Workspace "))
                else:
                    print("Failed to setup RXN project for this Workspace ")

            else:
                if GLOBAL_SETTINGS["display"] == "notebook":
                    from IPython.display import display, Markdown

                    display(Markdown("A new RXN Project has been setup for this Workspace"))
                else:
                    print("A new RXN Project has been setup for this Workspace ")
        else:
            if GLOBAL_SETTINGS["display"] == "notebook":
                from IPython.display import display, Markdown

                display(
                    Markdown("RXN Project has been set to " + cmd_pointer.settings["workspace"] + " for this Workspace")
                )

        return True

    def validate_project(self, cmd_pointer, project_name):
        # print("validate_project")
        # projects = self.get_all_projects(cmd_pointer)
        projects = self.get_all_projects(cmd_pointer)
        # print(projects)

        if project_name in projects["name"].values:
            result = projects[projects.name == project_name]
            # print(result)
            return {"name": result["name"], "id": result["id"]}
            # return {'name':result['name'][:1].item(),'id':result['id'][:1].item(),'description':result['description'][:1].item(),'attempts':result['attempts'][:1].item()}
        else:
            return False

    def get_current_project(self, cmd_pointer):
        # print("get_current_project")
        rxn_position = cmd_pointer.login_settings["toolkits"].index("RXN") - 1
        try:
            return (
                cmd_pointer.login_settings["session_vars"][rxn_position]["current_project"],
                cmd_pointer.login_settings["session_vars"][rxn_position]["current_project_id"],
            )
        except:
            return None, None

    # Note get
    def get_all_projects_old(self, cmd_pointer) -> pd.DataFrame:
        api_key = cmd_pointer.login_settings["toolkits_api"][cmd_pointer.login_settings["toolkits"].index("RXN")]

        rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][
            cmd_pointer.login_settings["toolkits"].index("RXN")
        ]
        # Prepare the data query
        source_list = []
        result = False
        retries = 0

        while result == False:
            try:
                x = rxn4chemistry_wrapper.list_all_projects(size=100)["response"]["payload"]["content"]
                # print(x)
                result = True
            except BaseException as e:
                retries = retries + 1
                sleep(2)
                if retries > 10:
                    raise BaseException("Unable to retrieve valid list of projects:" + str(e))

        df = self.pd.DataFrame(x)

        df = df[["name", "description", "id", "attempts"]]
        return df

    def output(self, text, cmd_pointer=None):
        if GLOBAL_SETTINGS["display"] == "notebook":
            from IPython.display import Markdown, display
        import re

        tags = {
            "h1": "\x1b[0m",  # Primary headers: default with yellow line below
            "h2": "\x1b[33m",  # Secondary headers: yellow
            "cmd": "\x1b[36m",  # Commands: cyan
            "error": "\x1b[31m",  # Errors: red
            "warning": "\x1b[33m",  # Warnings: yellow
            "success": "\x1b[32m",  # Success: green
            "link": "\x1b[4;94m",  # Links: undertline intense blue
            # 'edit': '\x1b[0;47;30m',  # Edit node: black on white
            "edit": "\x1b[7m",  # Edit node: reverse
            "editable": "\x1b[100m",  # Editable text: on_bright_black
            # Styles
            "reset": "\x1b[0m",
            "bold": "\x1b[1m",
            "soft": "\x1b[2m",
            "italic": "\x1b[3m",
            "underline": "\x1b[4m",
            "blink": "\x1b[5m",
            "reverse": "\x1b[7m",
            "hidden": "\x1b[8m",
            "strikethrough": "\x1b[9m",
            # Foreground colors - Regular
            "black": "\x1b[30m",
            "red": "\x1b[31m",
            "green": "\x1b[32m",
            "yellow": "\x1b[33m",
            "blue": "\x1b[34m",
            "magenta": "\x1b[35m",
            "cyan": "\x1b[36m",
            "white": "\x1b[37m",
            # Foreground colors - Intense (non-standard)
            "bright_black": "\x1b[90m",
            "bright_red": "\x1b[91m",
            "bright_green": "\x1b[92m",
            "bright_yellow": "\x1b[93m",
            "bright_blue": "\x1b[94m",
            "bright_magenta": "\x1b[95m",
            "bright_cyan": "\x1b[96m",
            "bright_white": "\x1b[97m",
            # Background colors - Regular
            "on_black": "\x1b[0;40m",
            "on_red": "\x1b[0;41m",
            "on_green": "\x1b[0;42m",
            "on_yellow": "\x1b[0;43m",
            "on_blue": "\x1b[0;44m",
            "on_magenta": "\x1b[0;45m",
            "on_cyan": "\x1b[0;46m",
            "on_white": "\x1b[0;47m",
            # Background colors - Intense (non-standard)
            "on_bright_black": "\x1b[100m",
            "on_bright_red": "\x1b[101m",
            "on_bright_green": "\x1b[102m",
            "on_bright_yellow": "\x1b[103m",
            "on_bright_blue": "\x1b[104m",
            "on_bright_magenta": "\x1b[105m",
            "on_bright_cyan": "\x1b[106m",
            "on_bright_white": "\x1b[107m",
            # Unused but useful as a reference
            # \x1b[7m - Reversed
        }

        def parse_tags(text: str):
            """
            Parse xml tags and return styled output:
            <red>lorem ipsum</red>
            """
            pattern = rf"(.*?)<({'|'.join(list(tags))})>(.*?)</\2>"
            return re.sub(pattern, lambda match: _replace(match, pattern), text)

        def _replace(match: object, pattern, parent_tag="reset"):
            """Replace regex matches with appropriate styling."""
            text_before = match.group(1)
            tag = match.group(2)
            inner_text = match.group(3)
            ansi_code_open = tags[tag]
            ansi_code_close = tags[parent_tag]
            ansi_code_reset = tags["reset"]

            # Replace any nested tags.
            if re.findall(pattern, inner_text):
                inner_text = re.sub(pattern, lambda match: _replace(match, pattern, tag), inner_text)

            return f"{text_before}{ansi_code_open}{inner_text}{ansi_code_reset}{ansi_code_close}"

        def strip_tags(text: str):
            """Recursively remove all XML tags."""

            # Strip any nested tags.
            def _strip(match: object, pattern):
                inner_text = match.group(2)
                if re.findall(pattern, inner_text):
                    inner_text = re.sub(pattern, lambda match: _strip(match, pattern), inner_text)
                return inner_text

            # Strip tags.
            text = text.replace("\n", "---LINEBREAK2---")
            pattern = rf"<({'|'.join(list(tags))})>(.*?)</\1>"
            text = re.sub(pattern, lambda match: _strip(match, pattern), text)
            return text.replace("---LINEBREAK2---", "\n")

        def tags_to_markdown(text: str):
            if text is None:
                return ""

            # Replace leading spaces with non-breaking spaces.
            text = re.sub(r"^( *)", "", text, flags=re.MULTILINE)

            # Replace line breaks so all text is parsed on one line.
            # Because html breaks (<br>) don't play well with headings,
            # and end of line characters don't play well with `code`
            # blocks, we have to do some trickery here.
            text = re.sub(
                r"(</h[123]>)(\n+)", lambda match: match.group(1) + len(match.group(2)) * "---LINEBREAKSOFT---", text
            )
            text = re.sub(
                r"(\n+)(<h[123]>)", lambda match: len(match.group(1)) * "---LINEBREAKSOFT---" + match.group(2), text
            )
            text = text.replace("\n", "---LINEBREAK3---")

            # Replace tags
            # We only replace <soft> and <underline> tags
            # when they don't appear inside <cmd> tags.
            text = re.sub(r"(?<!\<cmd\>)<soft>([^<]*)</soft>(?!\</cmd\>)", r'<span style="color: #ccc">\1</span>', text)
            text = re.sub(
                r"(?<!\<cmd\>)<underline>([^<]*)<\/underline>(?!\</cmd\>)",
                r'<span style="text-decoration: underline">\1</span>',
                text,
            )
            text = re.sub(r"<h1>(.*?)<\/h1>", r"## \1", text)
            text = re.sub(r"<h2>(.*?)<\/h2>", r"### \1", text)
            text = re.sub(r"<link>(.*?)<\/link>", r'<a target="_blank" href="\1">\1</a>', text)
            text = re.sub(r"<bold>(.*?)<\/bold>", r"**\1**", text)
            text = re.sub(r"<cmd>(.*?)<\/cmd>", r"`\1`", text)
            text = re.sub(r"<on_red>(.*?)<\/on_red>", r'<span style="background: #d00; color: #fff">\1</span>', text)
            text = re.sub(
                r"<on_green>(.*?)<\/on_green>", r'<span style="background: #0d0; color: #fff">\1</span>', text
            )
            text = re.sub(r"<success>(.*?)<\/success>", r'<span style=" color: #008000">\1</span>', text)
            text = re.sub(r"<fail>(.*?)<\/fail>", r'<span style=" color: #ff0000">\1</span>', text)

            # Escape quotes.
            text = text.replace("'", "'")

            # Replace all other tags
            text = strip_tags(text)

            # Restore line breaks.
            text = text.replace("---LINEBREAKSOFT---", "\n")
            text = text.replace("---LINEBREAK3---", "<br>")

            return text

        if GLOBAL_SETTINGS["display"] == "api":
            # API
            return strip_tags(text)
        elif GLOBAL_SETTINGS["display"] == "notebook":
            # Jupyter
            return Markdown(tags_to_markdown(text))
        else:
            # CLI
            print(parse_tags(str(text)))
