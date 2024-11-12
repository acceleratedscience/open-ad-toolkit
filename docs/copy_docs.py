"""
Function that copies over files to the openad-docs repo.
"""

import os
import shutil
from openad.helpers.output import output_error, output_warning, output_text, output_success

# Define source and destination directories
SOUCRE_DIR = "docs/output/markdown"
DEST_DIR = "../openad-docs"
DIR_STRUCTURE = [
    "/my-repos       <-- Parent folder",
    "  /openad-docs  <-- Destination",
    "  /openad       <-- Current repository",
    "    /docs",
    "      /output   <-- Output files to be copied",
]
FLAG_SUCCESS = f"<on_green> COPIED </on_green>"
FLAG_ERROR = f"<on_red> FAILED </on_red>"


def copy_docs(files, source_dir=SOUCRE_DIR, dest_dir=DEST_DIR):
    # Title
    output_text("\n\n\n<h1>Copying readme files to the <yellow>openad-docs</yellow> repository</h1>")

    # Source folder not found
    if not os.path.exists(source_dir):
        output_warning(
            [
                FLAG_ERROR
                + f" The source folder <error>{source_dir}</error> was not found, no files have been copied.\n"
                "<reset>To have the output files copied automatically, make sure to place the openad-docs repo adjacent to this repo:</reset>",
                "\n".join(DIR_STRUCTURE),
            ],
        )
        return

    # Destination folder not found
    if not os.path.exists(dest_dir):
        output_warning(
            [
                FLAG_ERROR + " The <error>openad-docs</error> repository was not found, no files have been copied.\n"
                "<reset>To have the output files copied automatically, make sure to place the openad-docs repo adjacent to this repo:</reset>",
                "\n".join(DIR_STRUCTURE),
            ],
        )
        return

    for [i, file_name] in enumerate(files):
        # Construct full file path
        source_file = os.path.join(source_dir, file_name)
        dest_file = os.path.join(dest_dir, file_name)

        if not os.path.exists(source_file):
            output_text(
                FLAG_ERROR + f" <red>{file_name}</red> could not be found at <yellow>{source_file}</yellow>.",
                pad_btm=1,
            )
            return

        # Copy the file
        shutil.copy2(source_file, dest_file)
        output_text(
            FLAG_SUCCESS + f" <yellow>{file_name}</yellow> to the <reset>openad-docs</reset> repository.",
            pad_btm=1,
        )
