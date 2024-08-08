# OpenAD Documentation

The folder contains the "single source of truth" for descriptions of OpenAD and its core concepts that are propagated across different aspects of the OpenAD help experience.

### How to regenerate documentation

    python3 docs/generate_docs.py

### What this does

When regenerating the docs, the text files in the `/source` folder are used to update:

-   The main [README.md](/) in this repository
-   The markdown files for [openad-docs] repository (they are copied over automatically if the openad-docs repo is found)
-   The intro text that is displayed when running the `intro` command (see [/openad/helpers/output_content.py])
-   The informational paragraphs that are displayed when running the `? workspace`, `? tookit`, `? context` and `? run` commands

### How it works

-   The `/source` folder contains a number of text files with descriptions for OpenAD and some of its core concepts.
-   The `/input` folder contains the markdown templates that are used to generate the markdown pages for the [openad-docs] repository
-   The `/output` folder contains all the files that are exported by the script:
    -   [csv/commands.csv](csv/commands.csv) - A CSV with all available commands that is not used anywhere but comes in handy to have a clean overview of available commands.
    -   [markdown/base-concepts.md](output/markdown/base-concepts.md)
    -   [markdown/commands.mdtext](output/markdown/commands.md)
    -   [markdown/index.md](output/markdown/index.md)\
    -   [markdown/installation.md](output/markdown/installation.md)

In addition to this, the script will also update the `llm_description.txt` for each individual plugin,

[openad-docs]: https://github.com/acceleratedscience/openad-docs
