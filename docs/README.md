# OpenAD Documentation

The documentation generation script assures consistency by propagating the source descriptions of OpenAD and its base concepts across different aspects of the OpenAD help experience, while also automatically updating the documentation website with the README content of the OpenAD repository.

### How to regenerate documentation

    python3 docs/generate_docs.py

### What this does

1. Use the descriptions from the [`/source`](/source) folder to update:
   - The `README.md`
   - The [intro text](/openad/helpers/output_content.py) that is displayed when running the `intro` command
   - The informational paragraphs that are displayed when running the `? workspace`, `? mws`, `? tookit`, `? context` and `? run` commands (see _output_content_ in [main.py](/openad/app/main.py))
   - The [Base Concepts](https://acceleratedscience.github.io/openad-docs/base-concepts.html) page on the documentation website
2. Duplicate all the README markdown pages, prepare them for just-the-docs consumption, then copy them over to the documentation repository if it's available

### How it works

-   The `/source` folder is the "single source of truth" and contains a number of text files with descriptions for OpenAD and some of its core concepts.
-   The `/input` folder contains the markdown templates that are used to generate the markdown pages for the [openad-docs] repository
-   The `/output` folder contains all the files that are exported by the script:
    -   **/markdown**<br>The markdown files for the [openad-docs] repository:
        -   [markdown/base-concepts.md](output/markdown/base-concepts.md)
        -   [markdown/commands.md](output/markdown/commands.md)
        -   [markdown/index.md](output/markdown/index.md)
        -   [markdown/installation.md](output/markdown/installation.md)
    -   **/csv**
        -   [csv/commands.csv](output/csv/commands.csv) - A CSV file with all available commands that is not used anywhere but comes in handy to have a clean overview of available commands

In addition to this, the script will also update the `llm_description.txt` for each individual plugin with the correct commands, used for LLM training.

[openad-docs]: https://github.com/acceleratedscience/openad-docs
