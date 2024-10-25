# OpenAD Automated Documentation

The documentation generation script assures consistency by propagating the source descriptions of OpenAD and its base concepts across different aspects of the OpenAD help experience, while also automatically updating the documentation website with the README content of the OpenAD repository.

### How to regenerate documentation

    python3 docs/generate_docs.py

<br>

### What this does

1. Use the descriptions from the [`/source`](./source) folder to update:
   - The `README.md` with OpenAD description
   - The `README_plugins.md` with plugins description
   - The [openad_intro] text (/openad/helpers/output_content.py) that is displayed when running the `intro` command
   - The informational paragraphs that are displayed when running the `? workspace`, `? mws`, `? tookit`, `? context` and `? run` commands (see _output_content_ in [main.py](/openad/app/main.py))
   - The [Base Concepts](https://acceleratedscience.github.io/openad-docs/base-concepts.html) page on the documentation website
1. Duplicate all the README markdown pages, prepare them for ***Just the Docs*** consumption, then copy them over to the documentation repository if it's available
1. Generate two additional bespoke pages for the documentation website:
   - [Base Concepts](https://acceleratedscience.github.io/openad-docs/base-concepts.html)
   - [Commands](https://acceleratedscience.github.io/openad-docs/commands.html)
2. Update the [llm_description.txt](/openad/user_toolkits/DS4SD/llm_description.txt) for each individual plugin with the correct commands, used for LLM training.

<br>

### How it works

-   The `/source` folder is the "single source of truth" and contains a number of text files with descriptions for OpenAD and some of its core concepts.
-   The `/input` folder contains the markdown templates that are used to generate the markdown pages for the [openad-docs] repository. The header info in the input files defines the order of pages in the side navigation.
-   The `/output` folder contains all the files that are exported by the script:
    -   `/markdown`<br>The markdown files for the [openad-docs] repository:
        - [index.md](output/markdown/index.md)
        - [installation.md](output/markdown/installation.md)
        - [getting-started.md](output/markdown/getting-started.md)
        - [base-concepts.md](output/markdown/base-concepts.md)
        - [commands.md](output/markdown/commands.md)
        - [models-service.md](output/markdown/models-service.md)
        - [plugins.md](output/markdown/plugins.md)
        - [ai-assistant.md](output/markdown/ai-assistant.md)
        - [developers.md](output/markdown/developers.md)
    -   `/csv`
        -   [csv/commands.csv](output/csv/commands.csv) - A CSV file with all available commands that is not used anywhere but comes in handy to have a clean overview of available commands

<br>

### Link mapping

A key part of the translation from **GitHub** markdown to **Just the Docs** and **PyPI** markdown files, is to update internal links.

- **GitHub:** Links and relative and point to other README files.
- **[Documentation](https://acceleratedscience.github.io/openad-docs):** Links and relative and point to the other pages on the site.
- **[PyPI](https://pypi.org/project/openad/):** Links are absolute and point to the documentation website.
    
**Source: GitHub**

    [Link](README.md)
    [Link](README.md#quick-install)
    [Link](README_developers.md)
    [Link](README_developers.md#testing-a-branch)
    [Link]: README.md
    [Link]: README.md#quick-install
    [Link]: README_developers.md
    [Link]: README_developers.md#testing-a-branch
    
    [Link](/README.md)
    [Link](/README.md#quick-install)
    [Link](/README_developers.md)
    [Link](/README_developers.md#testing-a-branch)
    [Link]: /README.md
    [Link]: /README.md#quick-install
    [Link]: /README_developers.md
    [Link]: /README_developers.md#testing-a-branch

**Translation: Just the Docs**

    [Link](index.html)
    [Link](index.html#quick-install)
    [Link](developers.html)
    [Link](developers.html#testing-a-branch)
    [Link]: index.html
    [Link]: index.html#quick-install
    [Link]: developers.html
    [Link]: developers.html#testing-a-branch

    [Link](index.html)
    [Link](index.html#quick-install)
    [Link](developers.html)
    [Link](developers.html#testing-a-branch)
    [Link]: index.html
    [Link]: index.html#quick-install
    [Link]: developers.html
    [Link]: developers.html#testing-a-branch

**Translation: PyPI**

    [Link](https://acceleratedscience.github.io/openad-docs/)
    [Link](https://acceleratedscience.github.io/openad-docs/#quick-install)
    [Link](https://acceleratedscience.github.io/openad-docs/developers.html)
    [Link](https://acceleratedscience.github.io/openad-docs/developers.html#testing-a-branch)
    [Link]: https://acceleratedscience.github.io/openad-docs/
    [Link]: https://acceleratedscience.github.io/openad-docs/#quick-install
    [Link]: https://acceleratedscience.github.io/openad-docs/developers.html
    [Link]: https://acceleratedscience.github.io/openad-docs/developers.html#testing-a-branch

    [Link](https://acceleratedscience.github.io/openad-docs/)
    [Link](https://acceleratedscience.github.io/openad-docs/#quick-install)
    [Link](https://acceleratedscience.github.io/openad-docs/developers.html)
    [Link](https://acceleratedscience.github.io/openad-docs/developers.html#testing-a-branch)
    [Link]: https://acceleratedscience.github.io/openad-docs/
    [Link]: https://acceleratedscience.github.io/openad-docs/#quick-install
    [Link]: https://acceleratedscience.github.io/openad-docs/developers.html
    [Link]: https://acceleratedscience.github.io/openad-docs/developers.html#testing-a-branch

[openad-docs]: https://github.com/acceleratedscience/openad-docs
