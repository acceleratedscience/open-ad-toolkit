<!--

DO NOT EDIT
-----------
This file is auto-generated with modified links for PyPI.
To update it, consult instructions:
https://github.com/acceleratedscience/open-ad-toolkit/tree/main/docs

-->



<img src="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/molecule-header.jpg" width="830" height="300" alt="OpenAD" />

# OpenAD

**Open-source framework for molecular and materials discovery**

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/openad)](https://pypi.org/project/openad/)
[![PyPI version](https://img.shields.io/pypi/v/openad)](https://pypi.org/project/openad/)
[![License MIT](https://img.shields.io/github/license/acceleratedscience/open-ad-toolkit)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Open Accelerated Discovery (aka OpenAD) is an open-source framework for molecular and materials discovery developed at IBM Research.

The OpenAD client is accessible from our command line interface, via Jupyter Notebook or our API. It provides unified access to a variety of tools and AI models for literature knowledge extraction, forward and retrosynthesis prediction, generative methods and property inference. OpenAD lets you train models on your own data, to then visualize and filter your candidate molecules.

**[INSTALLATION]**
&nbsp;&nbsp;&nbsp;
**[GETTING STARTED]**
&nbsp;&nbsp;&nbsp;
**[COMMANDS]**
&nbsp;&nbsp;&nbsp;
**[MODELS SERVICE]**

<sub>

**[PLUGINS]**
&nbsp;&nbsp;&nbsp;
**[AI ASSISTANT]**
&nbsp;&nbsp;&nbsp;
**[DEVELOPERS]**

</sub>

<br><br>

### Useful links

- **Homepage**<br>
  [accelerate.science/projects/openad](https://accelerate.science/projects/openad)
- **GitHub**<br>
  [github.com/acceleratedscience/open-ad-toolkit](https://github.com/acceleratedscience/open-ad-toolkit)
- **Pip install**<br>
  [pypi.org/project/openad](https://pypi.org/project/openad)
- **Documentation**<br>
  [acceleratedscience.github.io/openad-docs](https://acceleratedscience.github.io/openad-docs)

<br><br>

<div align="center"><img src="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/openad-ui.png" width="100%" alt="OpenAD interface" /></div>

<br><br>

## Quick Install

> [!IMPORTANT]
> This will install OpenAD in your global space. If you wish to use a virtual environment (recommended), please refer to the [Installation] page.

    pip install openad
    openad

Get started with Jupyter Notebook examples:

    init_magic
    init_examples
    jupyter lab ~/openad_notebooks/Table_of_Contents.ipynb

If you get an error when running `init_magic`, you may first need to setup the default iPython profile for magic commands.

    ipython profile create

<br><br>

## What's New

-   `%Openadd` has been added to the magic commands for commands that return data.
-   Upgraded SkyPilot to 0.6.0
-   Support for deploying in OpenShift AI/Open Data hub workbench or Podman/Docker image. [See the workbench repo](https://github.com/acceleratedscience/openad_workbench).
-   Support for application API
-   New property and dataset generation services. See OpenAD [Models Service]
    -   GT4SD Generation Services
    -   GT4SD Property Services
    -   GT4SD MoleR Generation
    -   GT4SD Molformer

    You'll need an AWS account and the ability to launch your own EC2 instances and catalog a remote service via URL. Instructions provided.

<br><br>

## Feedback

- **Feature requests, feedback and questions:**<br>
  [phil.downey1@ibm.com](mailto:phil.downey1@ibm.com)
- **Bug reports:**<br>
  [Create issue on GitHub](https://github.com/acceleratedscience/open-ad-toolkit/issues)

<br><br>

## Learn

- [OpenAD blog](https://blog.accelerate.science/)
- **Demo Notebooks:** See `init_examples` under installation

<br><br>

## Contribute
Stay tuned for detailed instructions on how to build your own OpenAD plugins.<br>
Check the [Developers] section for guidance with other contributions.


[INSTALLATION]: https://acceleratedscience.github.io/openad-docs/installation.html
[GETTING STARTED]: https://acceleratedscience.github.io/openad-docs/getting-started.html
[COMMANDS]: https://acceleratedscience.github.io/openad-docs/commands.html
[MODELS SERVICE]: https://acceleratedscience.github.io/openad-docs/models-service.html
[PLUGINS]: https://acceleratedscience.github.io/openad-docs/plugins.html
[AI ASSISTANT]: https://acceleratedscience.github.io/openad-docs/ai-assistant.html
[DEVELOPERS]: https://acceleratedscience.github.io/openad-docs/developers.html
