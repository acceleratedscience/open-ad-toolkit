<!--

DO NOT EDIT
-----------
This file is auto-generated with modified links for PyPI.
To update it, consult instructions:
https://github.com/acceleratedscience/open-ad-toolkit/tree/main/docs

-->



# OpenAD

**Open-source framework for molecular and materials discovery**

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/openad)](https://pypi.org/project/openad/)
[![PyPI version](https://img.shields.io/pypi/v/openad)](https://pypi.org/project/openad/)
[![License MIT](https://img.shields.io/github/license/acceleratedscience/open-ad-toolkit)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

<a href="assets/openad-cli.png" target="_blank"><img src="assets/openad-cli.png" width="500" height="328" /></a>

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

## Quick Install

> <div class='alert-icn-wrap' style='color:#8250df'><svg class="alert-icon" width="16" height="16" viewBox="0 0 16 16" fill="#8250df" xmlns="http://www.w3.org/2000/svg"><path d="M8 1C6.61553 1 5.26216 1.41054 4.11101 2.17971C2.95987 2.94888 2.06266 4.04213 1.53285 5.32122C1.00303 6.6003 0.86441 8.00776 1.13451 9.36563C1.4046 10.7235 2.07129 11.9708 3.05026 12.9497C4.02922 13.9287 5.2765 14.5954 6.63437 14.8655C7.99224 15.1356 9.3997 14.997 10.6788 14.4672C11.9579 13.9373 13.0511 13.0401 13.8203 11.889C14.5895 10.7378 15 9.38447 15 8C15 6.14348 14.2625 4.36301 12.9497 3.05025C11.637 1.7375 9.85652 1 8 1ZM8 14C6.81332 14 5.65328 13.6481 4.66658 12.9888C3.67989 12.3295 2.91085 11.3925 2.45673 10.2961C2.0026 9.19974 1.88378 7.99334 2.11529 6.82946C2.3468 5.66557 2.91825 4.59647 3.75736 3.75736C4.59648 2.91824 5.66558 2.3468 6.82946 2.11529C7.99335 1.88378 9.19975 2.0026 10.2961 2.45672C11.3925 2.91085 12.3295 3.67988 12.9888 4.66658C13.6481 5.65327 14 6.81331 14 8C14 9.5913 13.3679 11.1174 12.2426 12.2426C11.1174 13.3679 9.5913 14 8 14Z"/><path d="M8.5 4H7.5V9.5H8.5V4Z"/><path d="M8 11C7.85167 11 7.70666 11.044 7.58333 11.1264C7.45999 11.2088 7.36386 11.3259 7.30709 11.463C7.25033 11.6 7.23547 11.7508 7.26441 11.8963C7.29335 12.0418 7.36478 12.1754 7.46967 12.2803C7.57456 12.3852 7.7082 12.4567 7.85369 12.4856C7.99917 12.5145 8.14997 12.4997 8.28701 12.4429C8.42406 12.3861 8.54119 12.29 8.6236 12.1667C8.70602 12.0433 8.75 11.8983 8.75 11.75C8.75 11.5511 8.67098 11.3603 8.53033 11.2197C8.38968 11.079 8.19892 11 8 11Z"/></svg> IMPORTANT</div><span style='color: #8250df'>This will install OpenAD in your global space. If you wish to use a virtual environment (recommended), please refer to the [Installation](https://acceleratedscience.github.io/openad-docs/installation.html) page.</span>

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
-   New property and dataset generation services. See [OpenAD Models Service](https://acceleratedscience.github.io/openad-docs/models-service.html)
    -   GT4SD Generation Services
    -   GT4SD Property Services
    -   GT4SD MoleR Generation
    -   GT4SD Molformer

    Pre-Requisite is that you have a AWS Account and can launch your own EC2 Instances Or someone else can launch them for you and you can catalog a Remote Service via URL.

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