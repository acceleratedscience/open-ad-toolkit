<!--

Table of contents
-----------------
- Install the "Markdown All in One" plugin in VSCode
- The TOC should be automatically updated. If it's not:
  • Open the command palette (Press Cmd+Shift+P (macOS) or Ctrl+Shift+P (Windows))
  • Type "Markdown All in One: Update Table of Contents"


Screenshots
-----------
For screenshots to look good, they should be small and ideally
all the same size. The script below lets you open the URLs in
the right size. Just paste this into the browser console and
press enter.

To take the screenshots with browser UI included on Mac, press
cmd+shift+4 followed by the spacebar, then click the window.
For consistency, stick to Chrome, hide your bookmarks & extensions.

    urls = [
        'https://cps.foc-deepsearch.zurich.ibm.com',
        'https://rxn.app.accelerate.science',
        'https://sds.app.accelerate.science',

    ]
    for (var i=0; i< urls.length; i++) {
        window.open(urls[i], '_blank', 'width=1000,height=600');
    }

-->

# OpenAD

**Open-source framework for molecular and materials discovery**

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/openad)](https://pypi.org/project/openad/)
[![PyPI version](https://img.shields.io/pypi/v/openad)](https://pypi.org/project/openad/)
[![License MIT](https://img.shields.io/github/license/acceleratedscience/open-ad-toolkit)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Docs](https://img.shields.io/badge/website-live-brightgreen)](https://acceleratedscience.github.io/openad-docs/)

<a href="assets/openad-cli.png" target="_blank"><img src="assets/openad-cli.png" width="500" /></a>

<!-- description -->
Open Accelerated Discovery (aka OpenAD) is an open-source framework for molecular and materials discovery developed at IBM Research.

The OpenAD client is accessible from our command line interface, via Jupyter Notebook or our API. It provides unified access to a variety of tools and AI models for literature knowledge extraction, forward and retrosynthesis prediction, generative methods and property inference. OpenAD lets you train models on your own data, to then visualize and filter your candidate molecules.
<!-- /description -->

[Installation](https://acceleratedscience.github.io/openad-docs/installation.html)
&nbsp;&nbsp;&nbsp;
[Getting Started](https://acceleratedscience.github.io/openad-docs/getting-started.html)
&nbsp;&nbsp;&nbsp;
[Models Service](https://acceleratedscience.github.io/openad-docs/models-service.html)
&nbsp;&nbsp;&nbsp;
[Plugins](https://acceleratedscience.github.io/openad-docs/plugins.html)

<br><br>

## Quick Install

> [!IMPORTANT]  
> This will install OpenAD in your global space. If you wish to use a virtual environment (recommended), please refer to the [Installation](/readme/installation.md) page.

    pip install openad
    openad

Get started with Jupyter:

    init_magic
    init_examples
    jupyter lab ~/openad_notebooks/Table_of_Contents.ipynb

If you get an error when running `init_magic`, you may first need to setup the default iPython profile for magic commands.

    ipython profile create

<br><br>

## See What's New

-   `%Openadd` has been added to the magic commands for commands that return data.
-   Upgraded SkyPilot to 0.6.0
-   Support for deploying in OpenShift AI/Open Data hub workbench or Podman/Docker image. [See the workbench repo](https://github.com/acceleratedscience/openad_workbench).
-   Support for application API
-   New property and dataset generation services. See [OpenAD Models Service](/readme/models-service.md)
    -   GT4SD Generation Services
    -   GT4SD Property Services
    -   GT4SD MoleR Generation
    -   GT4SD Molformer

    Pre-Requisite is that you have a AWS Account and can launch your own EC2 Instances Or someone else can launch them for you and you can catalog a Remote Service via URL.

