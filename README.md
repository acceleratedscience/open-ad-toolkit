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

# Open Accelerated Discovery <!-- omit from toc -->

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/openad)](https://pypi.org/project/openad/)
[![PyPI version](https://img.shields.io/pypi/v/openad)](https://pypi.org/project/openad/)
[![License MIT](https://img.shields.io/github/license/acceleratedscience/open-ad-toolkit)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Docs](https://img.shields.io/badge/website-live-brightgreen)](https://acceleratedscience.github.io/openad-docs/)

<br>

<!-- description -->

OpenAD is an open-source framework for molecular and materials discovery developed by IBM Research.

The OpenAD client is accessible from a command line interface, Jupyter Notebook and an API. It provides unified access to a variety of tools and AI models for literature knowledge extraction, forward and retrosynthesis prediction, generative methods and property inference. You can train models on your own data as well as visualize and filter candidate molecules.

<!-- /description -->

[Installation](https://acceleratedscience.github.io/openad-docs/installation.html)
[Getting Started](https://acceleratedscience.github.io/openad-docs/getting-started.html)
[Models Service](https://acceleratedscience.github.io/openad-docs/models-service.html)
[Plugins](https://acceleratedscience.github.io/openad-docs/plugins.html)




<!-- [Documentation](https://acceleratedscience.github.io/openad-docs/) -->

<br>

### See What's New <!-- omit from toc -->

<details>
<summary>See what's new in OpenAD</summary>
<div markdown="block">

-   `%Openadd` has been added to the magic commands for commands that return data.
-   Upgraded SkyPilot to 0.6.0
-   Support for deploying in OpenShift AI/Open Data hub workbench or Podman/Docker image. [See the workbench repo](https://github.com/acceleratedscience/openad_workbench).
-   Support for application API
-   New property and dataset generation services.<br>We currently support the following model services:

    -   GT4SD Generation Services `git@github.com:acceleratedscience/generation_inference_service.git`
    -   GT4SD Property Services `git@github.com:acceleratedscience/property_inference_service.git`
    -   GT4SD MoleR Generation `git@github.com:acceleratedscience/moler_inference_service.git`
    -   GT4SD Molformer `git@github.com:acceleratedscience/molformer_inference_service.git`

    Pre-Requisite is that you have a AWS Account and can launch your own EC2 Instances Or someone else can launch them for you and you can catalog a Remote Service via URL.

    **Example:**

    -   Install a service

             catalog model service from 'git@github.com:acceleratedscience/property_inference_service.git' as prop

    -   Start the service

             model service up prop

    -   Wait until the service is ready

             model service status

    -   Once the service is ready, you can run the following commands to test:

        ```
        prop get molecule property [qed,esol] for [ C(C(C1C(=C(C(=O)O1)O)O)O)O ,[H-] ]
        ```

        ```
        prop get molecule property esol for C(C(C1C(=C(C(=O)O1)O)O)O)O
        ```

    -   Examples are supplied in the sample Notebooks.<br>See `init_examples` under the [Jupyter installation instructions](#setting-up-jupyter) below for more information.

    -   To shut down the service

            model service down prop

    -   Available commands for managing model services...

            model service status
            model service config <service_name>
            model catalog list
            uncatalog model service <service_name>
            catalog model service from (remote) '<path or github>' as <service_name>
            model service up <service_name> [no_gpu]
            model service local up <service_name>
            model service down <service_name>

</div>
</details>

<br>

## Table of Contents <!-- omit from toc -->

<!-- toc -->
- [Getting Started - CLI](#getting-started---cli)
- [Getting Started - Jupyter](#getting-started---jupyter)
  - [Setting up Jupyter](#setting-up-jupyter)
  - [Launching OpenAD in Jupyter](#launching-openad-in-jupyter)
<!-- tocstop -->

<br>


---

<br>

### Before You Start

-   OpenAD is available for Linux and MacOS
-   We support Windows 11 via WSL 2 (ubuntu 22.04) - see [Installing on Windows](#installing-on-windows)
-   When not installing into a virtual environment on MacOS, you may need to use `python3` and `pip3` instead of `python` and `pip` respectively
-   When updating to 0.4.0 or above, first remove all toolkits by runnning `list toolkits` and then `remove toolkit <toolkit_name>`.


<br>

## Quick Install

> **Note:** This will install OpenAD in your global space. If you wish to use a virtual environment (recommended), please see more [detailed instructions](#installation) below.

    pip install openad
    openad

Get started with Jupyter:

    init_magic
    init_examples
    jupyter lab ~/openad_notebooks/Table_of_Contents.ipynb

If you get an error when running `init_magic`, you may first need to setup the default iPython profile for magic commands.

    ipython profile create

<br>

---

<br>



# Getting Started - CLI

-   **Enter the virtual environment**

    > **Note:** If you just installed OpenAD, you probably already activated the virtual environment.

        source ~/ad-venv/bin/activate

-   **Enter the command shell**

        openad

    <!-- ![Landing](assets/screenshot-landing.png) -->
    <!-- <a href="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/screenshot-landing.png" target="_blank"><img src="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/screenshot-landing.png" /></a> -->

-   **Exit the command shell**<br>
    Hit `ctrl+c` or run:

        exit

-   **Run a single command from outside the command shell**

        openad <command>

-   **Exit the virtual environment**<br>

        deactivate

<br>

# Getting Started - Jupyter

## Setting up Jupyter

The following commands only need to be run once after installation:

1.  **Activate your virtual environment**

    > **Note:** If you just installed OpenAD, you probably already activated the virtual environment.

        source ~/ad-venv/bin/activate

1.  **Create an iPython kernel**<br>
    This ports your virtual environment to Jupyter.

        python -m ipykernel install --user --name=ad-venv

    > **Note:** To list your installed iPython kernels, you can run `jupyter kernelspec list`, and to remove the kernel you can run `jupyter kernelspec uninstall ad-venv`

1.  **Install the magic commands**<br>
    This enables OpenAD commands to be run within a Jupyter Notebook.

        init_magic

    <details>
    <summary><b>Alternative:</b> Manually add magic commands</summary>
    <div markdown="block">

    If you don't want to activate magic commands in all Notebooks, you can instead activate them for individual Notebooks.

    -   Run `init_examples`
    -   Copy the file `~/openad_notebooks/openad.ipynb` to the same directory as the Notebook you wish to activate.
    -   In your Notebook, run this inside a code cell: `!run openad.ipynb`

    </div>
    </details>

1.  **Install example Notebooks**<br>
    This installs our example Notebooks at `~/openad_notebooks`.

        init_examples

<br>

## Launching OpenAD in Jupyter

1.  **Open any Notebook**<br>
    The following command will open up the example Notebook:

        jupyter lab ~/openad_notebooks/Table_of_Contents.ipynb

1.  **Select the kernel**<br>
    Make sure to select the "ad-venv" iPython kernel. You can do this under _Kernel > Change Kernel_, or in the latest versions of Jupyter by clicking the kernel name in the top right hand corner. If you don't see your iPython kernel, make sure you followed the Jupyter Setup instructions listed above.

    <figure>
        <a href="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/jupyter-notebook-kernel.png" target="_blank"><img src="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/jupyter-notebook-kernel.png"></a>
    </figure>
    <figure>
        <a href="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/jupyter-lab-kernel.png" target="_blank"><img src="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/jupyter-lab-kernel.png"></a>
    </figure>

1.  **Magic Commands**<br>
    Magic commands let you run terminal commands from within Jupyter. They are invoked by the `%openad` prefix. All OpenAD CLI commands can be accessed like this. For example:<br>

        %openad list files