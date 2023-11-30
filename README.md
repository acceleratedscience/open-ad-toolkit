# OpenAD
[![PyPI version](https://img.shields.io/pypi/v/openad)](https://pypi.org/project/openad/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/openad)](https://pypi.org/project/openad/)
[![License MIT](https://img.shields.io/github/license/acceleratedscience/open-ad-toolkit)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Docs](https://img.shields.io/badge/website-live-brightgreen)](https://github.com/acceleratedscience/open-ad-toolkit)

**Open Accelerated Discovery Client**<br>
[Project homepage](https://acceleratedscience.github.io/open-ad-toolkit/)

---

<br>

## Table of Contents

- [Installation](#installation)
- [Getting Started - CLI](#getting-started---cli)
- [Getting Started - Jupyter ](#getting-started---cli)
- [Getting Access to RXN, DS4SD and Tell Me Functionality](#getting-access-to-rxn-ds4sd-and-tell-me-functionality)
- [Installation for Development](#installation-for-development)
- [Testing a Branch](#testing-a-branch)

<br>

## Notes

-   Only available for Linux (including Windows 11 WSL) and MacOS
-   Currently only the OpenAI API is available for the _Tell Me_ Function (WatsonX coming soon)
-   If you're on Mac and not installing into a virtual environment, you may need use `pip3` and `python3` instead of `pip` and `python` respectively.<br>

<br>

## Installation

> **Note:** If you're installing a development version of OpenAD, please jump to [Installation for Development](#installation-for-development)

1.  **Step 0: Before you start**<br>
    Ensure you're running Python 3.10.10 or 3.11.<br>
    We do not have a requirement on how to install Python on your local machine;however, we recommend Pyenv.

    ```shell
    git clone https://github.com/pyenv/pyenv.git ~/.pyenv
    pyenv install 3.10
    ```


1.  **Step 1: Set up virtual environment** (optional)<br>

        python -m venv ~/ad-venv
        source ~/ad-venv/bin/activate

    > **Note:** To exit the virtual environment, you can run `deactivate`
    > Also Openad Supports a Poetry install

1.  **Step 2: Installation**<br>
       To Install into a global or Python virtual environment<br>
       
               pip install openad

       To install into an Poetry Environment

                poetry add openad
    

1.  **Launch**<br>
    To enter the command shell, simply enter `openad` from the command line.

    > _**Note:** To see available commands, run `?`_

<br>

## Getting Started - CLI

-   **Entering the Shell Environment**<br>
    Run from any directory:

        openad

    ![screenshot-landing](https://github.com/acceleratedscience/open-ad-toolkit/assets/30096303/8c78414b-89fe-4662-8c23-5c5443e79084)

-   **Exiting the Shell Environment**<br>
    Hit `ctrl+c` or run:

        exit

-   **Installing Toolkits**<br>
    You can install the `DS4SD`, `GT4SD`, `ST4SD` and `RXN` toolkits, however please note that at this time, only `DS4SD` and `RXN` support experimental functionality while the others are meant as placeholders.

        add toolkit ds4sd
        add toolkit rxn

-   **Running Bash Commands**<br>
    To run any command in bash mode, prepend it with `openad` and make sure to prepend any quotes with `\`.

        openad show molecules using file \'base_molecules.sdf\'

<br>

## Getting Started - Jupyter

### Jupyter Setup

If you plan to use this application inside Jupyter Notebook of JupyterLab, you should set it up as follows:

1.  **Activate your virtual environment**<br>
    This should already be active, per [Installation](#installation) instructions on top.

        source ~/ad-venv/bin/activate

2.  **Create an iPython kernel**<br>
    This is used to run Notebook commands inside the virtual environment:

        python -m ipykernel install --user --name=ad-kernel

    if you are also using Poetry we suggest naming your kernel  `poetry-ad-kernel`

    > **Note:** To list your installed iPython kernels, you can run `jupyter kernelspec list`<br>
    > To remove the kernel you can run `jupyter kernelspec uninstall ad-kernel`

4.  **Initiate the magic commands.**<br>
    This copies the magic commands into the iPython startup directory for your created profile:

        init_magic

    > **Note:** You can also:<br>
    > • Initiate in a single notebook only: `init_magic .` followed by `run openad.py`<br>
    > • Initiate within another iPython profile: `init_magic <profile_name>`<br>
    > • Initiate without installing anything: `run openad.ipynb` after installing examples (see next bullet)<br>
    > --> This executes the file `~/openad_notebooks/openad.ipynb` and needs to be run every time after restarting the kernel

5.  **Install the Notebook examples**<br>
    This creates the `~/openad_notebooks` folder with a number of demonstration notbeooks.<br>
    Start with `Table_of_Contents.ipynb`.

         init_examples

<br>

### Jupyter Launch

-   After installing the Notebook examples, open the table of contents to get an introduction and be taken through step by step how to use the tool.

        jupyter lab ~/openad_notebooks/Table_of_Contents.ipynb

    > **NOTE:** By launching Jupyter this way it will automatically launch the trial notebooks.

-   Make sure to select your newly created "ad-kernel" iPython kernel. You can do this under _Kernel > Change Kernel_, or in the latest versions of Jupyter by clicking the kernel name in the top right hand corner. If you don't see your iPython kernel, make sure you followed the Jupyter Setup instructions listed above.

<br>
<figure>
    <figcaption align="center" style="font-size:0.9em;opacity:.6;margin-bottom:-15px"><i>Jupyter Notebook</i></figcaption>
    <img width="1112" alt="jupyter-notebook-kernel" src="https://github.com/acceleratedscience/open-ad-toolkit/assets/30096303/50dee4a9-80d9-4ddb-92f0-321d6c8328e5">

</figure>

<figure>
    <figcaption align="center" style="font-size:0.9em;opacity:.6;margin-bottom:-15px"><i>Jupyter Lab</i></figcaption>
    <img width="1112" alt="jupyter-lab-kernel" src="https://github.com/acceleratedscience/open-ad-toolkit/assets/30096303/99e9fda1-e9e3-4533-b375-267eaed5c657">

</figure>

-   Magic commands are implemented by the _openad.py_ or _openad.ipynb_ files, and are invoked by the `%openad` prefix. For example:<br>

        %openad list files

-   An example magic comands to play with DS4SD:<br>

        %openad display all collections

-   An example magic commands to play with RXN

        %openad list rxn models

<br>

## Getting Access to RXN, DS4SD and Tell Me Functionality

Below you find login instructions for RXN and DeepSearch. If you choose to use the `Tell Me` function, you will also need to obtain a OpenAI API account.<br>

<br>

### DS4SD (DeepSearch)

1. First, you'll need to generate an API key on the DeepSearch website.

    - Visit the DeepSearch website and create an account: [deepsearch-experience.res.ibm.com](https://deepsearch-experience.res.ibm.com)<br>
    - Once logged in, click the Toolkit/API icon in the top right hand corner, then open the HTTP section
    - Click the "Generate new API key" button<br>
      <br>
    
<img width="1112" alt="ds4sd-api-key" src="https://github.com/acceleratedscience/open-ad-toolkit/assets/30096303/60ed3108-b0f5-48eb-b0f9-eb5a9a4aa643">

1. Once inside the OpenAD client, you'll be prompted to authenticate when activating the DeepSearch (DS4SD) toolkit. When running `set context ds4sd` :

    - **Hostname:** [https://sds.app.accelerate.science](https://sds.app.accelerate.science)
    - **Email:** Your IBM email
    - **API_key:** The DS4SD API key you obtained following the instructions above.<br>
        
        > **Note:** Your DS4SD auth config file is saved as `~/.openad/ds-auth.ext-v2.json`. If you ever want to reset your DS4SD login information, simply delete this file.<br>

1. You should get a message saying you successfully logged in.

<br>

### RXN

1. First, you'll need to generate an API key on the RXN website.

    - Sign up for an RXN account at [rxn.app.accelerate.science](https://rxn.app.accelerate.science)
    - Obtain your API key by clicking the user profile icon in the top right hand corner and select "My profile".<br>
      <br>
      <img width="1112" alt="rxn-api-key" src="https://github.com/acceleratedscience/open-ad-toolkit/assets/30096303/d00199c2-3a94-4f11-a345-27cbde313732">


1. When setting the context to RXN using `set context rxn` you'll be prompted to create a new auth configuration file:

    - **Hostname:** [https://rxn.app.accelerate.science](https://rxn.app.accelerate.science)<br>
    - **API_key:** The RXN API key you obtained following the instructions above.

1. You should get a message saying you successfully logged in.<br>

    > **Note:** Your RXN auth config file is saved as `~/.openad/rxn-auth.ext-v2.json`. If you ever want to reset your RXN login information, simply delete this file. You can also do this by running `set context ds4sd reset`<br>

<br>

### OpenAI

In order to use the "Tell me" functionality, you will need to create an account with OpenAI. There is a one month free trial.

> **Note:** WatsonX coming soon

1. Go to [platform.openai.com](https://platform.openai.com) and create an account

1. Click on the profile icon in the top right and choose "View API keys"

1. Create a new key

1. Run `tell me` to be prompted for your OpenAI API credentials

1. Your hostname is [https://api.openai.com/v1/models](https://api.openai.com/v1/models)

<img width="1112" alt="openai-api-key" src="https://github.com/acceleratedscience/open-ad-toolkit/assets/30096303/b0ce5207-a1b7-4558-9e57-18aee87baaee">


## Installation for Development

Only follow these instructions if you're contributing to the codebase.

1.  **Step 0: Before you start**<br>
    Ensure you're running Python 3.10.10 or above.

1.  **Step 1: Set up virtual environment** (optional)<br>

        python -m venv ~/ad-venv
        source ~/ad-venv/bin/activate

    > **Note:** To exit the virtual environment, you can run `deactivate`

1.  **Step 2: Installation**<br>
    [Download](https://github.com/acceleratedscience/open-ad-toolkit) or clone the right branch from GitHub:

        git clone -b main https://github.com/acceleratedscience/open-ad-toolkit.git

    > **Note:** To clone a particular branch, replace `main` with your branch name.

    Then, enter the repo's top directory and install the requirements

        cd open-ad-opentoolkit
        pip install -e .

    > **Note:** The `-e` flag stands for "editable". This means that instead of copying the package's files to the Python site-packages directory as in a regular installation, pip creates a symbolic link (symlink) from your package's source code directory into your Python environment.<br><br>
    > This way you can make changes to the source code of the package, and those changes are immediately reflected in your Python environment. You don't need to reinstall the package every time you make a change.

<br>

## Testing a branch

To do a regular install from a particular branch, you can run:

    pip install git+https://github.com/acceleratedscience/open-ad-toolkit.git@<branch_name>

<!--

For screenshots to look good, they should be small and ideally
all the same size. The script below lets you open the URLs in
the right size. Just paste this into the browser console and
press enter.

To take the screenshots with browser UI included on Mac, press
cmd+shift+4 followed by the spacebar, then click the window.
For consistency, stick to Chrome.

- - -

urls = [
    'https://cps.foc-deepsearch.zurich.ibm.com',
    'https://rxn.app.accelerate.science',
    'https://sds.app.accelerate.science',
    'https://platform.openai.com/account/api-keys'
]
for (var i=0; i< urls.length; i++) {
    window.open(urls[i], '_blank', 'width=1000,height=600');
}

-->
