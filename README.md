# OpenAD

**Open Accelerated Discovery Client**<br>
[Project homepage](https://pages.github.ibm.com/Accelerated-Discovery/ad4e-opentoolkit/)

---

<br>

## Notes

-   Only available for Linux and MacOS
-   Currently only the OpenAi API is available for the _Tell Me_ Function (WatsonX coming soon)
-   If you're on Mac and not installing into a virtual environment, you may need use `pip3` and `python3` instead of `pip` and `python` respectively.<br>

<br>

## Installation

1.  **Step 0: Before you start**<br>
    Ensure you're running Python 3.10.10 or above.

1.  **Step 1: Set up virtual environment** (optional)<br>

        python -m venv ./openad
        source ./openad/bin/activate

    > **Note:** To exit the virtual environment, you can run `deactivate`

    ![](https://placehold.co/20x20/dd0000/dd0000.png) <span style="color:#d00">Internal note: we should change the virtual env name to something that's descriptive, because just calling it "openad" obfuscates what is happening here. I kept is as-is to avoid inconsistencies with other content. - Moenen</span>

1.  **Step 2: Installation**<br>

        pip install git+ssh://git@github.ibm.com/Accelerated-Discovery/ad4e-opentoolkit.git

    > _**Note:** Before pip installing from git, ensure you have ssh set up for git install, otherwise you can download the repository and run `pip install .` from the ad4e-opentoolkit top directory._

1.  **Launch**<br>
    To enter the command shell, simply enter `openad` from the command line.

    > _**Notes:**<br>
    > • Alternatively, you can run bash commands as such: `openad <command>`<br>
    > • To see available commands, run `?`_

<br>

## Installation for Development

Only follow these instructions if you're contributing to the codebase.

1.  **Step 0 & 1**<br>
    See main installation instructions above.

1.  **Step 2: Installation**<br>

    -   **Option A:** Install a particular branch:

            pip install git+ssh://git@github.ibm.com/Accelerated-Discovery/ad4e-opentoolkit.git@<branch_name>

    -   **Option B:** Install the entire git repository, ready for contributing:

        [Download](https://github.ibm.com/Accelerated-Discovery/ad4e-opentoolkit) or clone the right branch from GitHub:

            git clone -b <branch_name> git@github.ibm.com:Accelerated-Discovery/ad4e-opentoolkit.git

        Then, from the _ad4e-opentoolkit_ top directory, run:

            pip install -e .

        > The `-e` flag stands for "editable". This means that instead of copying the package's files to the Python site-packages directory as in a regular installation, pip creates a symbolic link (symlink) from your package's source code directory into your Python environment.

        > This way you can make changes to the source code of the package, and those changes are immediately reflected in your Python environment. You don't need to reinstall the package every time you make a change.

<br>

## Getting Started - CLI

-   **Entering the Shell Environment**<br>
    Run from any directory:

        openad

    ![Landing](readme/screenshot-landing.png)

-   **Exiting the Shell Environment**<br>
    Hit `ctrl+c` or run:

        exit

-   **Installing Toolkits**<br>
    You can install the `DS4SD`, `GT4SD`, `ST4SD` and `RXN` toolkits, however please note that at this time, only `DS4SD` and `RXN` support experimental functionality while the others are meant as placeholders.

        add toolkit ds4sd
        add toolkit rxn

-   **Running Bash Commands**<br>
    To run any commands as a bash command, make sure to prepend any quotes with `\`.

        openad show molecules using file \'base_molecules.sdf\'

<br>

## Getting Started - Jupyter

### Jupyter Setup

If you plan to use this application inside Jupyter Notebook of JupyterLab, you should set it up as follows:

1.  Make sure your virtual environment was activated (`source ./openad/bin/activate`), as described under [Installation](#installation) on top.

1.  Install ipykernel, which includes IPython:

        pip install ipykernel

1.  Create a kernel that can be used to run Notebook commands inside the virtual environment:

        python -m ipykernel install --user --name=openad

1.  Initiate the magic commands.

    -   **Option A: Across all notebooks** (recommended)<br>
        This copies the magic commands into the iPython startup directory for your created profile:

            init_magic

    -   **Option B: Within a single notebook**<br>
        If you wish to use the magic commands in a single notebook only, you can run instead:

            init_magic .

        Then to initiate the magic commands, run:

            run openad.py

    -   **Option C: Within a custom iPython profile**<br>
        This would install the magic commands within the iPhython custom profile called 'myprofile'

            init_magic myprofile

    > **Note:** If you don't want to install magic commands in your jupyter default profile, you can initiate the magic commands manually per notebook.<br>
    > • Run `init_examples` to copy the Jupyter example files to your home directory.<br>
    > • In each Notebook where you wish to use the magic commands, run `run openad.ipynb` first. This executes the file `~/openad_notebooks/openad.ipynb` which activates the magic commands.

<br>

### Jupyter Launch

> **Note:** The commands described below should be run from the regular CLI, not from openad command shell.

-   Make sure the `~/openad_notebooks` directory with example Notebooks exists (it was created by running `init_magic` or `init_examples` earlier).

-   Open the table of contents to get an introduction and be taken through step by step how to use the tool.

        jupyter lab ~/openad_notebooks/Table_of_Contents.ipynb

    > **NOTE:** By launching Jupyter this way it will automatically launch the trial notebooks.

-   Make sure to select your virtual environment as the Python kernel. In JupyterLab, you can select this in the top right-hand corner. In Juyter Notebook, you can select this under _Kernel > Change Kernel_. If you don't see your virtual environment, make sure you followed the Jupyter Setup instructions listed above.

<br>
<figure>
    <figcaption align="center" style="font-size:0.9em;opacity:.6;margin-bottom:-15px"><i>Jupyter Lab</i></figcaption>
    <img src="readme/jupyter-lab-kernel.png">
</figure>

<figure>
    <figcaption align="center" style="font-size:0.9em;opacity:.6;margin-bottom:-15px"><i>Jupyter Notebook</i></figcaption>
    <img src="readme/jupyter-notebook-kernel.png">
</figure>

-   Magic commands are implemented by the _openad.py_ or _openad.ipynb_ file and are invoked by the `%openad` prefix. For example:<br>

        %openad list files

    If you are using your virtual envrinment kernel as per instructions above, and if you've run the `init_magic` command, the magic commands should be enabled already.<br>

-   Some example magic comands to play with Deep Search:<br>

        %openad exec display_collection(domain='Material Science')

-   Some example magic commands to play with RXN

        %openad list rxn models

<br>

## Getting Access to RXN, DeepSearch and Tell Me Functionality

Below you find login instructions for RXN and DeepSearch. If you choose to use the `Tell Me` function, you will also need to obtain a OpenAI API account.<br>

<br>

### DeepSearch

1. First, you'll need to generate an API key on the DeepSearch website.

    - **IBM Users:** Activate the Cisco Secure Client VPN to access the IBM-internal edge version of DeepSearch: [cps.foc-deepsearch.zurich.ibm.com](https://cps.foc-deepsearch.zurich.ibm.com)
    - **Non-IBM Users:** Access the public version of DeepSearch: [deepsearch-experience.res.ibm.com](https://deepsearch-experience.res.ibm.com)
    - Once logged in, click the Toolkit/API icon, then open the HTTP section
    - Click the "Generate new API key" button<br>
      <br>
      ![Landing](readme/ds4sd-api-key.png)

1. When setting the context to DeepSearch using `set context ds4sd` you'll be prompted to create a new auth configuration file:

    - **Hostname:** [https://cps.foc-deepsearch.zurich.ibm.com](https://cps.foc-deepsearch.zurich.ibm.com)<br>
    - **Email:** Your IBM email<br>
    - **API_key:** The DS4SD API key you obtained following the instructions above.

1. You should get a message saying you successfully logged in.

    > **Note:** Your DS4SD auth config file is saved as `~/.openad/ds-auth.ext-v2.json`. If you ever want to reset your DS4SD login information, simply delete this file.<br>

<br>

### RXN

1. First, you'll need to generate an API key on the RXN website.

    - Sign up for an RXN account at [rxn.app.accelerate.science](https://rxn.app.accelerate.science)
    - Obtain your API key by clicking the user profile icon in the top right hand corner and select "My profile".<br>
      <br>
      ![Landing](readme/rxn-api-key.png)

1. When setting the context to RXN using `set context rxn` you'll be prompted to create a new auth configuration file:

    - **Hostname:** [https://rxn.app.accelerate.science](https://rxn.app.accelerate.science)<br>
    - **API_key:** The RXN API key you obtained following the instructions above.

1. You should get a message saying you successfully logged in.

    > **Note:** Your RXN auth config file is saved as `~/.openad/rxn-auth.ext-v2.json`. If you ever want to reset your RXN login information, simply delete this file.<br>

<br>

### OpenAI

In order to use the "Tell me" functionality, you will need to create an account with OpenAI. There is a one month free trial.

> **Note:** WatsonX coming soon

1. Go to [platform.openai.com](https://platform.openai.com) and create an account

1. Click on the profile icon in the top right and choose "View API keys"

1. Create a new key

1. Run `tell me` to be prompted for your OpenAI API credentials

1. Your hostname is [https://api.openai.com/v1/models](https://api.openai.com/v1/models)

![Landing](readme/openai-api-key.png)

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
