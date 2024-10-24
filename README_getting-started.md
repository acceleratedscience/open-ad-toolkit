<sub>[&larr; BACK](./README.md#openad)</sub>

# Getting Started with OpenAD

---

- If you haven't done so yet, [install OpenAD](README_installation.md) first.
- When installing on macOS without a virtual environment, you may need to use `python3` and `pip3` instead of `python` and `pip`.
- When updating to OpenAD `0.4.0` or above, first remove all toolkits by runnning `list toolkits` and then `remove toolkit <toolkit_name>`.

---

<br>

# Getting Started - CLI

-   **Enter the virtual environment**

    > **Note:** If you just installed OpenAD, you probably already activated the virtual environment.

        source ~/ad-venv/bin/activate

-   **Enter the command shell**

        openad

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