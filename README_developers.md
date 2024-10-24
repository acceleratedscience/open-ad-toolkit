<sub>[&larr; BACK](./README.md#openad)</sub>

# OpenAD for Developers

OpenAD is fully open source and we encourage contributions.

If you have any questions in the meantime, please [reach out]().

<br><br>

## Developing Plugins

Stay tuned for documentation on how to build plugins, enabling you to integrate your own tools into the OpenAD eco-system.

<br><br>

## Installation for Development

<details>
<summary>Install using the setup wizard (uses poetry)</summary>
<div markdown="block">

1.  **Step 1: Download the repo**

        git clone https://github.com/acceleratedscience/open-ad-toolkit.git

    > **Note:** To download a specific branch, you can run instead:<br> > `git clone -b <branch_name> https://github.com/acceleratedscience/open-ad-toolkit.git`

1.  **Step 2: Launch the setup wizard**

        cd open-ad-toolkit
        ./setup.sh

</div>
</details>

<details>
<summary>Install using pip</summary>
<div markdown="block">

<!-- Note: step 1 & 2 are repeated, make sure any updates are done in both places -->
1.  **Step 0: Before you start**<br>
    Ensure you're running Python 3.10.10 or 3.11. See [instructions below](#upgrading-python).

    To see what version you are running:

        python -V

    > **Note:** Due to an issue with one of our dependencies, Python 3.12 is not yet supported.

1.  **Step 1: Set up your virtual environment** (recommended)<br>

        python -m venv ~/ad-venv
        source ~/ad-venv/bin/activate

    > **Note:** Use `python3` on macOS.
    > **Note:** To exit the virtual environment, you can run `deactivate`

1.  **Step 2: Download the repo**

        git clone https://github.com/acceleratedscience/open-ad-toolkit.git

    > **Note:** To download a specific branch, you can run instead:<br> > `git clone -b <branch_name> https://github.com/acceleratedscience/open-ad-toolkit.git`

1.  **Step 2: Install OpenAD**

        cd open-ad-toolkit
        pip install -e .

    > **Note:** The -e flag stands for "editable". This means that instead of copying the package's files to the Python site-packages directory as in a regular installation, pip creates a symbolic link (symlink) from your package's source code directory into your Python environment.<br>This way you can make changes to the source code of the package, and those changes are immediately reflected in your Python environment. You don't need to reinstall the package every time you make a change.

</div>
</details>

<br><br>

## Testing a branch

To do a regular install from a particular branch, you can run:

    pip install git+https://github.com/acceleratedscience/open-ad-toolkit.git@<branch_name>