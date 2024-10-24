<sub>[&larr; BACK](../#openad)</sub>

# OpenAD Installation

#### Notes
- **Simplified installation:** For install without virtual environment, see [Quick Install](/README.md#quick-install)
- **Contributors:** Skip to [Installation for Development](#installation-for-development)<br>
- **Linux users:** Check the [Linux Notes](#linux-notes)<br>
- **Poetry:** If you prefer Poetry, you can run the setup wizard instead: `poetry add openad`

---

<!-- Note: step 1 & 2 are repeated, make sure any updates are done in both places -->
1.  **Step 0: Before you start**<br>
    Ensure you're running Python 3.10.10 or 3.11. See [instructions below](#upgrading-python).

    To see what version you are running:

        python -V

    > **Note:** Due to an issue with one of our dependencies, Python 3.12 is not yet supported.

2.  **Step 1: Set up your virtual environment** (recommended)<br>

        python -m venv ~/ad-venv
        source ~/ad-venv/bin/activate

    > **Note:** Use `python3` on macOS.
    > **Note:** To exit the virtual environment, you can run `deactivate`

3.  **Step 2: Install OpenAD**

        pip install openad

<br>

# Installing on Windows

In order to run OpenAD on Windows 11, you will need to install the Ubuntu WSL package ("Windows Subsystem for Linux").

<br>

## Before You Start

-   **Verify Windows version**<br>
    To check if you are running Windows 11 or later, press `Win` + `R`, type "winver", and press `Enter`. A window will open showing your Windows version.

-   **Verify WSL**<br>
    To check if you already have WSL installed, run `wsl -l -v` into the terminal. To see more information about your current version of Ubuntu, run `lsb_release -a`

<br>

## Installing WSL

Install WSL and create a user called 'openad' or one of your choosing.

    wsl --install Ubuntu-22.04

**Optional:** To setup an Ubuntu Python environment from scratch, continue to <a href="#linux-notes">Linux Notes</a>

<br>

# Appendix

## Upgrading Python

There's many ways to install or upgrade Python. We'll use `pyenv`.

1.  **Install pyenv**

        curl https://pyenv.run | bash
    
1.  **Set up your shell environment for Pyenv**<br>
    Detailed instructions cam be found [here](https://github.com/pyenv/pyenv?tab=readme-ov-file#set-up-your-shell-environment-for-pyenv). If you're using Zsh, you can run the commands below:

        echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.zshrc
        echo '[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.zshrc
        echo 'eval "$(pyenv init -)"' >> ~/.zshrc
    
1.  **Reboot your shell**<br>
    You can either open a new window or run:

        exec $SHELL

1. **Install Python**
    Please note that OpenAD requires **Python 3.10** or **3.11**. Due to an issue with one of our dependencies, Python 3.12 is not yet supported.
    
        pyenv install 3.11

2.  **Activate this version of Python**
    If you wish to set this version as the default:

        pyenv global 3.11
        
    Alternatively, if you only wish to activate it in the current shell:

        pyenv shell 3.11

## Linux Notes

If you wish to setup an Ubuntu Python environment from scratch, run:

    sudo add-apt-repository ppa:deadsnakes/ppa
    sudo apt update
    sudo apt install python3.11-full
    sudo apt install python3-pip
    sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 100
    sudo pip install pip --upgrade

You will need to restart your Linux session before running `pip install openad` so that the python libraries are in your path.

If you get an error when running `init_magic`, you may first need to setup the default iPython profile for magic commands.

    ipython profile create