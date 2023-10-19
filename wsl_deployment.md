# Installing OpenAD on Windows using WSL Ubuntu

Download the this respository and place it in your Downloads directory on Windows.

-   **Step 1**: Install Ubuntu 22.0.4 or higher as Windows WSL<br>
    To create the "openad" user, open the URL below and follow the steps.<br>
    [https://linuxconfig.org/ubuntu-22-04-on-wsl-windows-subsystem-for-linux](https://linuxconfig.org/ubuntu-22-04-on-wsl-windows-subsystem-for-linux)

-   **Step 2**: Set up the Python environment

        sudo add-apt-repository ppa:deadsnakes/ppa
        sudo apt update
        sudo apt install python3.11-full
        sudo apt install python3-pip
        sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 100

-   **Step 3**: Set up your virtual environment

        python3 -m venv ~/ad-venv
        source ~/ad-venv/bin/activate
        pip3 install jupyter
        pip3 install ipykernel
        python3 -m ipykernel install --user --name=ad-kernel

-   **Step 4**: Substitute Windows username

        ln -s /mnt/c/Users/`<windows_user>`/Downloads ~/downloads
        mkdir -p ~/.ipython/profile_default/startup
        cd ~/downloads/ad4e-opentoolkit/
        pip install .
        init_magic
        init_examples
        jupyter lab ~/openad_notebooks/Table_of_Contents.ipynb
