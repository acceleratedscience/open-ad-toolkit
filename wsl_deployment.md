
# Installing Open Ad on Windows using WSL Ubuntu
Download the attached zip file and place it in your Downloads directory on Windows


Step 1:  Install Ubuntu 22.0.4 or higher as Windows WSL 
    create "opened" as the user
    open the below URL and follow the steps
    https://linuxconfig.org/ubuntu-22-04-on-wsl-windows-subsystem-for-linux

Step 2: run the below to setup the python environment 

`sudo add-apt-repository ppa:deadsnakes/ppa`
`sudo apt update`
`sudo apt install python3.11-full`
`sudo apt install python3-pip`
`sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 100`


Step 3: setup python Virtual Environment
 `python3 -m venv ~/openad`
 `source ~/openad/bin/activate`
 `pip3 install jupyter`
 `pip3 install ipykernel`
 `python3 -m ipykernel install --user --name=openad`

Step 4:
 Substitute Windows User Name.
 `ln -s /mnt/c/Users/`<windows_user>`/Downloads ~/downloads`
 
 `mkdir -p ~/.ipython/profile_default/startup`

 `cd ~/downloads/ad4e-opentoolkit/`

 `pip install .`
 `init_magic`
 `init_examples`
 `jupyter lab ~/openad_notebooks/Table_of_Contents.ipynb`