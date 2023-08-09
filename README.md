# ADCCL

Accelerated Discovery Common Client<br>

_Note: if you're on Mac, please use `pip3` and `python3` instead of `pip` and `python` respectively._

<br>

## Preparation

1.  **Before You Start**<br>
    Ensure you're running Python 3.10.10 or above.

1.  **Step 1: Set up Virtual Environment** (optional)<br>

        python -m venv  ./myenv
        source ./myenv/bin/activate

1.  **Step 2: Installation**<br>
    Download [the repository](https://github.ibm.com/Accelerated-Discovery/ad4e-opentoolkit) or clone it using ssh:

         `https://github.ibm.com/Accelerated-Discovery/ad4e-opentoolkit.git`
     

    Install the requirements:
        cd into source directory then run 
        pip install -e .'

    If you plan to use it inside Jupyter Notebook of Jupyter Labs:<br>
    install ipykernel, which consists of IPython as well <br>
        `pip install ipykernel` <br>
    create a kernel that can be used to run notebook commands inside the virtual environment <br>
        `python -m ipykernel install --user --name=myenv`
    then you can Have the Magic Commands Loaded by default by copying the magic commands into the ipython startup directory for your created profile using 
        `cp ./magic/adccl.py ~/.ipython/profile_[name]/startup`

    if you have not created any profiles pr you just want to use the default you can just use the following command.
        `cp ./magic/adccl.py ~/.ipython/profile_default/startup`


##login details for RXN and Deepsearch, if you choose to use the Tell Me unction you will also need to obtain a Openai API account
deep Search Rep on IBM Network

###DeepSearch
URL: https://cps.foc-deepsearch.zurich.ibm.com/
login w3 id
obtain API key when logged into webpage
![Screenshot 2023-08-02 at 5 00 05 pm](https://media.github.ibm.com/user/225313/files/76807d43-262c-4ff0-969f-9086b15613ba)

###RXN
url  https://rxn.app.accelerate.science/rxn/home
login w3id

obtain API key by clicking the user profile in the top right hand corner
![Screenshot 2023-08-02 at 5 03 01 pm](https://media.github.ibm.com/user/225313/files/26d30714-f028-4f97-844c-82a434f9e0d8)


#OpenAI
you will need the details from the openai API loging. there is a free 1 monthtrial
https://platform.openai.com

![Screenshot 2023-08-02 at 5 04 29 pm](https://media.github.ibm.com/user/225313/files/50f34891-dd0f-4650-9548-45631606a0d1)

<br>

## Getting Started

-   **Enter the Shell Environment**

     Enter 'adccl' from any directory

    Alternatively for jupyter lab, you can:
    run jupyter lab or jupyter notebooks and select a notebook under you virtual environment or from the install directory run
    
    `./notebooks/jupyter lab Table_of_Contents.ipynb`


    NOTE: by launching jupyter this way it will automatically launch the trial notebooks.
    <br>
      
    ![Landing](readme/screenshot-landing.png)

-   **Exit**<br>
    Hit `ctrl+c` or type:

        exit

-   **Installing Toolkits**<br>
    You can install the `DS4SD`, `GT4SD`, `ST4SD` and `RXN` toolkits, however please note that at this time, they are meant as placeholders and only `DS4SD` and `RXN` supports experimental functionality.

        add toolkit ds4sd
        add toolkit rxn
        
        you will need to obtain login keys for both of these which cn be done at 
 
 For Deep Search, you will need to log onto the Zurich system on the IBM VPN and logon then click the programming icon in the top right corner to get your API key
        
        `https://cps.foc-deepsearch.zurich.ibm.com`
 
 ![Screenshot 2023-08-04 at 1 18 05 pm](https://media.github.ibm.com/user/225313/files/637e3cc6-6ae1-4d99-8b3c-294cb259df9e)

       
For RXN go to your Profile once you have created a logn and click in the top right hand corner to generate API key (profile)
         `https://rxn.app.accelerate.science/rxn/home`
        ![Screenshot 2023-08-04 at 1 14 26 pm](https://media.github.ibm.com/user/225313/files/14261abf-5839-4e6a-92f6-1dc6ed9803b7)


For the "Tell Me" functionality we currently have OPENAI , you will need to setup an Open AI API account and create and Save a API key
        
        ![Screenshot 2023-08-02 at 5 04 29 pm](https://media.github.ibm.com/user/225313/files/3a1a83cc-b07e-4328-9035-e956a45e629d)

        

-   **Running as a Bash Command**<br>
    To run any commands as a bash command, make sure to prepend any quotes with `\`.

        ./ad4e-opentoolkit/adccl show molecules using file \'base_molecules.sdf\'

-   **Working with Notebooks**


    -   Magic commands are implemented by the adccl.ipynb file and are invoked by the `%adccl` prefix. For example:

            %adccl list files

    -   Open the table of contents to get an introduction and be taken through step by step how to use the tool.
        ![Notebook table of contents](readme/notebook-toc.png)

    -   Play with Deep Search using the magic comands
        ![Notebook DS4SD](readme/notebook-ds4sd.png)

    -   You can also copy commands from a terminal version and work in your own terminal or the Jupyter Lab built-in Terminal
