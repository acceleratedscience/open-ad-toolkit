# ADCCL

Accelerated Discovery Common Client<br>

_Note: if you're on Mac and not installing into a virtual environment you may need use `pip3` and `python3` instead of `pip` and `python` respectively._

<br>

## Preparation

1.  **Before You Start**<br>
    Ensure you're running Python 3.10.10 or above.

1.  **Step 1: Set up Virtual Environment** (optional)<br>

       `python -m venv  ./myenv` <br> 
        `source ./myenv/bin/activate` <br>

1.  **Step 2: Installation**<br>
    Install the requirements:<br>

   `pip install   git+ssh://git@github.ibm.com/Accelerated-Discovery/ad4e-opentoolkit.git` <br>
        
   To enter the command Shell simply enter `adccl` from the command line or to run as a  bash command `adccl <ad4e toolkit command>` <br>
    
   Use `?` for help from the command line <br><br>
    

   If you plan to use it inside Jupyter Notebook of Jupyter Labs:<br>
    install ipykernel, which consists of IPython as well <br>
        `pip install ipykernel` <br>
        <br>
    Create a kernel that can be used to run notebook commands inside the virtual environment <br>
        `python -m ipykernel install --user --name=myenv`
    
   Then you can Have the Magic Commands Loaded by default by copying the magic commands into the ipython startup directory for your created profile using <br>
        `init_magic` <br>
        Or to copy it for use within only chose notebook do the following to copy it to the current directory and run with !run adccl.py to initiate Magic commands <br>
        `init_magic . ` <br>
        Or to copy it to a ipython custom profile <br>
        `init_magic myprofile`
        
        **Note:** if you do not run this command you will need to run it he init_examples and use `init_examples` command and retrieve the file  `~/adccl_notebooks/adccl.ipynb` and executer it in each notebook you run with `!run adccl.ipynb` to activate the magic commands.
    

   
## Getting access To RXN, DeepSearch and Tell Me founctionality

login details for RXN and Deepsearch, if you choose to use the `Tell Me` function you will also need to obtain a Openai API account
deep Search Rep on IBM Network <br>

### DeepSearch <br>
URL: https://cps.foc-deepsearch.zurich.ibm.com/ <br>
login w3 id <br>
to address this you will need to have your VPN available <br>
obtain API key when logged into webpage <br>
![Screenshot 2023-08-02 at 5 00 05 pm](https://media.github.ibm.com/user/225313/files/76807d43-262c-4ff0-969f-9086b15613ba)
<br>Note: to reset login delete the file `~/.adccl/ds-auth.ext-v2.json`<br>
### RXN <br>
for rxn you can sign up to <br>
url  https://rxn.app.accelerate.science/rxn/home <br>
login w3id <br>

obtain API key by clicking the user profile in the top right hand corner <br>
![Screenshot 2023-08-02 at 5 03 01 pm](https://media.github.ibm.com/user/225313/files/26d30714-f028-4f97-844c-82a434f9e0d8)
<br>Note: to reset login delete the file `~/.adccl/rxn-auth.ext-v2.json`<br>

### OpenAI
you will need the details from the openai API loging. there is a free 1 month trial <br>
https://platform.openai.com  <br>

![Screenshot 2023-08-02 at 5 04 29 pm](https://media.github.ibm.com/user/225313/files/50f34891-dd0f-4650-9548-45631606a0d1)

<br>

## Getting Started

-   **Enter the Shell Environment**

     Enter 'adccl' from any directory

    Alternatively for jupyter lab, you can:
    run `jupyter lab` or `jupyter notebook` and select a notebook under you virtual environment <br>

    
    Make sure the python kernel used is your virtual envrionment kernel, you can select this in the top right hand corner of the browser. <br>
    ![Screenshot 2023-08-10 at 1 06 07 pm](https://media.github.ibm.com/user/225313/files/f4ab9f61-dc34-4a33-9a8d-b5cc64b00dbe)

    NOTE: by launching jupyter this way it will automatically launch the trial notebooks.<br>
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
 
 ## Accelerated Discovery System Access
 For Deep Search, you will need to log onto the Zurich system on the IBM VPN and logon then click the programming icon in the top right corner to get your API key
        
        `https://cps.foc-deepsearch.zurich.ibm.com`
 
 ![Screenshot 2023-08-04 at 1 18 05 pm](https://media.github.ibm.com/user/225313/files/637e3cc6-6ae1-4d99-8b3c-294cb259df9e)

       
For RXN go to your Profile once you have created a logn and click in the top right hand corner to generate API key (profile)<br>

        `https://rxn.app.accelerate.science/rxn/home`
  ![Screenshot 2023-08-04 at 1 14 26 pm](https://media.github.ibm.com/user/225313/files/14261abf-5839-4e6a-92f6-1dc6ed9803b7)


For the "Tell Me" functionality we currently have OPENAI, you will need to setup an Open AI API account and create and Save a API key www.openai.com <br>
 ***Note*** WatsonX coming soon !<br>

   ![Screenshot 2023-08-02 at 5 04 29 pm](https://media.github.ibm.com/user/225313/files/3a1a83cc-b07e-4328-9035-e956a45e629d)

        

-   **Running as a Bash Command**<br>
    To run any commands as a bash command, make sure to prepend any quotes with `\`.

        adccl show molecules using file \'base_molecules.sdf\'

-   **Working with Notebooks**

    -   To Start using the default Notebooks run the following command <br>
        `init_notebooks`. this will create the directory ~/adccl_notebooks for you to work from <br>
        run the below to start playing ! <br>
        `jupyter lab ~/adccl_notebooks/Table_of_Contents.ipynb` <br>
       
    -   Magic commands are implemented by the adccl.py or adccl.ipynb file and are invoked by the `%adccl` prefix. For example:<br>

            %adccl list files
        
        if you use your virtual envrinment kernel , as per above, they and ran the init_magic command the magic commands should be enabled already.<br>
        
    -   Open the table of contents to get an introduction and be taken through step by step how to use the tool.<br>
        ![Notebook table of contents](readme/notebook-toc.png)

    -   Play with Deep Search using the magic comands<br>
        ![Notebook DS4SD](readme/notebook-ds4sd.png)

    
