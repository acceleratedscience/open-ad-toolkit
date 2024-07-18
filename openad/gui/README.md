# OpenAD GUI

The OpenAD GUI provides a visual window onto your data, helping you with evaluation and triage.


### How it works

The GUI can be launched from the CLI (in the browser) or a Jupyter Notebook (in an iFrame).
Any command that requires the GUI will start the server, which will then keep on running until the application or Notebook is shut down.

## Components

1. File browser
1. Molecule viewer
2. Molecule set viewer
1. Data viewer (to be implemented)
1. "My molecules" working set
1. Results

<br>

### 1. File Browser

|![test](readme/file-browser.png)|
|---|


|   |   |   |   |   |
|---|---|---|---|---|

<img src= "readme/file-browser.png" alt="File browser" style="border: 1px solid rgba(0,0,0,.1);">

> ![test](readme/file-browser.png)



The file browser lets us open our own proprietary file formats:

- **.mol.json** --> Individual molecule files
- **.molset.json** --> Sets of molecule files

As well as a number of commonly used file formats:

- **.mol**
- **.sdf**
- **.smi**
- **.json**
- **.csv**

Files can easily be opened in your default system app, which is the default for any unsupported file formats.
