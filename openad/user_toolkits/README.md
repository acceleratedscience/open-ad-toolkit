# Creating Your Own Toolkit <!-- no toc -->

Integrating your own workflows into OpenAD is relatively straightforward. If you would like to publish your toolkit for the entire OpenAD community, please [get in touch](https://acceleratedscience.github.io/openad-docs/about.html)

### Table of Contents <!-- no toc -->
- [Creating Your Own Toolkit ](#creating-your-own-toolkit-)
    - [Table of Contents ](#table-of-contents-)
  - [Setup](#setup)
    - [metadata.json](#metadatajson)
    - [login.py](#loginpy)
  - [Publishing a Toolkit](#publishing-a-toolkit)
  - [description.txt](#descriptiontxt)
  - [oneline\_desc.txt](#oneline_desctxt)


## Setup

The toolkit architecture depends on a few basic files to work. You can copy the [DEMO](./DEMO) toolkit directory to hit the ground running.


1. [`metadata.json`](#metadatajson)
2. [`login.py`](#loginpy) (optional)

<br>

### metadata.json

This mandatory file is responsible for generating the toolkit's splash screen. This is displayed when the user enters the toolkit context by running `set context <toolkit_name>`.

    {
        "banner": "DEMO",
        "title": "This is a Demo Toolkit",
        "author": "Jane Doe",
        "version": "0.0.1",
        "intro": "This toolkit is meant as a demonstration on how to set up a toolkit. This intro paragraph should contain a brief description of what the toolkit does, ideally not much longer than ~250 characters.",
        "commands": {
            "hello world": "Say hello",
            "demo docs": "Read the docs in your browser"
        }
    }

The splash screen generated from the JSON file above looks like this:

![demo-splash-page](readme/demo-splash.png)

<br>

### login.py

Here you can expose your authentication API. If this file is present, is will be called whenever the user enters the toolkit context by running `set context <toolit_name>`. If this file is not present, authentication will be skipped. (@Phil correct?)

The [`login.py`](./DEMO/login.py) template takes care of success and error handling and ensures a unified user experience across all tookits. Instructions are in the file. You may have to customize it a bit more if your authentication API doesn't follow jwt/host/email/api_key conventions. (@Phil what about other variables like username instead of email?)

<br>

## Publishing a Toolkit

If you wish to publish your toolkit to the OpenAD comminity, a few more files are required.
1. [`description.txt`](descriptiontxt) Used to train a large language model (LLM).
2. [`oneline_desc.txt`](oneline_desctxt) Used to describe your toolkit in the list of available toolkits.

<br>

## description.txt

The `description.txt` file is used to train the LLM with the toolkit functionality. It should contain a detailed description of how your toolkit works and what it is meant to achieve, written in an unambiguous way. You can look at the other toolkits for inspiration.

At the bottom of your file, on a separate line, you should include the following line, verbatim:

    The following commands are available for this toolkit:

Then you should run a script that gathers all your toolkit commands and lists them at the bottom of the description file. The script will look for this line described above and replace everything after with the updated commands. If this exact line is not present, the script will abort and throw an error.

    python openad/user_toolkits/<toolkit_name>/description_update.py

<br>

## oneline_desc.txt

This file contains a very brief description of the toolkit, using only 4-5 words. This wil be displayed when listing available toolkits.