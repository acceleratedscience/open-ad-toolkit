# Creating Your Own Toolkit

Integrating your own workflows into OpenAD is relatively straightforward. If you would like to publish your toolkit for the entire OpenAD community, please [get in touch](https://acceleratedscience.github.io/openad-docs/about.html)

## Setup

The toolkit architecture depends on a few basic files to work. You can copy the [DEMO](./DEMO) toolkit directory to hit the ground running.


1. [`metadata.json`](#metadatajson) This file controls the splash screen when the toolkit is launched.
1. [`login.py`](#loginpy) (optional) If your toolkit requires authentication, it can be exposed here.

If you wish to publish your toolkit, a few more files are required:
1. [`description.txt`](descriptiontxt) (optional) Used to train a large language model (LLM)
1. [`oneline_desc.txt`](oneline_desctxt) (optional) This file is required if you wish to publish your toolkit.

### metadata.json

This file is responsible for generating the splash screen. This is displayed when the user enters the toolkit context by running `set context demo`.

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


### login.py

Here you can expose your authentication API. If this file is present, is will be called whenever the user enters the toolkit context by running `set context <toolit_name>`. If this file is not present, authentication will be skipped. (@Phil correct?)

The [`login.py`](./DEMO/login.py) template takes care of success and error handling and ensures a unified user experience across all tookits. Instructions are in the file. You may have to customize it a bit more if your authentication API doesn't follow jwt / host / email conventions. (@Phil what about other variables like username instead of email?)


## LLM Training

Every toolkit should have a `description.txt` file which is used to train the LLM with the toolkit functionality.

All the toolkit commands are listed in the description file, but they are auto-generated. See `/docs/generate-docs.py` for instructions on how to update the description file with the latest version of the toolkit commands.

**Important:** The toolkit's `description.txt` file should contain the following line:

    The following commands are available for this toolkit:

The `generate-docs.py` script will look for this line and replace everything after with the updated commands. If this exact line is not present, the script will abort and throw an error.


<details>
<summary markdown="block"></summary>
<div markdown="block">
</div>
</details>