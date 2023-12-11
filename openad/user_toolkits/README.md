# Creating Your Own Toolkit

Integrating your own workflows into OpenAD is relatively straightforward. If you would like to publish your toolkit for the entire OpenAD community, please [get in touch](https://acceleratedscience.github.io/openad-docs/about.html)

## Setup

The toolkit architecture depends on a few basic files to work. You can copy the [DEMO](./DEMO) toolkit to hit the ground running.


- `metadata.json` This file controls the splash screen when the toolkit is launched.
- `login.py` (optional) If your toolkit requires authentication, it can be exposed here.
- `oneline_desc.txt` (optional) If you wish to publish your toolkit, this is how it will be listed.

### metadata.json

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


## LLM Training

Every toolkit should have a `description.txt` file which is used to train the LLM with the toolkit functionality.

All the toolkit commands are listed in the description file, but they are auto-generated. See `/docs/generate-docs.py` for instructions on how to update the description file with the latest version of the toolkit commands.

**Important:** The toolkit's `description.txt` file should contain the following line:

    The following commands are available for this toolkit:

The `generate-docs.py` script will look for this line and replace everything after with the updated commands. If this exact line is not present, the script will abort and throw an error.