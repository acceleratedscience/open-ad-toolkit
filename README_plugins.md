<sub>[&larr; BACK](./README.md#openad)</sub>

# OpenAD Plugins

<!-- about_plugin -->
Plugins are the way molecular tools and AI models are made available to the OpenAD client. Thanks to a unified language, accessing these tools through OpenAD lets you bypass a lot of complexity.

OpenAD comes preloaded with a number of plugins for literature knowledge extraction (DS4SD), forward and retrosynthesis prediction (RXN) as well as generative methods and property inference (GT4SD).

You can create your own plugins, and the publicly available plugins will soon include a much larger variety of open-source tools.

Note: Plugins are currently referred to as "toolkits" by the commands, however this language will be updated soon.
<!-- /about_plugin -->

> [!IMPORTANT]
> When running commands from Jupyter, prepend them with `%openad`

<br><br>

## Registration

Before you can interact with the plugins, you'll need to register with each individual plugin.

<details>
<summary>Register with DS4SD (Deep Search)</summary>
<div markdown="block">

1. First, you'll need to generate an API key on the Deep Search website.

    - Visit the Deep Search website and create an account:<br>
      [deepsearch-experience.res.ibm.com](https://deepsearch-experience.res.ibm.com)<br>
    - Once logged in, click the `Toolkit / API` icon in the top right hand corner, then open the HTTP section
    - Click the "Generate new API key" button<br>
      <br>
      <a href="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/ds4sd-api-key.png" target="_blank"><img src="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/ds4sd-api-key.png" /></a>

2. Once inside the OpenAD client, you'll be prompted to authenticate when activating the Deep Search (DS4SD) toolkit. When running `set context ds4sd` :

    - **Hostname:** Default: [https://sds.app.accelerate.science](https://sds.app.accelerate.science)
    - **Email:** Your email
    - **API_key:** The DS4SD API key you obtained following the instructions above.

3. You should get a message saying you successfully logged in.

    > **Note:** Your DS4SD auth config file is saved as `~/.openad/deepsearch_api.cred`. If you ever want to reset your DS4SD login information you can run `set context ds4sd reset`, or you can delete this file.<br>

</div>
</details>

<details>
<summary>Register with RXN</summary>
<div markdown="block">

1. First, you'll need to generate an API key on the RXN website.

    - Sign up for an RXN account at [rxn.app.accelerate.science](https://rxn.app.accelerate.science)
    - Obtain your API key by clicking the user profile icon in the top right hand corner and select "Account", then select the "My keys" tab.<br>
      <br>
      <a href="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/rxn-api-key.png" target="_blank"><img src="https://raw.githubusercontent.com/acceleratedscience/open-ad-toolkit/main/assets/rxn-api-key.png" /></a>

2. When setting the context to RXN using `set context rxn` you'll be prompted to create a new auth configuration file:

    - **Hostname:** Default: [https://rxn.app.accelerate.science](https://rxn.app.accelerate.science)<br>
    - **API_key:** The RXN API key you obtained following the instructions above.

3. You should get a message saying you successfully logged in.<br>

    > **Note:** Your RXN auth config file is saved as `~/.openad/rxn_api.cred`. If you ever want to reset your RXN login information you can run `set context rxn reset`, or you can delete this file.<br>

</div>
</details>

<br><br>

## Adding a Toolkit

First install the toolkit, then set the context to this toolkit.

    add toolkit ds4sd
    set context ds4sd

<br><br>

## Sample Commands

    # DS4SD
    display all collections

    # RXN
    list rxn models

<br><br>

## Running Bash Commands (CLI)

To run a command in bash mode, prepend it with `openad` and make sure to escape quotes.

    openad show molecules using file \'base_molecules.sdf\'