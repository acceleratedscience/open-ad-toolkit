<sub>[&larr; BACK](./README.md#openad)</sub>

# OpenAD Plugins

<!-- about_plugin -->
Plugins are how molecular tools and AI models are made available to the OpenAD client. They provide drastically simplified access to a series of advanced tools, and they make it easy for your own Python applications to interface with OpenAD.

Creating your own plugins is easy if you have a basic understanding of Python.<!-- /about_plugin --> Jump to the [Plugin Development Manual](README_plugins_development.md).

<br><br>

## Discovering Plugins

OpenAD comes preloaded with a number of plugins:

- **[DS4SD](#1-ds4sd-deep-search):** Literature knowledge extraction
- **[RXN](#2-rxn):** Forward and retrosynthesis prediction
- **[GT4SD](#3-guardian):** Generative methods and property inference <!-- https://github.com/GT4SD/gt4sd-core -->
- **[BMFM](#4-bmfm):** Biomedical Foundation Models

Additinal available plugins will we published on the OpenAD Portal (coming soon).

### 1. DS4SD (Deep Search)
The DS4SD plugin lets you search for similar molecules, substructures and other data from patents and industry literature.<br>
[ds4sd.github.io](https://ds4sd.github.io/)

    pip install @openad/ds4sd

<details>
<summary>Registration</summary>
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

### 2. RXN
Predict reactions, find retrosynthesis pathways and derive experimental procedures.<br>
[rxn.app.accelerate.science](https://rxn.app.accelerate.science/)

    pip install @openad/rxn

<details>
<summary>Registration</summary>
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

### 3. Guardian
Built on GT4SD, Guardian makes training and using state-of-the-art generative AI models a matter of minutes.<br>
[open.accelerator.cafe](https://open.accelerator.cafe/)

    pip install @openad/guardian



<br><br>

## Installing a Plugin

Installing any non-native plugins can be done with a simple pip install of the plugin's GitHup repository.

    pip install git+https://github.com/some_username/some_plugin.git

<br><br>

## Writing your Own Plugins

Creating your own plugins is quick and easy if you have a basic understanding of Python. Check the [Plugin Development Manual](README_plugins_development.md) to hit the ground running.