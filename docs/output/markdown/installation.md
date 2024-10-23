---
title: Installation
layout: home
nav_order: 2
---

<!--

DO NOT EDIT
-----------
This file auto-generated.
To update it, consult instructions:
https://github.com/acceleratedscience/open-ad-toolkit/tree/main/docs

-->

## Quick Install

> [!IMPORTANT]  
> This will install OpenAD in your global space. If you wish to use a virtual environment (recommended), please refer to the [Installation](/README_installation.md) page.

    pip install openad
    openad

Get started with Jupyter Notebook examples:

    init_magic
    init_examples
    jupyter lab ~/openad_notebooks/Table_of_Contents.ipynb

If you get an error when running `init_magic`, you may first need to setup the default iPython profile for magic commands.

    ipython profile create

<br><br>

## What's New

-   `%Openadd` has been added to the magic commands for commands that return data.
-   Upgraded SkyPilot to 0.6.0
-   Support for deploying in OpenShift AI/Open Data hub workbench or Podman/Docker image. [See the workbench repo](https://github.com/acceleratedscience/openad_workbench).
-   Support for application API
-   New property and dataset generation services. See [OpenAD Models Service](/README_models-service.md)
    -   GT4SD Generation Services
    -   GT4SD Property Services
    -   GT4SD MoleR Generation
    -   GT4SD Molformer

    Pre-Requisite is that you have a AWS Account and can launch your own EC2 Instances Or someone else can launch them for you and you can catalog a Remote Service via URL.
