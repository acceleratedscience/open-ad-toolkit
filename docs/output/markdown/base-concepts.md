---
title: Base Concepts
layout: home
nav_order: 3
---

<!--

DO NOT EDIT
-----------
This file auto-generated.
To update it, consult instructions:
https://github.com/acceleratedscience/open-ad-toolkit/tree/main/docs

-->

# Base Concepts

## Workspaces

A workspace represents a directory where all of your files, settings and runs are stored. All operations run within a workspace and each workspace comes with its own command history. This allows you to work on unrelated projects without them contaminating eachother.

A default workspace is created on startup, and you can create as many additional workspaces needed. They can live in the designated workspaces directory, or anywhere else on your file system if that's preferred.

To see how to work with workspaces:

    ? workspace

## Plugins

Toolkits are highly specialized applications that OpenAD interfaces with. By interacting with the toolkits through OpenAD, you can bypass a lot of complexity and different inconsistent APIs.

The available tookits for the OpenAD beta are **DS4SD** (DeepSearch) and **RXN** with support for **ST4SD** and **GT4SD** coming soon. For more information about the individual toolkits, click [here]({% link index.md %}#toolkits).

In the future, we hope to expand our toolkit offering with other opensource tools and users will be able to create and customize their own.

To see how to work with plugins:

    ? toolkit

## Context

After installing a toolkit, in order to interact with it you first need to set the context to that toolkit. This will make all commands for this toolkit available.

To see how to switch contexts:

    ? context

## Runs

When working in the terminal, you can record a series of consecutive commands we call runs. They can be replayed at any given time later, removing a lot of the repetitive work.

To see how to work with runs:

    ? run
