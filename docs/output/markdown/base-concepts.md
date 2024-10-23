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

Your workspace represents an isolated environment where your molecules, molecule sets and other files related to your research are stored.

Each workspace has its own command history and corresponds with a directory stored in ~/openad/workspaces, allowing you to work on multiple isolated projects at once. OpenAD comes loaded with a default workspace called 'default', and you can create as many additional workspaces as needed.


To see how to work with workspaces:

    ? workspace

## Plugins

Plugins are the way molecular tools and AI models are made available to the OpenAD client. Thanks to a unified language, accessing these tools through OpenAD lets you bypass a lot of complexity.

OpenAD comes preloaded with a number of plugins for literature knowledge extraction (DS4SD), forward and retrosynthesis prediction (RXN) as well as generative methods and property inference (GT4SD).

You can create your own plugins, and the publicly available plugins will soon include a much larger variety of open-source tools.

Note: Plugins are currently referred to as "toolkits" by the commands, however this language will be updated soon.

To see how to work with plugins:

    ? toolkit

## Context

Your context defines what plugin you are currently working with. It ensures you are logged onto the a plugin's corresponding system and all the related functionality is loaded.


To see how to switch contexts:

    ? context

## Runs

A run is a prerecorded chain of commands that can be reused to automate certain workflows and avoid unnecessary repetition.

To see how to work with runs:

    ? run
