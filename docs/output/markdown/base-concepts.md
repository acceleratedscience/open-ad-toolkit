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

You can have as many workspaces as needed, allowing you to work on multiple isolated projects. Each workspace has its own command history and corresponds with a directory stored in ~/openad/workspaces. OpenAD comes loaded with a default workspace called 'default'.

To see how to work with workspaces:

    ? workspace

## Plugins

The OpenAD client is an interface to interact with a variety of molecular tools and AI models, which are exposed through plugins. Thanks to a unified language, accessing these tools through OpenAD lets you bypass a lot of complexity.

OpenAD comes preloaded with a number of plugins for literature knowledge extraction (DS4SD), forward and retrosynthesis prediction (RXN) as well as generative methods and property inference (GT4SD).

You can create your own plugins, and the publicly available plugins will soon include a much larger variety of open-source tools.

Note: Plugins are currently referred to as "toolkits" by the commands, however this language will be updated soon.

To see how to work with plugins:

    ? toolkit

## Context

In order to interact with any plugin, you first need to set the context to that plugin. This will make all commands for the plugin available.

To see how to switch contexts:

    ? context

## Runs

A run is a prerecorded chain of commands that can be reused to automate certain workflows and avoid unnecessary repetition.

To see how to work with runs:

    ? run
