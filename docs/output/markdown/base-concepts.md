---
title: Base Concepts
layout: default
nav_order: 4
---

<!--

DO NOT EDIT
-----------
This file is auto-generated.
To update it, consult instructions:
https://github.com/acceleratedscience/open-ad-toolkit/tree/main/docs

-->

# Base Concepts

## Workspaces

Your workspace represents an isolated environment where your molecules, molecule sets and other files related to your research are stored.

Each workspace has its own command history and corresponds with a directory stored in ~/openad/workspaces, allowing you to work on multiple isolated projects at once. OpenAD comes loaded with a default workspace called 'default', and you can create as many additional workspaces as needed.

To see how to work with workspaces:

    ? workspace

## Molecule Working Set

Your molecule working set (MWS) is an in-memory list of small molecules that is the subject of your calculations and manipulations.

The MWS makes it easy to gather candidate molecules from various sources, either by adding them individually, by loading them in batch from a dataframe, CSV or SDF file, or by bookmarking molecules in the GUI. You can then calculate certain properties, perform triage and then store the results in your workspace as a molecule set, ready for further processing within OpenAD or elsewhere.

Please note that at this time the MWS does not support macromolecules, only small molecules.

To see how to work with your molecule working set:

    ? mws

## Plugins

Plugins are how molecular tools and AI models are made available to the OpenAD client. They provide drastically simplified access to a series of advanced tools, and they make it easy for your own Python applications to interface with OpenAD.

Creating your own plugins is easy if you have a basic understanding of Python.

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
