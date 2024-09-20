# Working with molecules

## Introduction

Working with molecules in OpenAD happens via the commands line interface (CLI) and the graphical user interface (GUI). The two interfaces serve different purposes and are meant to work together.

If you are working with Jupyter Notebook, you can access the CLI using Magic Commands, by prefixing commands with `%openad`. For example:

```sh
%openad show mol dopamine
```

To get started with the OpenAD Magic Commands in Jupyter, please consult [Launching OpenAD in Jupyter](../#launching-openad-in-jupyter).

In our examples we'll be using the CLI from within a terminal.

<br>

## Available molecule commands

To see all available commands, run `?` and look for the ones listed under `Small Molecules`, `Macromolecules`, `Molecule Sets` and `Molecules Working Set`. To see the documentation for each individual command, you can consult the [online documentation](https://acceleratedscience.github.io/openad-docs/) or simply paste (the first part of) a command in your terminal followed by a `?`.

For example:

```sh
show mol|molecule <name> ?
```

...will display:

```text
|    show mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>
|    ------------------------------------------------------------------
|
|    Inspect a molecule in the browser. If a molecule is not in the current Molecule Working set it will pull the result from Pubchem.
|
|    You can use the 'mol' shorthand instead of 'molecule'.
|
|    When you show a molecule by SMILES or InChI, we can display it immediately. When you show a molecule by name, InChIKey or PubChem CID, we
|    need to first retrieve it from PubChem, which can take a few seconds.
|
|    Examples:
|    - show mol aspirin
|    - show mol CC(=O)OC1=CC=CC=C1C(=O)O
|    - show mol InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)
|    - show mol 2244
```

<br>

## Visualizing & storing small molecules

To see details about a small molecule like eg. dopamine, you can run:

```sh
show molecule dopamine
```

This will launch the GUI. From there you can choose to add the molecule to your working set or store it in your workspace.

To see what files exist in your workspace, you can either run:

```sh
list files
```

Or you can launch the GUI to browse your workspace in a more user friendly manner:

```sh
launch gui
```

Our molecule viewer can open any `SDF`, `MOL` and `SMI` file, however we recommend your store individual small molecules as `.smol.json` files in your workspace, or grouped into a molecule set using the `.molset.json` format. How to do this is self explanatory when you are using the GUI.

<br>

## Visualizing & storing macromolecules

To see details about a macromolecule like eg. [2g64](https://www.rcsb.org/structure/2g64), you can run:

```sh
show mmol '2g64'
```

This will launch the GUI where can choose to store the molecule in your workspace.

> **Note:** For the time being, you can't store macromolecules into your molecules working set, nor can you add them to a molset.
