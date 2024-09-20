# Working with Molecules

<svg fill="none" viewBox="0 0 600 300" width="600" height="300" xmlns="http://www.w3.org/2000/svg">
  <foreignObject width="100%" height="100%">
    <div xmlns="http://www.w3.org/1999/xhtml">
      <style>
        h1 {
            color: red
        }
      </style>

      <div class="container">
        <h1>Hi there</h1>
      </div>
    </div>

  </foreignObject>
</svg>

## Introduction

Working with molecules in OpenAD happens via the commands line interface (CLI) and the graphical user interface (GUI). The two interfaces serve different purposes and are meant to work together.

If you are working with Jupyter Notebook, you can access the CLI using Magic Commands, by prefixing commands with `%openad`. For example:

```sh
%openad show mol dopamine
```

To install the OpenAD Magic Commands in Jupyter, check [Launching OpenAD in Jupyter](../#launching-openad-in-jupyter) in the main readme file.

In our examples we'll be using the CLI from within a terminal.

## Available Molecule Commands

To see all available commands, run `?` and look for the ones listed under `Small Molecules`, `Macromolecules`, `Molecule Sets` and `Molecules Working Set`. To see the documentation for each individual command, you can consult the [online documentation](https://acceleratedscience.github.io/openad-docs/) or simply paste (the first part of) a command in your terminal followed by a `?`.

For example:

```sh
show mol|molecule <name> ?
```

...will display:

```sh
show mol|molecule <name> | <smiles> | <inchi> | <inchikey> | <cid>
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

## Visualizing & Storing Small Molecules

To see details about a small molecule like dopamine for example, simply run:

```sh
show molecule dopamine
```

This will launch the GUI. From there you can add the molecule to your working set or store it in your workspace.

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

## Visualizing & Storing Macromolecules

To see details about a macromolecule, simply run:

```sh
show molecule dopamine
```
