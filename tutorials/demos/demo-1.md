OpenAD Demonstration

# Source PFAS molecules and find alternatives

In this demonstration we will:

-   Search for PFAS molecules on PubChem
-   Generate additional properties
-   Collate data between molecule sets
-   Generate similar molecules with higher soluability using IBM's open-source Regression Transformer
-   Use Deep Search to determine if any of our generated molecules are mentioned in a patent, and only proceed with the ones that do not
-   Take one of the molecules and use the RXN toolkits's retrosynthesis commands to generate a path to synthesis

# Step 1: Search for PFAS molecules on PubChem

We'll use the DS4SD (Deep Search) toolkit to crawl the PubChem molecule database.

First we need to install the toolkit. If it was previously installed, you will need to run `set context ds4sd` instead.

```sh
add toolkit ds4sd
```

Next we run our query for PFAS molecules.

```sh
search collection 'PubChem' for 'PFOA OR PFOS OR PFHxS OR PFNA OR HFPO-DA'
```

Now we can inspect the results in the GUI and save the molecules we want to our working set. This is done by clicking the bookmark icon on each individual molecule.

```sh
result open
```

# Step 2: Generate additional properties

With our candidate molecules loaded into the working set, we can use OpenAD's built-in property generation to generate missing properties.

# Step 3: Collate data between molecule sets
