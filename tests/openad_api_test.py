from openad import OpenadAPI


if __name__ == "__main__":

    myclass2 = OpenadAPI("class2")

    """print(myclass2.name)
    x = myclass2.request("prop get molecule property lipinski for CCO")
    print(x)
    properties_all = [
        "molecular_weight",
        "number_of_aromatic_rings",
        "number_of_h_acceptors",
        "number_of_atoms",
        "number_of_rings",
        "number_of_rotatable_bonds",
        "number_of_large_rings",
        "number_of_heterocycles",
        "number_of_stereocenters",
        "is_scaffold",
        "bertz",
        "tpsa",
        "logp",
        "qed",
        "plogp",
        "penalized_logp",
        "lipinski",
        "sas",
        "esol",
    ]

    a_molecule_list = [
        "O=C(O)C(F)(OC(O)(F)C(F)(F)C(F)(F)F)C(F)(F)F",
        "ON(O)C(F)(OC(F)(F)C(F)(F)C(F)(F)F)C(F)(F)F",
        "C(C(C1C(=C(C(=O)O1)O)O)O)O",
    ]
    for x in a_molecule_list:
        x = myclass2.request(f"add molecule {x} basic force")
    x = myclass2.request(f"prop get molecule property {properties_all} for  {a_molecule_list}")
    x = myclass2.request(f"merge molecules data using dataframe properties")
    x = myclass2.request("export molecules")
    print(x)
    x = myclass.request("export molecules")
    print(x)

    # print(myclass2.help_dump())

    # print(myclass2.request("? similar"))

    # print(myclass2.help_as_markdown("? similar"))

    # print(myclass2.request("export molecule CCO"))

    myclass2.request("set context DS4SD")

    df = myclass2.request("search collection 'PubChem' for 'PFOA OR PFOS OR PFHxS OR PFNA OR HFPO-DA'")

    myclass2.request("load molecules using dataframe df", **vars())
    a_list = list(set(df["SMILES"].to_list()))
    # Define list of Delta to be inferred properties
    properties = ["is_scaffold", "bertz", "tpsa", "logp", "qed", "plogp", "penalized_logp", "lipinski", "sas", "esol"]

    # Generate SMILES properties
    list_of_properties = myclass2.request(f" prop get molecule property {properties} for  {a_list}")
    myclass2.request(f"merge molecules data using dataframe list_of_properties", list_of_properties=list_of_properties)
    mol_list = myclass2.request("export molecules")

    datasets = []
    for row in mol_list.to_dict("records"):
        MY_SMILES = row["canonical_smiles"]
        esol = float(row["esol"]) + 2
        MY_PARAMS = {"fraction_to_mask": 0.1, "property_goal": {"<esol>": esol}}
        print("Generating Molecules for " + MY_SMILES + " with soluability:" + str(row["esol"]))
        command = f'gen generate with RegressionTransformerMolecules data for {MY_SMILES} sample 10 using(algorithm_version=solubility  search=sample temperature=1.5 tolerance=60.0 sampling_wrapper = "{str(MY_PARAMS)}" '

        result = myclass2.request(
            f'gen generate with RegressionTransformerMolecules data for {MY_SMILES} sample 10 using(algorithm_version=solubility  search=sample temperature=1.5 tolerance=60.0 sampling_wrapper = "{str(MY_PARAMS)}") '
        )
        datasets.append(result)

    x = 0
    from pandas import DataFrame

    patent_count = 0
    patents_to_search = []
    patented_molecules = []
    non_patented_molecules = []
    for result in datasets:
        for mol in result["0"].to_list():
            x = myclass2.request(f"search for patents containing molecule '{mol}'")
            if isinstance(x, DataFrame):
                patents_to_search.extend(x["PATENT ID"].to_list())
                patented_molecules.append(mol)
            else:
                if mol not in non_patented_molecules:
                    non_patented_molecules.append(mol)

    properties_all = [
        "molecular_weight",
        "number_of_aromatic_rings",
        "number_of_h_acceptors",
        "number_of_atoms",
        "number_of_rings",
        "number_of_rotatable_bonds",
        "number_of_large_rings",
        "number_of_heterocycles",
        "number_of_stereocenters",
        "is_scaffold",
        "bertz",
        "tpsa",
        "logp",
        "qed",
        "plogp",
        "penalized_logp",
        "lipinski",
        "sas",
        "esol",
    ]
    print(1)
    new_props = myclass2.request(" prop get molecule property {properties_all} for {patented_molecules} ")
    print(2)
    for x in patented_molecules:
        myclass2.request(f" add molecule {x} Force")
    print(3)
    myclass2.request(" merge molecules data using dataframe new_props", new_props=new_props)
    print(4)
    myclass2.request("enrich molecules with analysis")
    print(5)"""
    myclass2.request(" set context rxn ")
    print(6)
    print(
        myclass2.request(" predict retrosynthesis  'O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F' ")
    )

    df = myclass2.request(
        "bmfm generate with BmfmAntibodiesGenerator data sample 5 using ( algorithm_type=generation protein_id=5Nmv heavy_chain_id=H  with_properties=True  )"
    )
    print(df)
