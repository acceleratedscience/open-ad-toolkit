import json
from flask import render_template, request
from openad.molecules.mol_api import get_molecule_data


def fetchRoutesMolViewer(cmd_pointer, mol, mol_sdf, mol_svg):
    from openad.molecules.mol_functions import molformat_v2

    mol = molformat_v2(mol)
    mol_json = json.dumps(mol, indent="\t")

    def main():
        _stringify_prop_sources(mol)
        return render_template(
            "/molviewer/index.html",
            mol=mol,
            mol_json=mol_json,
            mol_sdf=mol_sdf,
            mol_svg=mol_svg,
        )

    # API endpoint to get the RDKit-enriched data.
    def enrich():
        inchi = request.data.decode("utf-8")
        mol = get_molecule_data(cmd_pointer, inchi)
        if mol:
            mol = molformat_v2(mol)
            mol_json = json.dumps(mol, indent="\t")
            _stringify_prop_sources(mol)
            html = render_template(
                "/molviewer/index.html",
                mol=mol,
                mol_json=mol_json,
                mol_sdf=mol_sdf,
                mol_svg=mol_svg,
            )
            return html
        else:
            return "", 500

    routes = {
        "/": {"func": main},
        "/enrich": {"func": enrich, "method": "POST"},
    }

    return routes


def _stringify_prop_sources(mol):
    # Turn the molecule's property_sources into strings that we can render.
    for prop in mol["property_sources"]:
        src_str = []
        if mol["property_sources"][prop] and type(mol["property_sources"][prop]) is dict:
            for key, val in mol["property_sources"][prop].items():
                src_str.append(f"{key}: {val}")
        mol["property_sources"][prop] = "\n".join(src_str)
