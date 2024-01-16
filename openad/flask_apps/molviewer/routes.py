import json
from flask import render_template


def fetchRoutesMolViewer(mol, mol_sdf, mol_svg):
    mol_json = json.dumps(mol, indent="\t")

    def main():
        return render_template("/molviewer/index.html", mol=mol, mol_json=mol_json, mol_sdf=mol_sdf, mol_svg=mol_svg)

    routes = {
        "/": {"func": main},
    }

    return routes
