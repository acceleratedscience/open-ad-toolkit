import json
from flask import render_template


def fetchRoutesMolViewer(mol_dict, mol_sdf, svg):
    mol_json = json.dumps(mol_dict, indent="\t")

    def main():
        return render_template("/molviewer/index.html", mol=mol_dict, mol_json=mol_json, mol_sdf=mol_sdf, svg=svg)

    routes = {
        "/": {"func": main},
    }

    return routes
