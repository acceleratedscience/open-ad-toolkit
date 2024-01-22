# Example commands:
# interpret recipe 'A solution of ((1S,2S)-1-{[(methoxymethyl-biphenyl-4-yl)-(2-pyridin-2-yl-cyclopropanecarbonyl)-amino]-methyl}-2-methyl-butyl)-carbamic acid tert-butyl ester (25 mg, 0.045 mmol) and dichloromethane (4 mL) was treated with a solution of HCl in dioxane (4 N, 0.5 mL) and the resulting reaction mixture was maintained at room temperature for 12 h. The reaction was then concentrated to dryness to afford (1R,2R)-2-pyridin-2-yl-cyclopropanecarboxylic acid ((2S,3S)-2-amino-3-methylpentyl)-(methoxymethyl-biphenyl-4-yl)-amide (18 mg, 95% yield) as a white solid.'
# interpret recipe 'recipe.txt'


import os
import importlib.util as ilu
from openad.helpers.output import output_text, output_error
from openad.helpers.general import remove_lines


def interpret_recipe(inputs: dict, cmd_pointer):
    """
    Interpret a free text paragraph into a recipe list of instructions.
    """

    recipe = inputs["recipe"]
    rxn4chemistry_wrapper = cmd_pointer.login_settings["client"][cmd_pointer.login_settings["toolkits"].index("RXN")]

    # Display processing...
    if (
        os.path.isfile(cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + recipe.strip())
        is True
    ):
        output_text(f"<soft>Processing the <yellow>{recipe}</yellow> file...</soft>", return_val=False, pad_top=1)
        with open(
            cmd_pointer.workspace_path(cmd_pointer.settings["workspace"].upper()) + "/" + recipe.strip(),
            "r",
            encoding="utf-8",
        ) as handle:
            recipe = handle.read()
    else:
        output_text(f"<soft>Processing paragraph...</soft>", return_val=False, pad_top=1)

    # Print result
    try:
        # raise Exception('This is a test error')
        return_result = []
        actios_from_procedure_results = rxn4chemistry_wrapper.paragraph_to_actions(recipe)
        return_result.append("<bold>Recipe steps:</bold>")
        for index, action in enumerate(actios_from_procedure_results["actions"], 1):
            return_result.append(f"{index}. {action}")
        remove_lines(2)
        return output_text("\n".join(return_result), pad=1, nowrap=True)
    except Exception as err:
        output_error(["Unable to parse the provided paragraph", err], return_val=False)
