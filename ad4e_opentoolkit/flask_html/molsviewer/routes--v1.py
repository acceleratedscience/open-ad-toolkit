import os
import signal
from rdkit.Chem import PandasTools
from flask import render_template, send_from_directory, request
from ad4e_opentoolkit.app.global_var_lib import _repo_dir
from ad4e_opentoolkit.helpers.output import msg, output_text, output_error, output_warning, output_success, output_table
from ad4e_opentoolkit.helpers.spinner import spinner


def _compile_default_m2g_params(mol_frame):
    available_params = mol_frame.columns.tolist()
    if 'NAME' in available_params:
        m2g_params = {
            'subset': ['NAME'],
            'tooltip': [x for i, x in enumerate(available_params) if (x not in ['NAME', 'mols2grid-id', 'IMG'])]
        }
        return m2g_params


def _kill_server():
    os.kill(os.getpid(), signal.SIGINT)


def fetchRoutes(cmd_pointer, parser, the_mols2grid, mol_frame, available_params, m2g_params_default, create_missing_dirs, dir_path, results_file, workspace_path):

    # @app.route('/')
    def home():
        # Parse URL arguments.
        args = dict(request.args)
        if len(args) > 0:
            # Create mol2grid display parameters object.
            m2g_params = {k: v.split(',') for k, v in args.items()}
        else:
            # No arguments are passed, set default display parameters.
            # m2g_params = _compile_default_m2g_params(mol_frame) # %%
            m2g_params = m2g_params_default

        # Render the grid.
        # We use a copy of the m2g_params because the_mols2grid.display modifies it.
        import copy
        m2g_params_copy = copy.deepcopy(m2g_params)
        m2g_instance = the_mols2grid.display(**m2g_params_copy)
        return render_template('/molsviewer/index.html', data=m2g_instance.data, available_params=available_params, m2g_params=m2g_params)

    # Make page-specific CSS files available.
    # @app.route('/molsviewer/<path>')
    def send_page_css(path):
        return send_from_directory(_repo_dir + '/../flask_html/molsviewer', f'{path}')

    # Submit
    # @app.route('/submit', methods=['POST'])
    def submit():
        import json
        indexes = json.loads(request.data.decode('utf-8'))
        if indexes is not None:
            filtered_df = the_mols2grid.dataframe.iloc[indexes]

        if 'results_file' in parser:
            # Store file

            # Create the directories if thet don't exist.
            if create_missing_dirs:
                os.makedirs(dir_path)

            # Save as SDF file.
            if results_file.split('.')[-1].lower() == 'sdf':
                # Add the ROMol object to the dataframe.
                # Either by overriding the ROMol column, or by creating it.
                # Without this, no sdf file can be created.
                for index, row in filtered_df.iterrows():
                    if 'ROMol' not in filtered_df.columns:
                        filtered_df['ROMol'] = None
                    filtered_df.at[index, 'ROMol'] = rdkit.Chem.MolFromSmiles(row['SMILES'])

                PandasTools.WriteSDF(
                    filtered_df,
                    workspace_path + results_file,
                    properties=list(filtered_df.columns),
                )

            # Save as CSV file.
            elif results_file.split('.')[-1].lower() == 'csv':
                # We remove the RoMol object as it's useless in the CSV file.
                # It can be regenerated from the SMILES string.
                filtered_df = filtered_df.drop(['ROMol'], axis=1, errors='ignore')

                # We remove the img column as it's useless in the CSV file.
                filtered_df = filtered_df.drop(['IMG'], axis=1, errors='ignore')
                filtered_df.to_csv(
                    workspace_path + results_file,
                    index=False
                )

            # Success message
            spinner.stop()
            output_success(msg('success_m2g_save', len(indexes), parser["results_file"]), cmd_pointer)
        else:
            # Display results
            note = None  # 'To see what you can do next, run <cmd>next ?</cmd>'
            spinner.stop()
            output_success(msg('success_m2g_select', len(indexes)), cmd_pointer, pad_top=1)
            output_table(filtered_df, cmd_pointer, note=note)

        # Shut down server
        _kill_server()
        return 'ok'

    # @app.route('/exit', methods=['POST'])
    def exit():
        output_error(msg('abort'), cmd_pointer)
        _kill_server()
        return 'ok'
