from setuptools import setup, find_packages, find_namespace_packages
import os


""" packages=['main', 'main.core'],
    package_dir={
        'project': '.',
        'project.submodule': './submodule1',
    }, """

setup(
    name="openad",
    version="0.1.0",
    # packages=find_packages(),
    packages=find_namespace_packages(include=["openad.*"]),
    include_package_data=True,
    install_requires=[
        "joblib",
        "matplotlib",
        "scikit_learn",
        "blessed",
        "flaky==3.7.0",
        "Flask==2.3.2",
        "Flask_Bootstrap==3.3.7.1",
        "geckodriver_autoinstaller==0.1.0",
        "halo==0.0.31",
        "imagehash==4.3.1",
        "ipython==8.12.0",
        "ipywidgets==7.7.5",
        "Jinja2==3.1.2",
        "nodejs==0.1.1",
        "npm",
        "numerize==0.12",
        "numpy==1.24.2",
        "pandas==2.0.0",
        "protobuf==3.20.3",
        "PubChemPy==1.0.4",
        "flask==2.3.2",
        "py3Dmol",
        "PyJWT==2.8.0",
        "pyparsing==3.0.9",
        "pytest==7.3.1",
        "rdkit==2023.3.2",
        "rdkit_pypi==2022.9.5",
        "streamlit==1.22.0",
        "tqdm==4.65.0",
        "typing",
        "traitlets",
        "pandas",
        "mols2grid",
        "deepsearch-toolkit==0.29.0",
        "tiktoken",
        "langchain==0.0.325",
        "traitlets==5.9",
        "deepsearch-toolkit",
        "rxn4chemistry",
        "tabulate",
        "faiss-cpu",
        "openai",
        "IPython==8.12.0",
        "ipykernel==6.22.0",
        "ipywidgets==7.7.5",
        "jupyter_client==8.2.0",
        "jupyter_core==5.3.0",
        "jupyter_server==2.7.2",
        "jupyterlab==3.5.3",
        "nbclient==0.7.3",
        "nbconvert==7.3.1",
        "nbformat==5.8.0",
        "notebook==6.5.4",
        "traitlets==5.9.0",
        "ipykernel",
        "pyperclip>=1.8.2",
    ],
    entry_points={
        "console_scripts": [
            "openad=openad.app.main:cmd_line",
            "init_magic=openad.app.init_magic:init_magic",
            "init_examples=openad.app.openad_examples:openad_create_examps",
        ],
    },
)
