# OpenAD Models Service

OpenAD lets you easily deploy different models to generate and manipulate molecule data sets.

## Available Models

<!--  -->

-   <details><summary>GT4SD - Generation inference</summary>

    <br>
    <div markdown="block">

        git@github.com:acceleratedscience/generation_inference_service.git

    [GitHub](https://github.com/acceleratedscience/generation_inference_service.git)

    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus non tellus vel arcu porttitor tincidunt. Curabitur efficitur sodales efficitur. Ut id dui ut mi sodales tempor non nec erat. Pellentesque consectetur, nibh quis tempor luctus, quam lacus sagittis libero, id viverra lectus est id ligula. Nam sem ex, molestie a arcu finibus, dapibus eleifend felis. Interdum et malesuada fames ac ante ipsum primis in faucibus. Sed laoreet elit vestibulum porta viverra.

    </div>
    </details>

<!--  -->

-   <details><summary>GT4SD - Property inference</summary>

    <br>
    <div markdown="block">

        git@github.com:acceleratedscience/property_inference_service.git

    [GitHub](https://github.com/acceleratedscience/property_inference_service.git)

    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus non tellus vel arcu porttitor tincidunt. Curabitur efficitur sodales efficitur. Ut id dui ut mi sodales tempor non nec erat. Pellentesque consectetur, nibh quis tempor luctus, quam lacus sagittis libero, id viverra lectus est id ligula. Nam sem ex, molestie a arcu finibus, dapibus eleifend felis. Interdum et malesuada fames ac ante ipsum primis in faucibus. Sed laoreet elit vestibulum porta viverra.

    </div>
    </details>

<!--  -->

-   <details><summary>GT4SD - MoleR inference</summary>

    <br>
    <div markdown="block">

        git@github.com:acceleratedscience/moler_inference_service.git

    [GitHub](https://github.com/acceleratedscience/moler_inference_service.git)

    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus non tellus vel arcu porttitor tincidunt. Curabitur efficitur sodales efficitur. Ut id dui ut mi sodales tempor non nec erat. Pellentesque consectetur, nibh quis tempor luctus, quam lacus sagittis libero, id viverra lectus est id ligula. Nam sem ex, molestie a arcu finibus, dapibus eleifend felis. Interdum et malesuada fames ac ante ipsum primis in faucibus. Sed laoreet elit vestibulum porta viverra.

    </div>
    </details>

<!--  -->

-   <details><summary>GT4SD - MoleR inference</summary>

    <br>
    <div markdown="block">

        git@github.com:acceleratedscience/molformer_inference_service.git

    [GitHub](https://github.com/acceleratedscience/molformer_inference_service.git)

    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus non tellus vel arcu porttitor tincidunt. Curabitur efficitur sodales efficitur. Ut id dui ut mi sodales tempor non nec erat. Pellentesque consectetur, nibh quis tempor luctus, quam lacus sagittis libero, id viverra lectus est id ligula. Nam sem ex, molestie a arcu finibus, dapibus eleifend felis. Interdum et malesuada fames ac ante ipsum primis in faucibus. Sed laoreet elit vestibulum porta viverra.

    </div>
    </details>

## Installation

### Requirements

1. <details><summary>AWS account</summary>

    <div markdown="block">

    - Head to [aws.com](https://aws.com/)
    - Click the **[Create an AWS Account]** button in the top right corner
    - Follow instructions, including setting up a root user

    </div>
    </details>

1. <details><summary>AWS user with correct permissions</summary>

    <div markdown="block">

    Starting from your AWS dashboard:

    - Search for "IAM" in the search bar
    - From your IAM dashboard, click "Users" in the lefthand sidebar
    - Click the "Create user" button in the top right hand corner
    - Leave the "Provide user access to the AWS Management Console" box unchecked
    - Up next on the "Set Permissions" screen, select the third option: "Attach policies directly"
    - In the box below click the "Create policy" button
    - Create a new policy with minimal permissions for Skypilot, following thye [Skypilot instructions](https://skypilot.readthedocs.io/en/latest/cloud-setup/cloud-permissions/aws.html)
    - On the next screen, search for the policy you just created, which would be called "minimal-skypilot-policy" per the instructions
    - Finish the process to attach the policy to your user

    - Go to user detail, click "create access key" in summary
    - Select CLI first option

    </div>
    </details>

1. <details><summary>AWS Access key</summary>

    <div markdown="block">

    Starting from the IAM dashboard:

    - Click "Users" in the lefthand sidebar
    - Click on the user you created in the previous step
    - Click "create access key" on the right side of the summary on top
    - Select the first option, "Command Line Interface (CLI)" as use case
    - Finish the process to create the access key

    </div>
    </details>

1. <details><summary>Conda</summary>

    <div markdown="block">

    Conda is a powerful command line tool for package and environment management. It is a prerequisite to work with SkyPilot.

    > **Note:** Conda provides [multiple installers](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). If you're just starting out, we recommend you go with Miniconda.

    - Install [Miniconda](https://docs.anaconda.com/miniconda/#quick-command-line-install) using the quick command line install

        > **Note:** After the installation, the installation script will ask you if you wish to activate Conda on start. You may not want to do this. If you choose to disable this, you can activate conda at any time like this:

        - Go to the miniconda installation folder, this is propbably

            cd ~/miniconda3

        - Activate the conda environment

            source ./bin/activate

    - After installing miniconda or an alternative, make sure to open a _new_ terminal window to access Conda functionality
    - To verify that conda is installed, run `conda --version`

    </div>
    </details>

1. <details><summary>SkyPilot</summary>

    <div markdown="block">

    SkyPilot is a framework for running AI and batch workloads on any infrastructure. We're using AWS.

    - https://skypilot.readthedocs.io/en/latest/getting-started/installation.html

    </div>
    </details>

1. <details><summary>AWS command line tool</summary>

    <div markdown="block">

    In a terminal window type:

    - `python -m pip install awscli`
    - For more nuanced instructions, please refer to [pypi](https://pypi.org/project/awscli/#getting-started)

    </div>
    </details>
