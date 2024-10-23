# OpenAD Models Service

OpenAD lets you easily deploy different models to generate and manipulate molecule data sets.

<br>

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

-   <details><summary>GT4SD - Molformer inference</summary>

    <br>
    <div markdown="block">

        git@github.com:acceleratedscience/molformer_inference_service.git

    [GitHub](https://github.com/acceleratedscience/molformer_inference_service.git)

    Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus non tellus vel arcu porttitor tincidunt. Curabitur efficitur sodales efficitur. Ut id dui ut mi sodales tempor non nec erat. Pellentesque consectetur, nibh quis tempor luctus, quam lacus sagittis libero, id viverra lectus est id ligula. Nam sem ex, molestie a arcu finibus, dapibus eleifend felis. Interdum et malesuada fames ac ante ipsum primis in faucibus. Sed laoreet elit vestibulum porta viverra.

    </div>
    </details>

<br>

## Installation Requirements

1.  <details><summary>AWS account</summary>

    <div markdown="block">

    -   Head to [aws.com](https://aws.com/)
    -   Click the _[Create an AWS Account]_ button in the top right corner
    -   Follow instructions, including setting up a root user

    </div>
    </details>

2.  <details><summary>AWS user with correct permissions</summary>

    <div markdown="block">

    Starting from your [AWS dashboard]:

    -   Search for _"IAM"_ in the search bar
    -   From your [IAM dashboard], click _"[Users]"_ in the lefthand sidebar
    -   Click the _[Create user]_ button in the top right hand corner
    -   Leave the _"Provide user access to the AWS Management Console"_ box unchecked
    -   Up next on the _"Set Permissions"_ screen, select the third option: _"Attach policies directly"_
    -   In the box below, click the _[Create policy]_ button
    -   Create a new policy with minimal permissions for Skypilot, following thye [Skypilot instructions](https://skypilot.readthedocs.io/en/latest/cloud-setup/cloud-permissions/aws.html)
    -   On the next screen, search for the policy you just created, which would be called `minimal-skypilot-policy` per the instructions
    -   Finish the process to attach the policy to your user

    </div>
    </details>

3.  <details><summary>AWS Access key</summary>

    <div markdown="block">

    Starting from the [IAM dashboard]:

    -   Click _"[Users]"_ in the lefthand sidebar
    -   Click on the user you created in the previous step
    -   Click _"Create access key"_ on the right side of the summary on top
    -   Select the first option, _"Command Line Interface (CLI)"_ as use case
    -   Finish the process to create the access key
    -   Store the secret access key in your password manager, as you will not be able to access it after creation

    </div>
    </details>

4.  <details><summary>AWS command line tool</summary>

    <div markdown="block">

    Starting from a terminal window:

    -   Install `awscli`

            python -m pip install awscli

        > **Note:** For more nuanced instructions, please refer to [pypi](https://pypi.org/project/awscli/#getting-started)

    -   Add the credentials for the AWS user you set up in step 3.

            aws configure

        -   Your user's access key can be found in your [IAM dashboard] > [Users], however the secret access key should have been stored in your password manager or elsewhere.
        -   The fields _"Default region name"_ and _"Default output format"_ can be left blank

    </div>
    </details>

5.  <details><summary>SkyPilot</summary>

    <div markdown="block">

    [SkyPilot](https://skypilot.readthedocs.io/en/latest/getting-started/installation.html) is a framework for running AI and batch workloads on any infrastructure. We're using AWS.

    Starting from a terminal window:

    -   If you are running OpenAD in a virtual environment, make sure your virtual environment is activated. If you followed the [default installation instructions](../#installation), you should be able to run:

            source ~/ad-venv/bin/activate

    -   Install Skypilot for AWS

            pip install "skypilot[aws]"

    -   After installation, verify if you have cloud access

            sky check

    </div>
    </details>

<br>

## Installation

1.  Install any service by its url (see [available models](#available-models) on top), for example the property inference service:

        catalog model service from 'git@github.com:acceleratedscience/property_inference_service.git' as prop

1.  Start the service – this can take up to 10 minutes

        model service up prop

1.  Check if the service is ready

        model service status

1.  Shut down the service

        model service down prop

1.  To see all available model commands, pull up the general help and look towards the bottom of the command list.

        ?

[AWS Dashboard]: https://console.aws.amazon.com
[IAM Dashboard]: https://console.aws.amazon.com/iam
[Users]: https://console.aws.amazon.com/iam/home#/users
