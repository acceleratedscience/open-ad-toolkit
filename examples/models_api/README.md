# Beginners Guide to using OpenAD model inference API

⚠️**Current Available Public Models. Under Construction** ⚠️

| Molformer                 | Moler         | Properties| Generation|
| --------                  | --------      | --------  | --------  |
| regression                |               |           |           |
| multitask classification  |               |           |           |
| classifciation            |               |           |           |


## What does this service do?

This service offered by Accelerated Discovery at IBM enables researchers to run model inference as a service (MIaaS). Public models are available to use by all by connecting to our service with an api key. Private models can also be hosted by our service to enable users to run inference without setting up machine learning models locally.

Our service helps users:
- Eliminate the need to run inference locally
- Easy deployment of ML models as a service
- Integrate an ML model to be usable directly in the Openad Toolkit CLI and notebooks
- Not worry about infrastructure just focus on making ML models


## Getting Started

### Create Account
To run inference with our models via the api you first need to create an account [[here](https://open.accelerator.cafe/)](https://open.accelerator.cafe/)

If you have any issues or inquiries please reach out to us via [phil.downey1@ibm.com](mailto:phil.downey1@ibm.com)

### Generate token
Upon account creation you will have access to the default publicly available groups. Now you need to get your access token to use the service. Once generated copy and proceed.

![alt text](/assets/proxy/access_token.png)

### Connect openad to api

Start up openad toolkit cli. Example connecting to the `molformer` model:
```shell
OpenAD:DEFAULT >> catalog model service from remote 'https://open.accelerator.cafe/proxy' as 'molformer' USING (Inference-Service=molformer)
```

Lets break this command down
- `catalog model service from remote` *command in openad that processes the connection*
- `https://open.accelerator.cafe/proxy` *URL endpoint for inference*
- `'molformer'` *name you want to give this service (could be anything you want)*
- `USING (Inference-Service=molformer)` *select the model you want to interface with (check your dashboard for available models)*

Create an authentication group for your api and pass the access token to the service to authenticate yourself. Replace `<api_key>` with your token from the dashboard.
```shell
OpenAD:DEFAULT >> model auth add group admin with  '<access_token>'

OpenAD:DEFAULT >> model auth add service molformer to group admin
```

The service should show as Connected
```shell
OpenAD:DEFAULT >>  model service status

Service    Status     Endpoint                             Host    Token Expires
---------  ---------  -----------------------------------  ------  --------------------------
molformer  Connected  https://open.accelerator.cafe/proxy  remote  Wed Sep  11, 2030
```

Take a peak at the whats available in the model
```shell
OpenAD:DEFAULT >>  molformer ?

Commands starting with "molformer"
- molformer get molecule property molformer_classification for [<list of SMILES>] | <SMILES>   USING (<parameter>=<value> <parameter>=<value>) (save_as '<filename.csv>')
- molformer get molecule property molformer_multitask_classification for [<list of SMILES>] | <SMILES>   USING (<parameter>=<value> <parameter>=<value>) (save_as '<filename.csv>')
- molformer get molecule property molformer_regression for [<list of SMILES>] | <SMILES>   USING (<parameter>=<value> <parameter>=<value>) (save_as '<filename.csv>')
```

### Inference in CLI

Run the following command to get a classification result
```shell
OpenAD:DEFAULT >>  molformer get molecule property molformer_classification for 'OC12COC3=NCC1C23'
✔ Request Returned

subject           property                  result
----------------  ------------------------  --------
OC12COC3=NCC1C23  molformer_classification  [1]
```
