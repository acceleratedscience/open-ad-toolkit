# AI Assistant

To enable our AI assistant, you'll either need access to [IBM BAM](https://bam.res.ibm.com/auth/signin) or use the free open source LLM [Ollama](https://ollama.com).

> **Note:** Ollama requires an 8GB GPU
> **Note:** Support for WatsonX is coming soon

<br>

## IBM BAM Setup

To use IBM BAM if you have access to it, simply provide your API key when prompted.

    set llm bam
    tell me <enter prompt>

<br>

## Ollama Setup

1.  Install [Ollama](https://ollama.com/download) onto your platform.

1.  Download the appropriate models.

        ollama pull llama3:latest
        ollama pull nomic-embed-text

1.  Start the server if not already started.

        ollama serve

That's it for local usage. If you want to run Ollama remotely, continue below.

### Ollama Remote Setup with SkyPilot

1.  Check out our configuration file to launch ollama on SkyPilot: [ollama_setup.yaml](./openad/ollama_setup.yaml)

        sky serve up ollama_setup.yaml

1.  Set up local environment variables

    -   For windows `setx OLLAMA_HOST=<sky-server-ip>:11434`
    -   For Linux and macOS `export OLLAMA_HOST=<sky-server-ip>:11434`
    -   To reset to local use `OLLAMA_HOST=0.0.0.0:11434`

### Run Ollama

> **Note:** If prompted for an API key and none was setup, just leave the input empty.

    set llm ollama
    tell me <enter prompt>

<br>