---
title: "Running nf-core Pipelines on Google Colab"
subtitle: "A guide to running and interacting with nf-core pipelines using Colab and VS Code"
---

## Running nf-core Pipelines in Google Colab

Running nf-core pipelines can be computationally intensive, requiring resources not easily available to students, newcomers, or participants in hands-on training workshops.
This tutorial shows two ways you can use Google Colab to run nf-core pipelines entirely in the cloud, either in a browser-based execution or connected to VS Code via the `vscode-colab` package for a full development experience.
While Colab has limitations, such as session timeouts and a lack of root access, it offers a free and accessible platform ideal for learning, teaching, and prototyping workflows in resource-constrained environments.
Also make sure to checkout the [blogpost](https://nf-co.re/blog/2025/nf-core-colab-guide) on this topic for deeper insights and more tips for running nf-core pipelines using Colab.

## Setting up the environment in Google Colab

If you are new to Google Colab, you can follow the official [Google Colab Getting Started guide](https://research.google.com/colaboratory/faq.html) or the [Colab Welcome Notebook](https://colab.research.google.com/notebooks/intro.ipynb) for instructions on how to register for an account and set up your first notebook.

After opening a Colab notebook with the machine type of your choice, the first step is to install Nextflow and its dependencies. Since Colab does not provide root access, you will need to install Nextflow in your home directory.

### Setting up Java for Nextflow

Before installing Nextflow, you must install Java. In Colab, run the following commands in a notebook code cell using the `%` and `!` prefixes, or run them directly in your terminal (without the prefixes).

```python
%cd .. # use '%'before cd to make a permanent change of directories
!sudo apt update
!sudo apt install openjdk-17-jdk
!export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
!export PATH=$JAVA_HOME/bin:$PATH
!source ~/.bashrc
```

### Installing and Setting up Nextflow

Next, install Nextflow and make it executable by running the following commands in a code cell as before:

```python
!wget -qO- https://get.nextflow.io | bash # Download Nextflow
!mv nextflow /usr/bin/nextflow # Move to a path Colab can access
!chmod +x /usr/bin/nextflow # Make it executable
!nextflow -v # Test it
```

### Setting Up Conda for Google Colab

Due to a lack of root access, it is not possible to run nf-core, or any Nextflow pipelines, using the `-profile docker` or `-profile singularity` container based configuration profiles.

That leaves the conda profile as the only way to run pipelines.
However, as Google Colab does not support native conda functionality, you need to install the [condacolab](https://pypi.org/project/condacolab/) package, which will set up Conda for you by running the following commands:

```python
!pip install -q condacolab
import condacolab
condacolab.install()
```

Make sure to set the conda defaults as follows for smooth pipeline execution:

```python
!conda config --add channels defaults
!conda config --add channels bioconda
!conda config --add channels conda-forge
!conda config --set channel_priority strict
```

## Running a test pipeline

Yes, the setup is that simple!
Now that the main parts of the environment are all set up, we can run our demo pipeline [nf-core/demo](https://nf-co.re/demo/).
To install the pipeline, run:

```python
! nextflow pull nf-core/demo
```

One last important step before running your pipeline, or any others in Colab, is to set the `MPLBACKEND` environment variable to `Agg`.
While you might not expect to use Matplotlib directly in nf-core pipelines, some dependencies or scripts may use it behind the scenes for rendering or saving plots, which can cause crashes if not configured for a headless environment like Colab.

You can do this either by running the following in a code cell:

```python
%env MPLBACKEND=Agg
```

Or alternatively, by running the following command in the terminal:

```bash title="Set MPLBACKEND to Agg in the terminal"
export MPLBACKEND=Agg
```

Now you can finally run your pipeline!

```python
! nextflow run nf-core/demo -profile conda,test --outdir demo-results
```

## Running and Editing Pipelines in VS Code via Colab

While you could get away with editing existing pipelines inside Colab's built-in terminal using editors like vim or nano, a VSCode IDE offers a more robust and richer environment for development.
Thankfully, the [vscode-colab](https://github.com/EssenceSentry/vscode-colab) Python library provides just the toolkit you need to take advantage of Colab's hardware in the comfort of the popular VSCode software suite.

The library makes use of the official [VS Code Remote Tunnels](https://code.visualstudio.com/docs/remote/tunnels) to securely and reliably connect Google Colab as well as Kaggle notebooks to a local or browser-based instance of VS Code.
You can read more about the library and even help contribute to new features on its [GitHub repository](https://github.com/EssenceSentry/vscode-colab).

:::note
If you decide to use this approach, it is best to run this after setting up the other dependencies that require installation by running inside Colab cells to avoid any connection issues.
:::

The first step is to install and import the library, which can easily be done by running the following command in your Colab code cell:

```python title="Install vscode-colab"
!pip install vscode-colab
import vscode_colab
```

Next, authenticate the connection using your GitHub credentials:

```python
vscode_colab.login()
```

![interactive login popup in cell output](/images/tutorials/google_colab/login.png)

Follow the displayed instructions to authorize the connection.

To start the VS Code tunnel, optionally configure Git, set up a Python version, or create a new project:

```python
vscode_colab.connect(
    name="my-tunnel",
    git_user_name="Your Name",
    git_user_email="you@example.com",
    setup_python_version="3.13",  # Optional: Specify Python version to install with pyenv
    create_new_project="my_new_project"  # Optional: Create a new project directory
)
```

![interactive tunnel connection popup in cell output](/images/tutorials/google_colab/connect.png)

The above step allows for a lot more customization than shown here. Please refer to the library's [vignette](https://github.com/EssenceSentry/vscode-colab) to see if any of these options will suit your needs.

Finally, in your local VS Code:

1. Ensure the [Remote Tunnels extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server) is installed.
2. Sign in with the same GitHub account used in the notebook.
3. Open Command Palette (`Ctrl+Shift+P` or `Cmd+Shift+P`).
4. Run `Remote Tunnels: Connect to Tunnel...` and select your notebook's tunnel.

You're now seamlessly connected to Colab through VS Code and can develop Nextflow pipelines more flexibly from just about any computer.

## Final tips for a Smooth Experience

Google Colab's storage is temporary and limited to around 100GB in most cases.
It's important to regularly back up your results to avoid data loss—mounting your personal Google Drive is convenient for small to moderate outputs, but may not be suitable for large workflow results, which can reach hundreds of gigabytes.
For larger datasets, consider syncing to external cloud storage or transferring results to institutional or project-specific storage solutions.
Additionally, if you plan on writing and developing your pipelines exclusively in Google Colab, make sure to use `git` and regularly commit and push your code, or alternatively, test in Colab but save changes from your local PC and commit to prevent loss of your progress after the notebook instance shuts down.

For long-term storage of results and access to your own datasets, you can mount your Google Drive in Colab:

```python
from google.colab import drive
drive.mount('/content/drive')
```

### Limitations to Note

And remember, while Google Colab is a powerful and accessible platform, it does have some constraints:

- Session timeouts and limited runtime duration
- No root access for system-level changes
- Limited resources compared to dedicated cloud VMs

## Conclusion

In this guide, we covered:

- Setting up Nextflow and dependencies in Google Colab
- Running nf-core pipelines in a cloud environment
- Connecting Colab to VS Code for a richer development experience
