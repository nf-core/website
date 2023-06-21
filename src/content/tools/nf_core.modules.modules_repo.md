<!-- markdownlint-disable -->

<a href="../../nf_core/modules/modules_repo.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.modules.modules_repo`




**Global Variables**
---------------
- **NFCORE_CACHE_DIR**
- **NFCORE_DIR**
- **NF_CORE_MODULES_NAME**
- **NF_CORE_MODULES_REMOTE**
- **NF_CORE_MODULES_DEFAULT_BRANCH**


---

<a href="../../nf_core/modules/modules_repo.py#L25"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `ModulesRepo`
An object to store details about the repository being used for modules. 

Used by the `nf-core modules` top-level command with -r and -b flags, so that this can be used in the same way by all sub-commands. 

We keep track of the pull-status of the different installed repos in the static variable local_repo_status. This is so we don't need to pull a remote several times in one command. 

<a href="../../nf_core/modules/modules_repo.py#L40"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(remote_url=None, branch=None, no_pull=False, hide_progress=False)
```

Initializes the object and clones the git repository if it is not already present 




---

<a href="../../nf_core/modules/modules_repo.py#L74"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `setup_local_repo`

```python
setup_local_repo(remote, branch, hide_progress=True, in_cache=False)
```

Sets up the local git repository. If the repository has been cloned previously, it returns a git.Repo object of that clone. Otherwise it tries to clone the repository from the provided remote URL and returns a git.Repo of the new clone. 



**Args:**
 
 - <b>`remote`</b> (str):  git url of remote 
 - <b>`branch`</b> (str):  name of branch to use Sets self.repo 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
