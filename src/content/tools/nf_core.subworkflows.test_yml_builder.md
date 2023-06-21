<!-- markdownlint-disable -->

<a href="../../nf_core/subworkflows/test_yml_builder.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nf_core.subworkflows.test_yml_builder`
The ModulesTestYmlBuilder class handles automatic generation of the modules test.yml file along with running the tests and creating md5 sums 



---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L36"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `SubworkflowTestYmlBuilder`




<a href="../../nf_core/subworkflows/test_yml_builder.py#L37"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(
    subworkflow=None,
    directory='.',
    run_tests=False,
    test_yml_output_path=None,
    force_overwrite=False,
    no_prompts=False,
    remote_url=None,
    branch=None
)
```








---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L144"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `build_all_tests`

```python
build_all_tests()
```

Go over each entry point and build structure 

---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L153"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `build_single_test`

```python
build_single_test(entry_point)
```

Given the supplied cli flags, prompt for any that are missing. 

Returns: Test command 

---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L224"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_if_empty_file`

```python
check_if_empty_file(fname)
```

Check if the file is empty, or compressed empty 

---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L79"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `check_inputs`

```python
check_inputs()
```

Do more complex checks about supplied flags. 

---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L247"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `create_test_file_dict`

```python
create_test_file_dict(results_dir, is_repeat=False)
```

Walk through directory and collect md5 sums 

---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L275"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `get_md5_sums`

```python
get_md5_sums(command, results_dir=None, results_dir_repeat=None)
```

Recursively go through directories and subdirectories and generate tuples of (<file_path>, <md5sum>) returns: list of tuples 

---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L202"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `parse_module_tags`

```python
parse_module_tags()
```

Parse the subworkflow test main.nf file to retrieve all imported modules for adding tags. 

---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L372"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `print_test_yml`

```python
print_test_yml()
```

Generate the test yml file. 

---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L65"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `run`

```python
run()
```

Run build steps 

---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L316"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `run_tests_workflow`

```python
run_tests_workflow(command)
```

Given a test workflow and an entry point, run the test workflow 

---

<a href="../../nf_core/subworkflows/test_yml_builder.py#L133"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `scrape_workflow_entry_points`

```python
scrape_workflow_entry_points()
```

Find the test workflow entry points from main.nf 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
