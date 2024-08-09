# nf_core.create

A Textual app to create a pipeline.

### _`class{:python}`_`nf_core.pipelines.create.PipelineCreateApp(driver_class: Type[Driver] | None = None, css_path: str | PurePath | List[str | PurePath] | None = None, watch_css: bool = False){:python}`

Bases: `App`\[`CreateConfig`]

A Textual app to manage stopwatches.

#### `BINDINGS{:python}`_: ClassVar\[list\[BindingType]]_ _= \[('d', 'toggle_dark', 'Toggle dark mode'), ('q', 'quit', 'Quit')]_

#### `CSS_PATH{:python}`_: ClassVar\[CSSPathType | None]_ _= 'create.tcss'_

File paths to load CSS from.

#### `LOGGING_STATE{:python}`_= None_

#### `LOG_HANDLER{:python}`_= \<CustomLogHandler (INFO)>_

#### `NFCORE_PIPELINE{:python}`_= True_

#### `SCREENS{:python}`_: ClassVar\[dict\[str, Screen\[Any] | Callable\[\[], Screen\[Any]]]]_ _= {'basic_details': BasicDetails(), 'choose_type': ChoosePipelineType(), 'final_details': FinalDetails(), 'github_exit': GithubExit(), 'github_repo': GithubRepo(), 'github_repo_question': GithubRepoQuestion(), 'logging': LoggingScreen(), 'type_custom': CustomPipeline(), 'type_nfcore': NfcorePipeline(), 'welcome': WelcomeScreen()}_

Screens associated with the app for the lifetime of the app.

#### `SUB_TITLE{:python}`_: str | None_ _= 'Create a new pipeline with the nf-core pipeline template'_

A class variable to set the default sub-title for the application.

To update the sub-title while the app is running, you can set the \[sub_title]\[textual.app.App.sub_title] attribute.
See also \[the Screen.SUB_TITLE attribute]\[textual.screen.Screen.SUB_TITLE].

#### `TEMPLATE_CONFIG{:python}`_= CreateConfig(org=None, name=None, description=None, author=None, version=None, force=True, outdir=None, skip_features=None, is_nfcore=None)_

#### `TITLE{:python}`_: str | None_ _= 'nf-core create'_

A class variable to set the _default_ title for the application.

To update the title while the app is running, you can set the \[title]\[textual.app.App.title] attribute.
See also \[the Screen.TITLE attribute]\[textual.screen.Screen.TITLE].

#### `_computes{:python}`_: ClassVar\[frozenset\[str]]_ _= frozenset({})_

#### `_css_type_name{:python}`_: str_ _= 'PipelineCreateApp'_

#### `_css_type_names{:python}`_: ClassVar\[frozenset\[str]]_ _= frozenset({'App', 'DOMNode', 'PipelineCreateApp'})_

#### `_decorated_handlers{:python}`_: dict\[type\[Message], list\[tuple\[Callable, str | None]]]_ _= {}_

#### `_inherit_bindings{:python}`_: ClassVar\[bool]_ _= True_

#### `_inherit_component_classes{:python}`_: ClassVar\[bool]_ _= True_

#### `_inherit_css{:python}`_: ClassVar\[bool]_ _= True_

#### `_merged_bindings{:python}`_: ClassVar\[\_Bindings | None]_ _= \_Bindings({'ctrl+c': Binding(key='ctrl+c', action='quit', description='Quit', show=False, key_display=None, priority=True), 'ctrl+backslash': Binding(key='ctrl+backslash', action='command_palette', description='', show=False, key_display=None, priority=True), 'd': Binding(key='d', action='toggle_dark', description='Toggle dark mode', show=True, key_display=None, priority=False), 'q': Binding(key='q', action='quit', description='Quit', show=True, key_display=None, priority=False)})_

#### `_reactives{:python}`_: ClassVar\[dict\[str, Reactive]]_ _= {'ansi_theme_dark': Reactive(\<rich.terminal_theme.TerminalTheme object>, layout=False, repaint=True, init=False, always_update=False, compute=True, recompose=False), 'ansi_theme_light': Reactive(\<rich.terminal_theme.TerminalTheme object>, layout=False, repaint=True, init=False, always_update=False, compute=True, recompose=False), 'app_focus': Reactive(True, layout=False, repaint=True, init=False, always_update=False, compute=False, recompose=False), 'dark': Reactive(True, layout=False, repaint=True, init=False, always_update=False, compute=False, recompose=False), 'sub_title': Reactive('', layout=False, repaint=True, init=False, always_update=False, compute=False, recompose=False), 'title': Reactive('', layout=False, repaint=True, init=False, always_update=False, compute=False, recompose=False)}_

#### `action_toggle_dark(){:python}`

An action to toggle dark mode.

#### `on_button_pressed(event: Pressed){:python}`

Handle all button pressed events.

#### `on_mount(){:python}`
