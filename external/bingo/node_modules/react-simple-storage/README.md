# React Simple Storage

A simple component and helper functions for using web storage with React.

[Check out the demo app and basic example. https://ryanjyost.github.io/react-simple-storage-example-project](https://ryanjyost.github.io/react-simple-storage-example-project/)

You may also want to check out my related Hacker Noon article, [How to take advantage of Local Storage in your React projects](https://hackernoon.com/how-to-take-advantage-of-local-storage-in-your-react-projects-a895f2b2d3f2), for some context and logic behind this project.

#### Good use cases for react-simple-storage
* Persist and experiment with a component's state while developing.
* Save form data across user sessions.
* A simple, quick fake backend for a practice or portfolio project.
* More I can't think of... 

## Install

[Install via yarn.](https://www.npmjs.com/package/yarn)
```
yarn add react-simple-storage
```
#### Using on IE11
For react-simple-storage to work on IE11, you'll need to use [babel-polyfill](https://babeljs.io/docs/usage/polyfill/).
```
yarn add babel-polyfill
```
Then import in your project.
```
import "babel-polyfill";
```

## Usage

### Component

Import and include an instance of react-simple-storage in a component whose state you want to save to web storage.
```javascript
import React, { Component } from "react";
import SimpleStorage from "react-simple-storage";

export default class ParentComponent extends Component {
  constructor(props) { 
    super(props)
    this.state = {
      text: "",
    }
  }

  render() {
    return ( 
      <div>
      
        // include the component somewhere in the parent to save the parent's state in web storage
        <SimpleStorage parent={this} />

        // the value of this input will be saved in web storage
        <input
          type="text"
          value={this.state.text}
          onChange={e => this.setState({ text: e.target.value })}
        />
        
      </div>
    ) 
  }
}
```

### Props
| Name             | Type            |Required? | Default      | Description
| ---------------- |:--------------- |:-------- | ------------ |-------------
| parent           | *object*        | Yes      | **none**     | reference to the parent component, i.e. `this`
| prefix           | *string*        | No       | ""           | prefix added to storage keys to avoid name clashes across instances     
| blacklist        | *array*         | No       | []           | a list of parent component's `state` names/keys to ignore when saving to storage
| onParentStateHydrated        | *func*         | No       | none           | fires after the parent component's `state` has been updated with storage items. Basically a callback for working with the parent component's `state` once updated with storage.




## Helper Functions
### `clearStorage(prefix)`
Clears items in `storage` with the given `prefix`, or all items if no `prefix` is given.
* `prefix: String | optional` - Corresponds to `prefix` prop passed to an instance of the `react-simple-storage` 
component.

#### Example
```javascript
import React, { Component } from "react";
import SimpleStorage, { clearStorage } from "react-simple-storage";

export default class ParentComponent extends Component {
  constructor(props) { 
    super(props)
    this.state = {
      text: "",
    }
  }

  render() {
    return ( 
      <div>
      
        // provide a prefix prop to be able to clear just the storage items 
        // created by this instance of the react-simple-storage component
        <SimpleStorage parent={this} prefix={"ParentComponent"} />

        <input
          type="text"
          value={this.state.text}
          onChange={e => this.setState({ text: e.target.value })}
        />
        
        // removes only storage items related to the ParentComponent
        <button onClick={() => clearStorage("ParentComponent")}>
          Clear storage for ParentComponent
        </button>
        
         // removes all items from storage
        <button onClick={() => clearStorage()}>
          Clear all storage
        </button>
        
      </div>
    ) 
  }
}
```


### `resetParentState(parent, initialState, keysToIgnore)`
Resets the parent's state to given `initialState`. 
* `parent: Object | required` - Reference to the parent component, allowing `react-simple-storage` to access and update 
the parent component's state. If called within the parent component, simply pass `this`.
* `initialState: Object | required` - The `state` of the parent component after the function executes.
* `keysToIgnore: Array | optional` - A list of keys in the parent component's `state` to ignore on `resetParentState
`. These pieces of that parent's state will NOT be reset.

#### Example

```javascript
import React, { Component } from "react";
import SimpleStorage, { resetParentState } from "react-simple-storage";

export default class ParentComponent extends Component {
  constructor(props) { 
    super(props)
    this.state = {
      text: "Initial Text",
    }
    
    // store the component's initial state to reset it
    this.initialState = this.state;
  }

  render() {
    return ( 
      <div>
      
        <SimpleStorage parent={this} />

        <input
          type="text"
          value={this.state.text}
          onChange={e => this.setState({ text: e.target.value })}
        />
        
        // will set "text" in state to "Initial Text"
         <button onClick={() => resetParentState(this, this.initialState)}>
           Reset parent state
         </button>
        
        // ignores "text" on reset, so will have no effect here
        <button onClick={() => resetParentState(this, this.initialState, ['text'])}>
          Do NOT reset text
        </button>
        
      </div>
    ) 
  }
}
```

## Built with
* [store.js](https://github.com/marcuswestin/store.js) - Cross-browser storage for all use cases, used across the web.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
