# Frontend Website for Data Analysis

This package contains the frontend website for analyzing simulations.

## Requirements
#### NPM
You need to install the node package manager (npm). It is part of Node.js. 

You can download it from here: https://nodejs.org/en/

We recommend using the long time support (LTS) version.

#### Emscripten
We use a WebAssembly module, which requires compilation from C++. 
This is done with the Emscripten SDK. 
Installation instructions can be found here:

https://emscripten.org/docs/getting_started/downloads.html

#### CMake
The WASM compilation also needs CMake. You can get the latest version here: 

https://cmake.org/download/

## Development and Release Builds
#### Getting the dependencies
**NPM**

To download the dependencies you have to run `npm install` or `npm i`.
This can take a few minutes.
[In case of any warnings to required peer dependencies of react, 
add the flag --save-dev react@^VERSION to npm install. Please also report
the bug to us.]

**Building the WebAssembly module**
1. Enable your Emscripten SDK by executing their `emsdk_dev` script.
2. Run `npm run build:wasm`

> Note: If you want to debug the WebAssembly module you can enable debug mode by executing
>
> `npm run build:wasm:debug`
>
> in step 2. You can then also check for memory leaks from JavaScript by executing the `doLeakCheck()` function of the WASM module.

#### Development
Run `npm run start` in this folder. 
Your default browser will open with the URL http://localhost:3000/. 
Every time you save a project file the project rebuilds, and the website automatically refreshes. 
If your default browser is Internet Explorer you have to open the URL (http://localhost:3000/) in another browser (Firefox, Edge or Chrome).

By default, the development version starts in `development` mode. 
To override this you can create a file `.env.development.local` and insert `REACT_APP_MODE=production`.
You have to restart the development server for this to have an effect.

**Useful browser plugins for development:**
- React Developer Tools: [Firefox](https://addons.mozilla.org/de/firefox/addon/react-devtools/), [Chrome](https://chrome.google.com/webstore/detail/react-developer-tools/fmkadmapgofadopljbjfkapdkoienihi)
- Redux Developer Tools: [Firefox](https://addons.mozilla.org/de/firefox/addon/reduxdevtools/), [Chrome](https://chrome.google.com/webstore/detail/redux-devtools/lmhkpmbekcpmknklioeibfkpmmfibljd?hl=de)

#### Release Builds
Run `npm run build` to create a release build. The resulting web page can be found in the build folder.