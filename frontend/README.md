# Frontend Website for Data Analysis

This package contains the frontend website for analyzing simulations.

## Requirements
You need to install the node package manager (npm). It is part of Node.js. 

You can download it from here: https://nodejs.org/en/

We recommend using the long time support (LTS) version.

## Development and Release Builds
#### Getting the dependencies
To download the dependencies you have to run `npm install` or `npm i`.
This can take a few minutes.

#### Development
Run `npm run start` in this folder. 
Your default browser will open with the URL http://localhost:3000/. 
Every time you save a project file the project rebuilds, and the website automatically refreshes. 
If your default browser is Internet Explorer you have to open the URL (http://localhost:3000/) in another browser (Firefox, Edge or Chrome).

**Useful browser plugins for development:**
- React Developer Tools: [Firefox](https://addons.mozilla.org/de/firefox/addon/react-devtools/), [Chrome](https://chrome.google.com/webstore/detail/react-developer-tools/fmkadmapgofadopljbjfkapdkoienihi)
- Redux Developer Tools: [Firefox](https://addons.mozilla.org/de/firefox/addon/reduxdevtools/), [Chrome](https://chrome.google.com/webstore/detail/redux-devtools/lmhkpmbekcpmknklioeibfkpmmfibljd?hl=de)

#### Release Builds
Run `npm run build` to create a release build. The resulting web page can be found in the build folder.