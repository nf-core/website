"use strict";
module.exports = {
    database: require("./database"),
    firestore: require("./firestore").doSetup,
    functions: require("./functions"),
    hosting: require("./hosting"),
    storage: require("./storage").doSetup,
    emulators: require("./emulators").doSetup,
    project: require("./project").doSetup,
};
