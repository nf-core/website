"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const _ = require("lodash");
const path = require("path");
const clc = require("cli-color");
const semver = require("semver");
const error_1 = require("./error");
const utils = require("./utils");
const cjson = require("cjson");
const MESSAGE_FRIENDLY_RUNTIMES = {
    nodejs6: "Node.js 6 (Deprecated)",
    nodejs8: "Node.js 8",
    nodejs10: "Node.js 10 (Beta)",
};
const ENGINE_RUNTIMES = {
    6: "nodejs6",
    8: "nodejs8",
    10: "nodejs10",
};
exports.ENGINES_FIELD_REQUIRED_MSG = clc.bold("Engines field is required in package.json but none was found.");
exports.UNSUPPORTED_NODE_VERSION_MSG = clc.bold(`package.json in functions directory has an engines field which is unsupported. ` +
    `The only valid choices are: ${clc.bold('{"node": "8"}')} and ${clc.bold('{"node": "10"}')}. ` +
    `Note that Node.js 6 is now deprecated.`);
exports.DEPRECATION_WARNING_MSG = clc.bold.yellow("functions: ") +
    "Deploying functions to Node 6 runtime, which is deprecated. Node 8 is available " +
    "and is the recommended runtime.";
exports.FUNCTIONS_SDK_VERSION_TOO_OLD_WARNING = clc.bold.yellow("functions: ") +
    "You must have a " +
    clc.bold("firebase-functions") +
    " version that is at least 2.0.0. Please run " +
    clc.bold("npm i --save firebase-functions@latest") +
    " in the functions folder.";
function getHumanFriendlyRuntimeName(runtime) {
    return _.get(MESSAGE_FRIENDLY_RUNTIMES, runtime, runtime);
}
exports.getHumanFriendlyRuntimeName = getHumanFriendlyRuntimeName;
function getRuntimeChoice(sourceDir) {
    const packageJsonPath = path.join(sourceDir, "package.json");
    const loaded = cjson.load(packageJsonPath);
    const engines = loaded.engines;
    if (!engines || !engines.node) {
        return null;
    }
    const runtime = ENGINE_RUNTIMES[engines.node];
    if (!runtime) {
        throw new error_1.FirebaseError(exports.UNSUPPORTED_NODE_VERSION_MSG, { exit: 1 });
    }
    if (runtime === "nodejs6") {
        utils.logWarning(exports.DEPRECATION_WARNING_MSG);
    }
    else {
        if (functionsSDKTooOld(loaded)) {
            utils.logWarning(exports.FUNCTIONS_SDK_VERSION_TOO_OLD_WARNING);
        }
    }
    return runtime;
}
exports.getRuntimeChoice = getRuntimeChoice;
function functionsSDKTooOld(loaded) {
    const SDKRange = _.get(loaded, "dependencies.firebase-functions");
    try {
        if (!semver.intersects(SDKRange, ">=2")) {
            return true;
        }
    }
    catch (e) {
    }
    return false;
}
