"use strict";
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
Object.defineProperty(exports, "__esModule", { value: true });
const _ = require("lodash");
const clc = require("cli-color");
const stream_1 = require("stream");
const configstore_1 = require("./configstore");
const error_1 = require("./error");
const logger = require("./logger");
const IS_WINDOWS = process.platform === "win32";
const SUCCESS_CHAR = IS_WINDOWS ? "+" : "✔";
const WARNING_CHAR = IS_WINDOWS ? "!" : "⚠";
exports.envOverrides = [];
function consoleUrl(project, path) {
    const api = require("./api");
    return `${api.consoleOrigin}/project/${project}${path}`;
}
exports.consoleUrl = consoleUrl;
function getInheritedOption(options, key) {
    let target = options;
    while (target) {
        if (_.has(target, key)) {
            return target[key];
        }
        target = target.parent;
    }
}
exports.getInheritedOption = getInheritedOption;
function envOverride(envname, value, coerce) {
    const currentEnvValue = process.env[envname];
    if (currentEnvValue && currentEnvValue.length) {
        exports.envOverrides.push(envname);
        if (coerce) {
            try {
                return coerce(currentEnvValue, value);
            }
            catch (e) {
                return value;
            }
        }
        return currentEnvValue;
    }
    return value;
}
exports.envOverride = envOverride;
function addSubdomain(origin, subdomain) {
    return origin.replace("//", `//${subdomain}.`);
}
exports.addSubdomain = addSubdomain;
function logSuccess(message, type = "info") {
    logger[type](clc.green.bold(`${SUCCESS_CHAR} `), message);
}
exports.logSuccess = logSuccess;
function logLabeledSuccess(label, message, type = "info") {
    logger[type](clc.green.bold(`${SUCCESS_CHAR}  ${label}:`), message);
}
exports.logLabeledSuccess = logLabeledSuccess;
function logBullet(message, type = "info") {
    logger[type](clc.cyan.bold("i "), message);
}
exports.logBullet = logBullet;
function logLabeledBullet(label, message, type = "info") {
    logger[type](clc.cyan.bold(`i  ${label}:`), message);
}
exports.logLabeledBullet = logLabeledBullet;
function logWarning(message, type = "warn") {
    logger[type](clc.yellow.bold(`${WARNING_CHAR} `), message);
}
exports.logWarning = logWarning;
function logLabeledWarning(label, message, type = "warn") {
    logger[type](clc.yellow.bold(`${WARNING_CHAR}  ${label}:`), message);
}
exports.logLabeledWarning = logLabeledWarning;
function reject(message, options) {
    return Promise.reject(new error_1.FirebaseError(message, options));
}
exports.reject = reject;
function explainStdin() {
    if (IS_WINDOWS) {
        throw new error_1.FirebaseError("STDIN input is not available on Windows.", {
            exit: 1,
        });
    }
    if (process.stdin.isTTY) {
        logger.info(clc.bold("Note:"), "Reading STDIN. Type JSON data and then press Ctrl-D");
    }
}
exports.explainStdin = explainStdin;
function stringToStream(text) {
    if (!text) {
        return undefined;
    }
    const s = new stream_1.Readable();
    s.push(text);
    s.push(null);
    return s;
}
exports.stringToStream = stringToStream;
function makeActiveProject(projectDir, newActive) {
    const activeProjects = configstore_1.configstore.get("activeProjects") || {};
    if (newActive) {
        activeProjects[projectDir] = newActive;
    }
    else {
        _.unset(activeProjects, projectDir);
    }
    configstore_1.configstore.set("activeProjects", activeProjects);
}
exports.makeActiveProject = makeActiveProject;
function endpoint(parts) {
    return `/${_.join(parts, "/")}`;
}
exports.endpoint = endpoint;
function getFunctionsEventProvider(eventType) {
    const parts = eventType.split("/");
    if (parts.length > 1) {
        const provider = _.last(parts[1].split("."));
        return _.capitalize(provider);
    }
    if (eventType.match(/google.pubsub/)) {
        return "PubSub";
    }
    else if (eventType.match(/google.storage/)) {
        return "Storage";
    }
    else if (eventType.match(/google.analytics/)) {
        return "Analytics";
    }
    else if (eventType.match(/google.firebase.database/)) {
        return "Database";
    }
    else if (eventType.match(/google.firebase.auth/)) {
        return "Auth";
    }
    else if (eventType.match(/google.firebase.crashlytics/)) {
        return "Crashlytics";
    }
    else if (eventType.match(/google.firestore/)) {
        return "Firestore";
    }
    return _.capitalize(eventType.split(".")[1]);
}
exports.getFunctionsEventProvider = getFunctionsEventProvider;
function promiseAllSettled(promises) {
    const wrappedPromises = _.map(promises, (p) => __awaiter(this, void 0, void 0, function* () {
        try {
            const val = yield Promise.resolve(p);
            return { state: "fulfilled", value: val };
        }
        catch (err) {
            return { state: "rejected", reason: err };
        }
    }));
    return Promise.all(wrappedPromises);
}
exports.promiseAllSettled = promiseAllSettled;
function promiseWhile(action, check, interval = 2500) {
    return __awaiter(this, void 0, void 0, function* () {
        return new Promise((resolve, promiseReject) => {
            const run = () => __awaiter(this, void 0, void 0, function* () {
                try {
                    const res = yield action();
                    if (check(res)) {
                        return resolve(res);
                    }
                    setTimeout(run, interval);
                }
                catch (err) {
                    return promiseReject(err);
                }
            });
            run();
        });
    });
}
exports.promiseWhile = promiseWhile;
function promiseProps(obj) {
    return __awaiter(this, void 0, void 0, function* () {
        const resultObj = {};
        const promises = _.keys(obj).map((key) => __awaiter(this, void 0, void 0, function* () {
            const r = yield Promise.resolve(obj[key]);
            resultObj[key] = r;
        }));
        return Promise.all(promises).then(() => resultObj);
    });
}
exports.promiseProps = promiseProps;
