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
const logger = require("../logger");
const parseTriggers = require("../parseTriggers");
const utils = require("../utils");
const os = require("os");
const path = require("path");
const fs = require("fs");
var EmulatedTriggerType;
(function (EmulatedTriggerType) {
    EmulatedTriggerType["BACKGROUND"] = "BACKGROUND";
    EmulatedTriggerType["HTTPS"] = "HTTPS";
})(EmulatedTriggerType = exports.EmulatedTriggerType || (exports.EmulatedTriggerType = {}));
const memoryLookup = {
    "128MB": 128,
    "256MB": 256,
    "512MB": 512,
    "1GB": 1024,
    "2GB": 2048,
};
class HttpConstants {
}
exports.HttpConstants = HttpConstants;
HttpConstants.CALLABLE_AUTH_HEADER = "x-callable-context-auth";
class EmulatedTrigger {
    constructor(definition, module) {
        this.definition = definition;
        this.module = module;
    }
    get memoryLimitBytes() {
        return memoryLookup[this.definition.availableMemoryMb || "128MB"] * 1024 * 1024;
    }
    get timeoutMs() {
        if (typeof this.definition.timeout === "number") {
            return this.definition.timeout * 1000;
        }
        else {
            return parseInt((this.definition.timeout || "60s").split("s")[0], 10) * 1000;
        }
    }
    getRawFunction() {
        if (!this.module) {
            throw new Error("EmulatedTrigger has not been provided a module.");
        }
        const func = _.get(this.module, this.definition.entryPoint);
        return func.__emulator_func || func;
    }
}
exports.EmulatedTrigger = EmulatedTrigger;
function getTriggersFromDirectory(projectId, functionsDir, firebaseConfig) {
    return __awaiter(this, void 0, void 0, function* () {
        let triggerDefinitions;
        try {
            triggerDefinitions = yield parseTriggers(projectId, functionsDir, {}, JSON.stringify(firebaseConfig));
        }
        catch (e) {
            utils.logWarning(`Failed to load functions source code.`);
            logger.info(e.message);
            return {};
        }
        return getEmulatedTriggersFromDefinitions(triggerDefinitions, functionsDir);
    });
}
exports.getTriggersFromDirectory = getTriggersFromDirectory;
function getEmulatedTriggersFromDefinitions(definitions, module) {
    return definitions.reduce((obj, definition) => {
        obj[definition.name] = new EmulatedTrigger(definition, module);
        return obj;
    }, {});
}
exports.getEmulatedTriggersFromDefinitions = getEmulatedTriggersFromDefinitions;
function getTemporarySocketPath(pid, cwd) {
    if (process.platform === "win32") {
        return path.join("\\\\?\\pipe", cwd, pid.toString());
    }
    else {
        return path.join(os.tmpdir(), `fire_emu_${pid.toString()}.sock`);
    }
}
exports.getTemporarySocketPath = getTemporarySocketPath;
function getFunctionRegion(def) {
    if (def.regions && def.regions.length > 0) {
        return def.regions[0];
    }
    return "us-central1";
}
exports.getFunctionRegion = getFunctionRegion;
function getFunctionService(def) {
    if (def.eventTrigger) {
        return def.eventTrigger.service;
    }
    return "unknown";
}
exports.getFunctionService = getFunctionService;
function waitForBody(req) {
    let data = "";
    return new Promise((resolve) => {
        req.on("data", (chunk) => {
            data += chunk;
        });
        req.on("end", () => {
            resolve(data);
        });
    });
}
exports.waitForBody = waitForBody;
function findModuleRoot(moduleName, filepath) {
    const hierarchy = filepath.split(path.sep);
    for (let i = 0; i < hierarchy.length; i++) {
        try {
            let chunks = [];
            if (i) {
                chunks = hierarchy.slice(0, -i);
            }
            else {
                chunks = hierarchy;
            }
            const packagePath = path.join(chunks.join(path.sep), "package.json");
            const serializedPackage = fs.readFileSync(packagePath, "utf8").toString();
            if (JSON.parse(serializedPackage).name === moduleName) {
                return chunks.join("/");
            }
            break;
        }
        catch (err) {
        }
    }
    return "";
}
exports.findModuleRoot = findModuleRoot;
