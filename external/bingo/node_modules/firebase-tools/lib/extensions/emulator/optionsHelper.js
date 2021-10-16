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
const fs = require("fs-extra");
const _ = require("lodash");
const path = require("path");
const paramHelper = require("./paramHelper");
const specHelper = require("./specHelper");
const triggerHelper = require("./triggerHelper");
const extensionsHelper = require("../extensionsHelper");
const Config = require("../../config");
const error_1 = require("../../error");
const emulatorLogger_1 = require("../../emulator/emulatorLogger");
const getProjectId = require("../../getProjectId");
function buildOptions(options) {
    return __awaiter(this, void 0, void 0, function* () {
        const extensionDir = yield specHelper.findExtensionYaml(process.cwd());
        options.extensionDir = extensionDir;
        const extensionYaml = yield specHelper.readExtensionYaml(extensionDir);
        extensionsHelper.validateSpec(extensionYaml);
        const params = yield paramHelper.readParamsFile(options.testParams);
        extensionsHelper.validateCommandLineParams(params, extensionYaml.params || []);
        params["PROJECT_ID"] = getProjectId(options, false);
        params["EXT_INSTANCE_ID"] = params["EXT_INSTANCE_ID"] || extensionYaml.name;
        params["DATABASE_INSTANCE"] = params["PROJECT_ID"];
        params["DATABASE_URL"] = `https://${params["DATABASE_INSTANCE"]}.firebaseio.com`;
        params["STORAGE_BUCKET"] = `${params["PROJECT_ID"]}.appspot.com`;
        const functionResources = specHelper.getFunctionResourcesWithParamSubstitution(extensionYaml, params);
        let testConfig;
        if (options.testConfig) {
            testConfig = readTestConfigFile(options.testConfig);
            checkTestConfig(testConfig, functionResources);
        }
        options.config = buildConfig(functionResources, testConfig);
        options.extensionEnv = params;
        const functionEmuTriggerDefs = functionResources.map((r) => triggerHelper.functionResourceToEmulatedTriggerDefintion(r));
        options.extensionTriggers = functionEmuTriggerDefs;
        return options;
    });
}
exports.buildOptions = buildOptions;
function checkTestConfig(testConfig, functionResources) {
    if (!testConfig.functions && functionResources.length) {
        emulatorLogger_1.EmulatorLogger.log("WARN", "This extension uses functions," +
            "but 'firebase.json' provided by --test-config is missing a top-level 'functions' object." +
            "Functions will not be emulated.");
    }
    if (!testConfig.firestore && shouldEmulateFirestore(functionResources)) {
        emulatorLogger_1.EmulatorLogger.log("WARN", "This extension interacts with Cloud Firestore," +
            "but 'firebase.json' provided by --test-config is missing a top-level 'firestore' object." +
            "Cloud Firestore will not be emulated.");
    }
    if (!testConfig.database && shouldEmulateDatabase(functionResources)) {
        emulatorLogger_1.EmulatorLogger.log("WARN", "This extension interacts with Realtime Database," +
            "but 'firebase.json' provided by --test-config is missing a top-level 'database' object." +
            "Realtime Database will not be emulated.");
    }
}
function readTestConfigFile(testConfigPath) {
    try {
        const buf = fs.readFileSync(path.resolve(testConfigPath));
        return JSON.parse(buf.toString());
    }
    catch (err) {
        throw new error_1.FirebaseError(`Error reading --test-config file: ${err.message}\n`, {
            original: err,
        });
    }
}
function buildConfig(functionResources, testConfig) {
    const config = new Config(testConfig || {}, { projectDir: process.cwd(), cwd: process.cwd() });
    const emulateFunctions = shouldEmulateFunctions(functionResources);
    if (!testConfig) {
        if (emulateFunctions) {
            config.set("functions", {});
        }
        if (shouldEmulateFirestore(functionResources)) {
            config.set("firestore", {});
        }
        if (shouldEmulateDatabase(functionResources)) {
            config.set("database", {});
        }
        if (shouldEmulatePubsub(functionResources)) {
            config.set("pubsub", {});
        }
    }
    if (config.get("functions")) {
        const sourceDirectory = getFunctionSourceDirectory(functionResources);
        config.set("functions.source", sourceDirectory);
    }
    return config;
}
function getFunctionSourceDirectory(functionResources) {
    let sourceDirectory;
    for (let r of functionResources) {
        let dir = _.get(r, "properties.sourceDirectory");
        if (!dir) {
            emulatorLogger_1.EmulatorLogger.log("INFO", `No sourceDirectory was specified for function ${r.name}, defaulting to 'functions'`);
            dir = "functions";
        }
        if (!sourceDirectory) {
            sourceDirectory = dir;
        }
        else if (sourceDirectory != dir) {
            throw new error_1.FirebaseError(`Found function resources with different sourceDirectories: '${sourceDirectory}' and '${dir}'. The extensions emulator only supports a single sourceDirectory.`);
        }
    }
    return sourceDirectory;
}
function shouldEmulateFunctions(resources) {
    return resources.length > 0;
}
function shouldEmulate(emulatorName, resources) {
    for (const r of resources) {
        const eventType = _.get(r, "properties.eventTrigger.eventType", "");
        if (eventType.includes(emulatorName)) {
            return true;
        }
    }
    return false;
}
function shouldEmulateFirestore(resources) {
    return shouldEmulate("cloud.firestore", resources);
}
function shouldEmulateDatabase(resources) {
    return shouldEmulate("google.firebase.database", resources);
}
function shouldEmulatePubsub(resources) {
    return shouldEmulate("google.pubsub", resources);
}
