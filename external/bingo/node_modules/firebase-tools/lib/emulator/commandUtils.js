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
const clc = require("cli-color");
const childProcess = require("child_process");
const controller = require("../emulator/controller");
const Config = require("../config");
const utils = require("../utils");
const logger = require("../logger");
const requireAuth = require("../requireAuth");
const requireConfig = require("../requireConfig");
const types_1 = require("../emulator/types");
const error_1 = require("../error");
const databaseEmulator_1 = require("../emulator/databaseEmulator");
const registry_1 = require("../emulator/registry");
const firestoreEmulator_1 = require("../emulator/firestoreEmulator");
const getProjectId = require("../getProjectId");
const hub_1 = require("./hub");
exports.FLAG_ONLY = "--only <emulators>";
exports.DESC_ONLY = "only specific emulators. " +
    "This is a comma separated list of emulator names. " +
    "Valid options are: " +
    JSON.stringify(types_1.ALL_SERVICE_EMULATORS);
exports.FLAG_INSPECT_FUNCTIONS = "--inspect-functions [port]";
exports.DESC_INSPECT_FUNCTIONS = "emulate Cloud Functions in debug mode with the node inspector on the given port (9229 if not specified)";
exports.FLAG_IMPORT = "--import [dir]";
exports.DESC_IMPORT = "import emulator data from a previous export (see emulators:export)";
exports.FLAG_TEST_CONFIG = "--test-config <firebase.json file>";
exports.DESC_TEST_CONFIG = "A firebase.json style file. Used to configure the Firestore and Realtime Database emulators.";
exports.FLAG_TEST_PARAMS = "--test-params <params.env file>";
exports.DESC_TEST_PARAMS = "A .env file containing test param values for your emulated extension.";
const DEFAULT_CONFIG = new Config({ database: {}, firestore: {}, functions: {}, hosting: {} }, {});
function beforeEmulatorCommand(options) {
    return __awaiter(this, void 0, void 0, function* () {
        const optionsWithDefaultConfig = Object.assign(Object.assign({}, options), { config: DEFAULT_CONFIG });
        const optionsWithConfig = options.config ? options : optionsWithDefaultConfig;
        const canStartWithoutConfig = options.only &&
            !controller.shouldStart(optionsWithConfig, types_1.Emulators.FUNCTIONS) &&
            !controller.shouldStart(optionsWithConfig, types_1.Emulators.HOSTING);
        try {
            yield requireAuth(options);
        }
        catch (e) {
            logger.debug(e);
            utils.logLabeledWarning("emulators", `You are not currently authenticated so some features may not work correctly. Please run ${clc.bold("firebase login")} to authenticate the CLI.`);
        }
        if (canStartWithoutConfig && !options.config) {
            utils.logWarning("Could not find config (firebase.json) so using defaults.");
            options.config = DEFAULT_CONFIG;
        }
        else {
            yield requireConfig(options);
        }
    });
}
exports.beforeEmulatorCommand = beforeEmulatorCommand;
function parseInspectionPort(options) {
    let port = options.inspectFunctions;
    if (port === true) {
        port = "9229";
    }
    const parsed = Number(port);
    if (isNaN(parsed) || parsed < 1024 || parsed > 65535) {
        throw new error_1.FirebaseError(`"${port}" is not a valid port for debugging, please pass an integer between 1024 and 65535.`);
    }
    return parsed;
}
exports.parseInspectionPort = parseInspectionPort;
function runScript(script, extraEnv) {
    return __awaiter(this, void 0, void 0, function* () {
        utils.logBullet(`Running script: ${clc.bold(script)}`);
        const env = Object.assign(Object.assign({}, process.env), extraEnv);
        const databaseInstance = registry_1.EmulatorRegistry.get(types_1.Emulators.DATABASE);
        if (databaseInstance) {
            const info = databaseInstance.getInfo();
            const address = `${info.host}:${info.port}`;
            env[databaseEmulator_1.DatabaseEmulator.DATABASE_EMULATOR_ENV] = address;
        }
        const firestoreInstance = registry_1.EmulatorRegistry.get(types_1.Emulators.FIRESTORE);
        if (firestoreInstance) {
            const info = firestoreInstance.getInfo();
            const address = `${info.host}:${info.port}`;
            env[firestoreEmulator_1.FirestoreEmulator.FIRESTORE_EMULATOR_ENV] = address;
            env[firestoreEmulator_1.FirestoreEmulator.FIRESTORE_EMULATOR_ENV_ALT] = address;
        }
        const hubInstance = registry_1.EmulatorRegistry.get(types_1.Emulators.HUB);
        if (hubInstance) {
            const info = hubInstance.getInfo();
            const address = `${info.host}:${info.port}`;
            env[hub_1.EmulatorHub.EMULATOR_HUB_ENV] = address;
        }
        const proc = childProcess.spawn(script, {
            stdio: ["inherit", "inherit", "inherit"],
            shell: true,
            windowsHide: true,
            env,
        });
        logger.debug(`Running ${script} with environment ${JSON.stringify(env)}`);
        return new Promise((resolve, reject) => {
            proc.on("error", (err) => {
                utils.logWarning(`There was an error running the script: ${JSON.stringify(err)}`);
                reject();
            });
            const exitDelayMs = 500;
            proc.once("exit", (code, signal) => {
                if (signal) {
                    utils.logWarning(`Script exited with signal: ${signal}`);
                    setTimeout(reject, exitDelayMs);
                    return;
                }
                const exitCode = code || 0;
                if (code === 0) {
                    utils.logSuccess(`Script exited successfully (code 0)`);
                }
                else {
                    utils.logWarning(`Script exited unsuccessfully (code ${code})`);
                }
                setTimeout(() => {
                    resolve(exitCode);
                }, exitDelayMs);
            });
        });
    });
}
function emulatorExec(script, options) {
    return __awaiter(this, void 0, void 0, function* () {
        const projectId = getProjectId(options, true);
        const extraEnv = {};
        if (projectId) {
            extraEnv.GCLOUD_PROJECT = projectId;
        }
        let exitCode = 0;
        try {
            yield controller.startAll(options, true);
            exitCode = yield runScript(script, extraEnv);
        }
        finally {
            yield controller.cleanShutdown();
        }
        if (exitCode !== 0) {
            throw new error_1.FirebaseError(`Script "${clc.bold(script)}" exited with code ${exitCode}`, {
                exit: exitCode,
            });
        }
    });
}
exports.emulatorExec = emulatorExec;
