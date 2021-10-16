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
const fs = require("fs");
const path = require("path");
const tcpport = require("tcp-port-used");
const pf = require("portfinder");
const logger = require("../logger");
const utils = require("../utils");
const track = require("../track");
const registry_1 = require("../emulator/registry");
const types_1 = require("../emulator/types");
const constants_1 = require("../emulator/constants");
const functionsEmulator_1 = require("../emulator/functionsEmulator");
const databaseEmulator_1 = require("../emulator/databaseEmulator");
const firestoreEmulator_1 = require("../emulator/firestoreEmulator");
const hostingEmulator_1 = require("../emulator/hostingEmulator");
const error_1 = require("../error");
const getProjectId = require("../getProjectId");
const pubsubEmulator_1 = require("./pubsubEmulator");
const commandUtils = require("./commandUtils");
const hub_1 = require("./hub");
const hubExport_1 = require("./hubExport");
const gui_1 = require("./gui");
const previews = require("../previews");
function checkPortOpen(port, host) {
    return __awaiter(this, void 0, void 0, function* () {
        try {
            const inUse = yield tcpport.check(port, host);
            return !inUse;
        }
        catch (e) {
            logger.debug(`port check error: ${e}`);
            return false;
        }
    });
}
exports.checkPortOpen = checkPortOpen;
function waitForPortClosed(port, host) {
    return __awaiter(this, void 0, void 0, function* () {
        const interval = 250;
        const timeout = 30000;
        try {
            yield tcpport.waitUntilUsedOnHost(port, host, interval, timeout);
        }
        catch (e) {
            throw new error_1.FirebaseError(`TIMEOUT: Port ${port} on ${host} was not active within ${timeout}ms`);
        }
    });
}
exports.waitForPortClosed = waitForPortClosed;
function startEmulator(instance) {
    return __awaiter(this, void 0, void 0, function* () {
        const name = instance.getName();
        const { host, port } = instance.getInfo();
        track("emulators:start", name);
        const portOpen = yield checkPortOpen(port, host);
        if (!portOpen) {
            yield cleanShutdown();
            const description = constants_1.Constants.description(name);
            utils.logLabeledWarning(name, `Port ${port} is not open on ${host}, could not start ${description}.`);
            utils.logLabeledBullet(name, `To select a different host/port for the emulator, specify that host/port in a firebase.json config file:
    {
      // ...
      "emulators": {
        "${name}": {
          "host": "${clc.yellow("HOST")}",
          "port": "${clc.yellow("PORT")}"
        }
      }
    }`);
            return utils.reject(`Could not start ${name} emulator, port taken.`, {});
        }
        yield registry_1.EmulatorRegistry.start(instance);
    });
}
exports.startEmulator = startEmulator;
function cleanShutdown() {
    return __awaiter(this, void 0, void 0, function* () {
        utils.logLabeledBullet("emulators", "Shutting down emulators.");
        for (const name of registry_1.EmulatorRegistry.listRunning()) {
            utils.logLabeledBullet(name, `Stopping ${constants_1.Constants.description(name)}`);
            yield registry_1.EmulatorRegistry.stop(name);
        }
        return true;
    });
}
exports.cleanShutdown = cleanShutdown;
function filterEmulatorTargets(options) {
    let targets = types_1.ALL_SERVICE_EMULATORS.filter((e) => {
        return options.config.has(e) || options.config.has(`emulators.${e}`);
    });
    if (options.only) {
        targets = _.intersection(targets, options.only.split(","));
    }
    return targets;
}
exports.filterEmulatorTargets = filterEmulatorTargets;
function shouldStart(options, name) {
    if (name === types_1.Emulators.HUB) {
        return !!options.project;
    }
    const targets = filterEmulatorTargets(options);
    if (name === types_1.Emulators.GUI) {
        return (previews.emulatorgui &&
            !!options.project &&
            targets.some((target) => types_1.EMULATORS_SUPPORTED_BY_GUI.indexOf(target) >= 0));
    }
    return targets.indexOf(name) >= 0;
}
exports.shouldStart = shouldStart;
function startAll(options, noGui = false) {
    return __awaiter(this, void 0, void 0, function* () {
        const targets = filterEmulatorTargets(options);
        options.targets = targets;
        const projectId = getProjectId(options, true);
        utils.logLabeledBullet("emulators", `Starting emulators: ${targets.join(", ")}`);
        if (options.only) {
            const requested = options.only.split(",");
            const ignored = _.difference(requested, targets);
            for (const name of ignored) {
                utils.logLabeledWarning(name, `Not starting the ${clc.bold(name)} emulator, make sure you have run ${clc.bold("firebase init")}.`);
            }
        }
        if (shouldStart(options, types_1.Emulators.HUB)) {
            const hubAddr = constants_1.Constants.getAddress(types_1.Emulators.HUB, options);
            const hubPort = yield pf.getPortPromise({
                host: hubAddr.host,
                port: hubAddr.port,
                stopPort: hubAddr.port + 100,
            });
            if (hubPort != hubAddr.port) {
                utils.logLabeledWarning("emulators", `${constants_1.Constants.description(types_1.Emulators.HUB)} unable to start on port ${hubAddr.port}, starting on ${hubPort}`);
            }
            const hub = new hub_1.EmulatorHub({
                projectId,
                host: hubAddr.host,
                port: hubPort,
            });
            yield startEmulator(hub);
        }
        let exportMetadata = {
            version: "unknown",
        };
        if (options.import) {
            const importDir = path.resolve(options.import);
            exportMetadata = JSON.parse(fs.readFileSync(path.join(importDir, hubExport_1.HubExport.METADATA_FILE_NAME), "utf8").toString());
        }
        if (shouldStart(options, types_1.Emulators.FUNCTIONS)) {
            const functionsAddr = constants_1.Constants.getAddress(types_1.Emulators.FUNCTIONS, options);
            const projectId = getProjectId(options, false);
            const functionsDir = path.join(options.extensionDir || options.config.projectDir, options.config.get("functions.source"));
            let inspectFunctions;
            if (options.inspectFunctions) {
                inspectFunctions = commandUtils.parseInspectionPort(options);
                utils.logLabeledWarning("functions", `You are running the functions emulator in debug mode (port=${inspectFunctions}). This means that functions will execute in sequence rather than in parallel.`);
            }
            const functionsEmulator = new functionsEmulator_1.FunctionsEmulator({
                projectId,
                functionsDir,
                host: functionsAddr.host,
                port: functionsAddr.port,
                debugPort: inspectFunctions,
                env: options.extensionEnv,
                predefinedTriggers: options.extensionTriggers,
            });
            yield startEmulator(functionsEmulator);
        }
        if (shouldStart(options, types_1.Emulators.FIRESTORE)) {
            const firestoreAddr = constants_1.Constants.getAddress(types_1.Emulators.FIRESTORE, options);
            const args = {
                host: firestoreAddr.host,
                port: firestoreAddr.port,
                projectId,
                auto_download: true,
            };
            if (exportMetadata.firestore) {
                const importDirAbsPath = path.resolve(options.import);
                const exportMetadataFilePath = path.join(importDirAbsPath, exportMetadata.firestore.metadata_file);
                utils.logLabeledBullet("firestore", `Importing data from ${exportMetadataFilePath}`);
                args.seed_from_export = exportMetadataFilePath;
            }
            const rulesLocalPath = options.config.get("firestore.rules");
            const foundRulesFile = rulesLocalPath && fs.existsSync(rulesLocalPath);
            if (rulesLocalPath) {
                const rules = path.join(options.projectRoot, rulesLocalPath);
                if (fs.existsSync(rules)) {
                    args.rules = rules;
                }
                else {
                    utils.logLabeledWarning("firestore", `Cloud Firestore rules file ${clc.bold(rules)} specified in firebase.json does not exist.`);
                }
            }
            else {
                utils.logLabeledWarning("firestore", "Did not find a Cloud Firestore rules file specified in a firebase.json config file.");
            }
            if (!foundRulesFile) {
                utils.logLabeledWarning("firestore", "The emulator will default to allowing all reads and writes. Learn more about this option: https://firebase.google.com/docs/emulator-suite/install_and_configure#security_rules_configuration.");
            }
            const firestoreEmulator = new firestoreEmulator_1.FirestoreEmulator(args);
            yield startEmulator(firestoreEmulator);
            utils.logLabeledBullet(types_1.Emulators.FIRESTORE, `For testing set ${clc.bold(`${firestoreEmulator_1.FirestoreEmulator.FIRESTORE_EMULATOR_ENV}=${firestoreAddr.host}:${firestoreAddr.port}`)}`);
        }
        if (shouldStart(options, types_1.Emulators.DATABASE)) {
            const databaseAddr = constants_1.Constants.getAddress(types_1.Emulators.DATABASE, options);
            const args = {
                host: databaseAddr.host,
                port: databaseAddr.port,
                projectId,
                auto_download: true,
            };
            if (shouldStart(options, types_1.Emulators.FUNCTIONS)) {
                const functionsAddr = constants_1.Constants.getAddress(types_1.Emulators.FUNCTIONS, options);
                args.functions_emulator_host = functionsAddr.host;
                args.functions_emulator_port = functionsAddr.port;
            }
            const rulesLocalPath = options.config.get("database.rules");
            const foundRulesFile = rulesLocalPath && fs.existsSync(rulesLocalPath);
            if (rulesLocalPath) {
                const rules = path.join(options.projectRoot, rulesLocalPath);
                if (fs.existsSync(rules)) {
                    args.rules = rules;
                }
                else {
                    utils.logLabeledWarning("database", `Realtime Database rules file ${clc.bold(rules)} specified in firebase.json does not exist.`);
                }
            }
            else {
                utils.logLabeledWarning("database", "Did not find a Realtime Database rules file specified in a firebase.json config file.");
            }
            if (!foundRulesFile) {
                utils.logLabeledWarning("database", "The emulator will default to allowing all reads and writes. Learn more about this option: https://firebase.google.com/docs/emulator-suite/install_and_configure#security_rules_configuration.");
            }
            const databaseEmulator = new databaseEmulator_1.DatabaseEmulator(args);
            yield startEmulator(databaseEmulator);
            utils.logLabeledBullet(types_1.Emulators.DATABASE, `For testing set ${clc.bold(`${databaseEmulator_1.DatabaseEmulator.DATABASE_EMULATOR_ENV}=${databaseAddr.host}:${databaseAddr.port}`)}`);
        }
        if (shouldStart(options, types_1.Emulators.HOSTING)) {
            const hostingAddr = constants_1.Constants.getAddress(types_1.Emulators.HOSTING, options);
            const hostingEmulator = new hostingEmulator_1.HostingEmulator({
                host: hostingAddr.host,
                port: hostingAddr.port,
                options,
            });
            yield startEmulator(hostingEmulator);
        }
        if (shouldStart(options, types_1.Emulators.PUBSUB)) {
            if (!projectId) {
                throw new error_1.FirebaseError("Cannot start the Pub/Sub emulator without a project: run 'firebase init' or provide the --project flag");
            }
            const pubsubAddr = constants_1.Constants.getAddress(types_1.Emulators.PUBSUB, options);
            const pubsubEmulator = new pubsubEmulator_1.PubsubEmulator({
                host: pubsubAddr.host,
                port: pubsubAddr.port,
                projectId,
                auto_download: true,
            });
            yield startEmulator(pubsubEmulator);
        }
        if (!noGui && shouldStart(options, types_1.Emulators.GUI)) {
            const guiAddr = constants_1.Constants.getAddress(types_1.Emulators.GUI, options);
            const guiPort = yield pf.getPortPromise({
                host: guiAddr.host,
                port: guiAddr.port,
                stopPort: guiAddr.port + 100,
            });
            if (guiPort != guiAddr.port) {
                utils.logLabeledWarning(types_1.Emulators.GUI, `${constants_1.Constants.description(types_1.Emulators.GUI)} unable to start on port ${guiAddr.port}, starting on ${guiPort}`);
            }
            const gui = new gui_1.EmulatorGUI({
                projectId,
                host: guiAddr.host,
                port: guiPort,
                auto_download: true,
            });
            yield startEmulator(gui);
        }
        const running = registry_1.EmulatorRegistry.listRunning();
        for (const name of running) {
            const instance = registry_1.EmulatorRegistry.get(name);
            if (instance) {
                yield instance.connect();
            }
        }
    });
}
exports.startAll = startAll;
