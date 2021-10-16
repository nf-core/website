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
const types_1 = require("./types");
const constants_1 = require("./constants");
const error_1 = require("../error");
const childProcess = require("child_process");
const utils = require("../utils");
const logger = require("../logger");
const clc = require("cli-color");
const fs = require("fs-extra");
const path = require("path");
const os = require("os");
const downloadEmulator = require("../emulator/download");
const EMULATOR_INSTANCE_KILL_TIMEOUT = 2000;
const CACHE_DIR = process.env.FIREBASE_EMULATORS_PATH || path.join(os.homedir(), ".cache", "firebase", "emulators");
const DownloadDetails = {
    database: {
        downloadPath: path.join(CACHE_DIR, "firebase-database-emulator-v4.4.1.jar"),
        version: "4.4.1",
        opts: {
            cacheDir: CACHE_DIR,
            remoteUrl: "https://storage.googleapis.com/firebase-preview-drop/emulator/firebase-database-emulator-v4.4.1.jar",
            expectedSize: 27926960,
            expectedChecksum: "ca39f25810a0943caec07fe6b8c1eb3e",
            namePrefix: "firebase-database-emulator",
        },
    },
    firestore: {
        downloadPath: path.join(CACHE_DIR, "cloud-firestore-emulator-v1.11.1.jar"),
        version: "1.11.1",
        opts: {
            cacheDir: CACHE_DIR,
            remoteUrl: "https://storage.googleapis.com/firebase-preview-drop/emulator/cloud-firestore-emulator-v1.11.1.jar",
            expectedSize: 63439953,
            expectedChecksum: "aa9a62f7b586731ed7664ab42fd20038",
            namePrefix: "cloud-firestore-emulator",
        },
    },
    gui: {
        version: "0.0.1",
        downloadPath: path.join(CACHE_DIR, "gui-v0.0.1.zip"),
        unzipDir: path.join(CACHE_DIR, "gui-v0.0.1"),
        binaryPath: path.join(CACHE_DIR, "gui-v0.0.1", `server.bundle.js`),
        opts: {
            cacheDir: CACHE_DIR,
            remoteUrl: "https://storage.googleapis.com/firebase-preview-drop/emulator/gui-v0.0.1.zip",
            expectedSize: 1204964,
            expectedChecksum: "52847b962bb66de639487d96a07a69d3",
            namePrefix: "gui",
        },
    },
    pubsub: {
        downloadPath: path.join(CACHE_DIR, "pubsub-emulator-0.1.0.zip"),
        version: "0.1.0",
        unzipDir: path.join(CACHE_DIR, "pubsub-emulator-0.1.0"),
        binaryPath: path.join(CACHE_DIR, "pubsub-emulator-0.1.0", `pubsub-emulator/bin/cloud-pubsub-emulator${process.platform === "win32" ? ".bat" : ""}`),
        opts: {
            cacheDir: CACHE_DIR,
            remoteUrl: "https://storage.googleapis.com/firebase-preview-drop/emulator/pubsub-emulator-0.1.0.zip",
            expectedSize: 36623622,
            expectedChecksum: "81704b24737d4968734d3e175f4cde71",
            namePrefix: "pubsub-emulator",
        },
    },
};
const EmulatorDetails = {
    database: {
        name: types_1.Emulators.DATABASE,
        instance: null,
        stdout: null,
    },
    firestore: {
        name: types_1.Emulators.FIRESTORE,
        instance: null,
        stdout: null,
    },
    pubsub: {
        name: types_1.Emulators.PUBSUB,
        instance: null,
        stdout: null,
    },
    gui: {
        name: types_1.Emulators.GUI,
        instance: null,
        stdout: null,
    },
};
const Commands = {
    database: {
        binary: "java",
        args: ["-Duser.language=en", "-jar", getExecPath(types_1.Emulators.DATABASE)],
        optionalArgs: ["port", "host", "functions_emulator_port", "functions_emulator_host"],
        joinArgs: false,
    },
    firestore: {
        binary: "java",
        args: ["-Duser.language=en", "-jar", getExecPath(types_1.Emulators.FIRESTORE)],
        optionalArgs: [
            "port",
            "webchannel_port",
            "host",
            "rules",
            "functions_emulator",
            "seed_from_export",
        ],
        joinArgs: false,
    },
    pubsub: {
        binary: getExecPath(types_1.Emulators.PUBSUB),
        args: [],
        optionalArgs: ["port", "host"],
        joinArgs: true,
    },
    gui: {
        binary: "node",
        args: [getExecPath(types_1.Emulators.GUI)],
        optionalArgs: [],
        joinArgs: false,
    },
};
function getExecPath(name) {
    const details = getDownloadDetails(name);
    return details.binaryPath || details.downloadPath;
}
function _getLogFileName(name) {
    return `${name}-debug.log`;
}
function _getCommand(emulator, args) {
    const baseCmd = Commands[emulator];
    const defaultPort = constants_1.Constants.getDefaultPort(emulator);
    if (!args.port) {
        args.port = defaultPort;
    }
    const cmdLineArgs = baseCmd.args.slice();
    Object.keys(args).forEach((key) => {
        if (baseCmd.optionalArgs.indexOf(key) < 0) {
            logger.debug(`Ignoring unsupported arg: ${key}`);
            return;
        }
        const argKey = "--" + key;
        const argVal = args[key];
        if (argVal === undefined) {
            logger.debug(`Ignoring empty arg for key: ${key}`);
            return;
        }
        if (baseCmd.joinArgs) {
            cmdLineArgs.push(`${argKey}=${argVal}`);
        }
        else {
            cmdLineArgs.push(argKey, argVal);
        }
    });
    return {
        binary: baseCmd.binary,
        args: cmdLineArgs,
        optionalArgs: baseCmd.optionalArgs,
        joinArgs: baseCmd.joinArgs,
    };
}
function _fatal(emulator, errorMsg) {
    if (emulator.instance) {
        emulator.instance.kill("SIGINT");
    }
    throw new error_1.FirebaseError(emulator.name + ": " + errorMsg, { exit: 1 });
}
function _runBinary(emulator, command, extraEnv) {
    return __awaiter(this, void 0, void 0, function* () {
        return new Promise((resolve) => {
            emulator.stdout = fs.createWriteStream(_getLogFileName(emulator.name));
            try {
                emulator.instance = childProcess.spawn(command.binary, command.args, {
                    env: Object.assign(Object.assign({}, process.env), extraEnv),
                    stdio: ["inherit", "pipe", "pipe"],
                });
            }
            catch (e) {
                if (e.code === "EACCES") {
                    utils.logLabeledWarning(emulator.name, `Could not spawn child process for emulator, check that java is installed and on your $PATH.`);
                }
                _fatal(emulator, e);
            }
            const description = constants_1.Constants.description(emulator.name);
            if (emulator.instance == null) {
                utils.logLabeledWarning(emulator.name, `Could not spawn child process for ${description}.`);
                return;
            }
            utils.logLabeledBullet(emulator.name, `${description} logging to ${clc.bold(_getLogFileName(emulator.name))}`);
            emulator.instance.stdout.on("data", (data) => {
                logger.debug(data.toString());
                emulator.stdout.write(data);
            });
            emulator.instance.stderr.on("data", (data) => {
                logger.debug(data.toString());
                emulator.stdout.write(data);
            });
            emulator.instance.on("error", (err) => {
                if (err.path === "java" && err.code === "ENOENT") {
                    _fatal(emulator, `${description} has exited because java is not installed, you can install it from https://openjdk.java.net/install/`);
                }
                else {
                    _fatal(emulator, `${description} has exited: ${err}`);
                }
            });
            emulator.instance.once("exit", (code, signal) => {
                if (signal) {
                    utils.logWarning(`${description} has exited upon receiving signal: ${signal}`);
                }
                else if (code && code !== 0 && code !== 130) {
                    _fatal(emulator, `${description} has exited with code: ${code}`);
                }
            });
            resolve();
        });
    });
}
function getDownloadDetails(emulator) {
    return DownloadDetails[emulator];
}
exports.getDownloadDetails = getDownloadDetails;
function get(emulator) {
    return EmulatorDetails[emulator];
}
exports.get = get;
function stop(targetName) {
    return __awaiter(this, void 0, void 0, function* () {
        const emulator = EmulatorDetails[targetName];
        return new Promise((resolve, reject) => {
            if (emulator.instance) {
                const killTimeout = setTimeout(() => {
                    const pid = emulator.instance ? emulator.instance.pid : -1;
                    const errorMsg = constants_1.Constants.description(emulator.name) + ": Unable to terminate process (PID=" + pid + ")";
                    logger.debug(errorMsg);
                    reject(new error_1.FirebaseError(emulator.name + ": " + errorMsg));
                }, EMULATOR_INSTANCE_KILL_TIMEOUT);
                emulator.instance.once("exit", () => {
                    clearTimeout(killTimeout);
                    resolve();
                });
                emulator.instance.kill("SIGINT");
            }
            else {
                resolve();
            }
        });
    });
}
exports.stop = stop;
function downloadIfNecessary(targetName) {
    return __awaiter(this, void 0, void 0, function* () {
        const hasEmulator = fs.existsSync(getExecPath(targetName));
        if (hasEmulator) {
            return;
        }
        yield downloadEmulator(targetName);
    });
}
exports.downloadIfNecessary = downloadIfNecessary;
function start(targetName, args, extraEnv = {}) {
    return __awaiter(this, void 0, void 0, function* () {
        const downloadDetails = DownloadDetails[targetName];
        const emulator = EmulatorDetails[targetName];
        const hasEmulator = fs.existsSync(getExecPath(targetName));
        if (!hasEmulator) {
            if (args.auto_download) {
                if (process.env.CI) {
                    utils.logWarning(`It appears you are running in a CI environment. You can avoid downloading the ${constants_1.Constants.description(targetName)} repeatedly by caching the ${downloadDetails.opts.cacheDir} directory.`);
                }
                yield downloadEmulator(targetName);
            }
            else {
                utils.logWarning("Setup required, please run: firebase setup:emulators:" + targetName);
                throw new error_1.FirebaseError("emulator not found");
            }
        }
        const command = _getCommand(targetName, args);
        logger.debug(`Starting ${constants_1.Constants.description(targetName)} with command ${JSON.stringify(command)}`);
        return _runBinary(emulator, command, extraEnv);
    });
}
exports.start = start;
