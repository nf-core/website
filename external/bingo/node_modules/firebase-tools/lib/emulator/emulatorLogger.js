"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const clc = require("cli-color");
const utils = require("../utils");
const logger = require("../logger");
const TYPE_VERBOSITY = {
    DEBUG: 0,
    INFO: 1,
    BULLET: 1,
    SUCCESS: 1,
    USER: 2,
    WARN: 2,
    WARN_ONCE: 2,
};
var Verbosity;
(function (Verbosity) {
    Verbosity[Verbosity["DEBUG"] = 0] = "DEBUG";
    Verbosity[Verbosity["INFO"] = 1] = "INFO";
    Verbosity[Verbosity["QUIET"] = 2] = "QUIET";
})(Verbosity = exports.Verbosity || (exports.Verbosity = {}));
class EmulatorLogger {
    static log(type, text) {
        if (EmulatorLogger.shouldSupress(type)) {
            logger.debug(`${type}: ${text}`);
            return;
        }
        switch (type) {
            case "DEBUG":
                logger.debug(text);
                break;
            case "INFO":
                logger.info(text);
                break;
            case "USER":
                logger.info(text);
                break;
            case "BULLET":
                utils.logBullet(text);
                break;
            case "WARN":
                utils.logWarning(text);
                break;
            case "WARN_ONCE":
                if (!this.warnOnceCache.has(text)) {
                    utils.logWarning(text);
                    this.warnOnceCache.add(text);
                }
                break;
            case "SUCCESS":
                utils.logSuccess(text);
                break;
        }
    }
    static handleRuntimeLog(log, ignore = []) {
        if (ignore.indexOf(log.level) >= 0) {
            return;
        }
        switch (log.level) {
            case "SYSTEM":
                EmulatorLogger.handleSystemLog(log);
                break;
            case "USER":
                EmulatorLogger.log("USER", `${clc.blackBright("> ")} ${log.text}`);
                break;
            case "DEBUG":
                if (log.data && Object.keys(log.data).length > 0) {
                    EmulatorLogger.log("DEBUG", `[${log.type}] ${log.text} ${JSON.stringify(log.data)}`);
                }
                else {
                    EmulatorLogger.log("DEBUG", `[${log.type}] ${log.text}`);
                }
                break;
            case "INFO":
                EmulatorLogger.logLabeled("BULLET", "functions", log.text);
                break;
            case "WARN":
                EmulatorLogger.logLabeled("WARN", "functions", log.text);
                break;
            case "WARN_ONCE":
                EmulatorLogger.logLabeled("WARN_ONCE", "functions", log.text);
                break;
            case "FATAL":
                EmulatorLogger.logLabeled("WARN", "functions", log.text);
                break;
            default:
                EmulatorLogger.log("INFO", `${log.level}: ${log.text}`);
                break;
        }
    }
    static handleSystemLog(systemLog) {
        switch (systemLog.type) {
            case "runtime-status":
                if (systemLog.text === "killed") {
                    EmulatorLogger.log("WARN", `Your function was killed because it raised an unhandled error.`);
                }
                break;
            case "googleapis-network-access":
                EmulatorLogger.log("WARN", `Google API requested!\n   - URL: "${systemLog.data.href}"\n   - Be careful, this may be a production service.`);
                break;
            case "unidentified-network-access":
                EmulatorLogger.log("WARN", `External network resource requested!\n   - URL: "${systemLog.data.href}"\n - Be careful, this may be a production service.`);
                break;
            case "functions-config-missing-value":
                EmulatorLogger.log("WARN", `Non-existent functions.config() value requested!\n   - Path: "${systemLog.data.valuePath}"\n   - Learn more at https://firebase.google.com/docs/functions/local-emulator`);
                break;
            case "non-default-admin-app-used":
                EmulatorLogger.log("WARN", `Non-default "firebase-admin" instance created!\n   ` +
                    `- This instance will *not* be mocked and will access production resources.`);
                break;
            case "missing-module":
                EmulatorLogger.log("WARN", `The Cloud Functions emulator requires the module "${systemLog.data.name}" to be installed as a ${systemLog.data.isDev ? "development dependency" : "dependency"}. To fix this, run "npm install ${systemLog.data.isDev ? "--save-dev" : "--save"} ${systemLog.data.name}" in your functions directory.`);
                break;
            case "uninstalled-module":
                EmulatorLogger.log("WARN", `The Cloud Functions emulator requires the module "${systemLog.data.name}" to be installed. This package is in your package.json, but it's not available. \
You probably need to run "npm install" in your functions directory.`);
                break;
            case "out-of-date-module":
                EmulatorLogger.log("WARN", `The Cloud Functions emulator requires the module "${systemLog.data.name}" to be version >${systemLog.data.minVersion} so your version is too old. \
You can probably fix this by running "npm install ${systemLog.data.name}@latest" in your functions directory.`);
                break;
            case "missing-package-json":
                EmulatorLogger.log("WARN", `The Cloud Functions directory you specified does not have a "package.json" file, so we can't load it.`);
                break;
            case "function-code-resolution-failed":
                EmulatorLogger.log("WARN", systemLog.data.error);
                const helper = ["We were unable to load your functions code. (see above)"];
                if (systemLog.data.isPotentially.wrong_directory) {
                    helper.push(`   - There is no "package.json" file in your functions directory.`);
                }
                if (systemLog.data.isPotentially.typescript) {
                    helper.push("   - It appears your code is written in Typescript, which must be compiled before emulation.");
                }
                if (systemLog.data.isPotentially.uncompiled) {
                    helper.push(`   - You may be able to run "npm run build" in your functions directory to resolve this.`);
                }
                utils.logWarning(helper.join("\n"));
            default:
        }
    }
    static logLabeled(type, label, text) {
        if (EmulatorLogger.shouldSupress(type)) {
            logger.debug(`[${label}] ${text}`);
            return;
        }
        switch (type) {
            case "DEBUG":
                logger.debug(`[${label}] ${text}`);
                break;
            case "BULLET":
                utils.logLabeledBullet(label, text);
                break;
            case "SUCCESS":
                utils.logLabeledSuccess(label, text);
                break;
            case "WARN":
                utils.logLabeledWarning(label, text);
                break;
            case "WARN_ONCE":
                if (!this.warnOnceCache.has(text)) {
                    utils.logLabeledWarning(label, text);
                    this.warnOnceCache.add(text);
                }
                break;
        }
    }
    static shouldSupress(type) {
        const typeVerbosity = TYPE_VERBOSITY[type];
        return EmulatorLogger.verbosity > typeVerbosity;
    }
}
exports.EmulatorLogger = EmulatorLogger;
EmulatorLogger.verbosity = Verbosity.DEBUG;
EmulatorLogger.warnOnceCache = new Set();
