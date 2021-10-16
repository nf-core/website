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
const express = require("express");
const os = require("os");
const fs = require("fs");
const path = require("path");
const bodyParser = require("body-parser");
const utils = require("../utils");
const logger = require("../logger");
const constants_1 = require("./constants");
const types_1 = require("./types");
const hubExport_1 = require("./hubExport");
const registry_1 = require("./registry");
const pkg = require("../../package.json");
class EmulatorHub {
    constructor(args) {
        this.args = args;
        this.hub = express();
        this.hub.use(bodyParser.json());
        this.hub.get("/", (req, res) => __awaiter(this, void 0, void 0, function* () {
            res.json(this.getLocator());
        }));
        this.hub.get(EmulatorHub.PATH_EMULATORS, (req, res) => __awaiter(this, void 0, void 0, function* () {
            const body = {};
            registry_1.EmulatorRegistry.listRunning().forEach((name) => {
                body[name] = registry_1.EmulatorRegistry.get(name).getInfo();
            });
            res.json(body);
        }));
        this.hub.post(EmulatorHub.PATH_EXPORT, (req, res) => __awaiter(this, void 0, void 0, function* () {
            const exportPath = req.body.path;
            utils.logLabeledBullet("emulators", `Received export request. Exporting data to ${exportPath}.`);
            try {
                yield new hubExport_1.HubExport(this.args.projectId, exportPath).exportAll();
                utils.logLabeledSuccess("emulators", "Export complete.");
                res.status(200).send({
                    message: "OK",
                });
            }
            catch (e) {
                const errorString = e.message || JSON.stringify(e);
                utils.logLabeledWarning("emulators", `Export failed: ${errorString}`);
                res.status(500).json({
                    message: errorString,
                });
            }
        }));
    }
    static readLocatorFile(projectId) {
        const locatorPath = this.getLocatorFilePath(projectId);
        if (!fs.existsSync(locatorPath)) {
            return undefined;
        }
        const data = fs.readFileSync(locatorPath, "utf8").toString();
        const locator = JSON.parse(data);
        if (locator.version !== this.CLI_VERSION) {
            logger.debug(`Found locator with mismatched version, ignoring: ${JSON.stringify(locator)}`);
            return undefined;
        }
        return locator;
    }
    static getLocatorFilePath(projectId) {
        const dir = os.tmpdir();
        const filename = `hub-${projectId}.json`;
        return path.join(dir, filename);
    }
    start() {
        return __awaiter(this, void 0, void 0, function* () {
            const { host, port } = this.getInfo();
            this.server = this.hub.listen(port, host);
            yield this.writeLocatorFile();
        });
    }
    connect() {
        return __awaiter(this, void 0, void 0, function* () {
        });
    }
    stop() {
        return __awaiter(this, void 0, void 0, function* () {
            this.server && this.server.close();
            yield this.deleteLocatorFile();
        });
    }
    getInfo() {
        const host = this.args.host || constants_1.Constants.getDefaultHost(types_1.Emulators.HUB);
        const port = this.args.port || constants_1.Constants.getDefaultPort(types_1.Emulators.HUB);
        return {
            host,
            port,
        };
    }
    getName() {
        return types_1.Emulators.HUB;
    }
    getLocator() {
        const { host, port } = this.getInfo();
        const version = pkg.version;
        return {
            version,
            host,
            port,
        };
    }
    writeLocatorFile() {
        return __awaiter(this, void 0, void 0, function* () {
            const projectId = this.args.projectId;
            const locatorPath = EmulatorHub.getLocatorFilePath(projectId);
            const locator = this.getLocator();
            if (fs.existsSync(locatorPath)) {
                utils.logLabeledWarning("emulators", `It seems that you are running multiple instances of the emulator suite for project ${projectId}. This may result in unexpected behavior.`);
            }
            logger.debug(`[hub] writing locator at ${locatorPath}`);
            return new Promise((resolve, reject) => {
                fs.writeFile(locatorPath, JSON.stringify(locator), (e) => {
                    if (e) {
                        reject(e);
                    }
                    else {
                        resolve();
                    }
                });
            });
        });
    }
    deleteLocatorFile() {
        return __awaiter(this, void 0, void 0, function* () {
            const locatorPath = EmulatorHub.getLocatorFilePath(this.args.projectId);
            return new Promise((resolve, reject) => {
                fs.unlink(locatorPath, (e) => {
                    if (e) {
                        reject(e);
                    }
                    else {
                        resolve();
                    }
                });
            });
        });
    }
}
exports.EmulatorHub = EmulatorHub;
EmulatorHub.EMULATOR_HUB_ENV = "FIREBASE_EMULATOR_HUB";
EmulatorHub.CLI_VERSION = pkg.version;
EmulatorHub.PATH_EXPORT = "/_admin/export";
EmulatorHub.PATH_EMULATORS = "/emulators";
