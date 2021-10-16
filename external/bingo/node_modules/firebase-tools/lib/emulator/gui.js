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
const downloadableEmulators = require("./downloadableEmulators");
const registry_1 = require("./registry");
const hub_1 = require("./hub");
const error_1 = require("../error");
const constants_1 = require("./constants");
class EmulatorGUI {
    constructor(args) {
        this.args = args;
    }
    start() {
        if (!registry_1.EmulatorRegistry.isRunning(types_1.Emulators.HUB)) {
            throw new error_1.FirebaseError(`Cannot start ${constants_1.Constants.description(types_1.Emulators.GUI)} without ${constants_1.Constants.description(types_1.Emulators.HUB)}!`);
        }
        const hubInfo = registry_1.EmulatorRegistry.get(types_1.Emulators.HUB).getInfo();
        const { auto_download, host, port, projectId } = this.args;
        const env = {
            HOST: host.toString(),
            PORT: port.toString(),
            GCLOUD_PROJECT: projectId,
            [hub_1.EmulatorHub.EMULATOR_HUB_ENV]: `${hubInfo.host}:${hubInfo.port}`,
        };
        return downloadableEmulators.start(types_1.Emulators.GUI, { auto_download }, env);
    }
    connect() {
        return __awaiter(this, void 0, void 0, function* () {
            return;
        });
    }
    stop() {
        return downloadableEmulators.stop(types_1.Emulators.GUI);
    }
    getInfo() {
        return {
            host: this.args.host,
            port: this.args.port,
        };
    }
    getName() {
        return types_1.Emulators.GUI;
    }
}
exports.EmulatorGUI = EmulatorGUI;
