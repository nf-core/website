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
const path = require("path");
const functionsEmulator_1 = require("../emulator/functionsEmulator");
const emulatorServer_1 = require("../emulator/emulatorServer");
const getProjectId = require("../getProjectId");
module.exports = {
    emulatorServer: undefined,
    start(options, args) {
        return __awaiter(this, void 0, void 0, function* () {
            const projectId = getProjectId(options, false);
            const functionsDir = path.join(options.config.projectDir, options.config.get("functions.source"));
            args = Object.assign({ projectId,
                functionsDir }, args);
            if (!args.disabledRuntimeFeatures) {
                args.disabledRuntimeFeatures = {
                    functions_config_helper: true,
                    network_filtering: true,
                    timeout: true,
                    memory_limiting: true,
                    admin_stubs: true,
                    pubsub_emulator: true,
                };
            }
            if (options.host) {
                args.host = options.host;
            }
            if (options.port) {
                const hostingRunning = options.targets && options.targets.indexOf("hosting") >= 0;
                if (hostingRunning) {
                    args.port = options.port + 1;
                }
                else {
                    args.port = options.port;
                }
            }
            this.emulatorServer = new emulatorServer_1.EmulatorServer(new functionsEmulator_1.FunctionsEmulator(args));
            yield this.emulatorServer.start();
        });
    },
    connect() {
        return __awaiter(this, void 0, void 0, function* () {
            yield this.emulatorServer.connect();
        });
    },
    stop() {
        return __awaiter(this, void 0, void 0, function* () {
            yield this.emulatorServer.stop();
        });
    },
    get() {
        return this.emulatorServer.get();
    },
};
