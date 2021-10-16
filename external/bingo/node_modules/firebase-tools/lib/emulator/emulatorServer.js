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
const registry_1 = require("./registry");
const controller = require("./controller");
const error_1 = require("../error");
class EmulatorServer {
    constructor(instance) {
        this.instance = instance;
    }
    start() {
        return __awaiter(this, void 0, void 0, function* () {
            const { port, host } = this.instance.getInfo();
            const portOpen = yield controller.checkPortOpen(port, host);
            if (!portOpen) {
                throw new error_1.FirebaseError(`Port ${port} is not open on ${host}, could not start ${this.instance.getName()} emulator.`);
            }
            yield registry_1.EmulatorRegistry.start(this.instance);
        });
    }
    connect() {
        return __awaiter(this, void 0, void 0, function* () {
            yield this.instance.connect();
        });
    }
    stop() {
        return __awaiter(this, void 0, void 0, function* () {
            yield registry_1.EmulatorRegistry.stop(this.instance.getName());
        });
    }
    get() {
        return this.instance;
    }
}
exports.EmulatorServer = EmulatorServer;
