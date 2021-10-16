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
const serveHosting = require("../serve/hosting");
const types_1 = require("../emulator/types");
const constants_1 = require("./constants");
class HostingEmulator {
    constructor(args) {
        this.args = args;
    }
    start() {
        return __awaiter(this, void 0, void 0, function* () {
            this.args.options.host = this.args.host;
            this.args.options.port = this.args.port;
            return serveHosting.start(this.args.options);
        });
    }
    connect() {
        return __awaiter(this, void 0, void 0, function* () {
            return;
        });
    }
    stop() {
        return __awaiter(this, void 0, void 0, function* () {
            return serveHosting.stop();
        });
    }
    getInfo() {
        const host = this.args.host || constants_1.Constants.getDefaultHost(types_1.Emulators.HOSTING);
        const port = this.args.port || constants_1.Constants.getDefaultPort(types_1.Emulators.HOSTING);
        return {
            host,
            port,
        };
    }
    getName() {
        return types_1.Emulators.HOSTING;
    }
}
exports.HostingEmulator = HostingEmulator;
