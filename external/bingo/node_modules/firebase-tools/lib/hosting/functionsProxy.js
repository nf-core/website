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
const lodash_1 = require("lodash");
const proxy_1 = require("./proxy");
const getProjectId = require("../getProjectId");
const registry_1 = require("../emulator/registry");
const types_1 = require("../emulator/types");
const functionsEmulator_1 = require("../emulator/functionsEmulator");
function default_1(options) {
    return (rewrite) => __awaiter(this, void 0, void 0, function* () {
        const projectId = getProjectId(options, false);
        let url = `https://us-central1-${projectId}.cloudfunctions.net/${rewrite.function}`;
        let destLabel = "live";
        if (lodash_1.includes(options.targets, "functions")) {
            destLabel = "local";
            const functionsEmu = registry_1.EmulatorRegistry.get(types_1.Emulators.FUNCTIONS);
            if (functionsEmu) {
                url = functionsEmulator_1.FunctionsEmulator.getHttpFunctionUrl(functionsEmu.getInfo().host, functionsEmu.getInfo().port, projectId, rewrite.function, "us-central1");
            }
        }
        return yield proxy_1.proxyRequestHandler(url, `${destLabel} Function ${rewrite.function}`);
    });
}
exports.default = default_1;
