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
const api = require("./api");
const getProjectId = require("./getProjectId");
const configstore_1 = require("./configstore");
const CONFIGSTORE_KEY = "webconfig";
function setCachedWebSetup(projectId, config) {
    const allConfigs = configstore_1.configstore.get(CONFIGSTORE_KEY) || {};
    allConfigs[projectId] = config;
    configstore_1.configstore.set(CONFIGSTORE_KEY, allConfigs);
}
function getCachedWebSetup(options) {
    const projectId = getProjectId(options, false);
    const allConfigs = configstore_1.configstore.get(CONFIGSTORE_KEY) || {};
    return allConfigs[projectId];
}
exports.getCachedWebSetup = getCachedWebSetup;
function fetchWebSetup(options) {
    return __awaiter(this, void 0, void 0, function* () {
        const projectId = getProjectId(options, false);
        const response = yield api.request("GET", `/v1beta1/projects/${projectId}/webApps/-/config`, {
            auth: true,
            origin: api.firebaseApiOrigin,
        });
        const config = response.body;
        setCachedWebSetup(config.projectId, config);
        return config;
    });
}
exports.fetchWebSetup = fetchWebSetup;
