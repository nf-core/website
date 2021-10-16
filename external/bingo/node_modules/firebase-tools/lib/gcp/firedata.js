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
const api = require("../api");
const logger = require("../logger");
const utils = require("../utils");
function _handleErrorResponse(response) {
    if (response.body && response.body.error) {
        return utils.reject(response.body.error, { code: 2 });
    }
    logger.debug("[firedata] error:", response.status, response.body);
    return utils.reject("Unexpected error encountered with FireData.", {
        code: 2,
    });
}
function createDatabaseInstance(projectNumber, instanceName) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("POST", `/v1/projects/${projectNumber}/databases`, {
            auth: true,
            origin: api.firedataOrigin,
            json: {
                instance: instanceName,
            },
        });
        if (response.status === 200) {
            return response.body.instance;
        }
        return _handleErrorResponse(response);
    });
}
exports.createDatabaseInstance = createDatabaseInstance;
function listDatabaseInstances(projectNumber) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("GET", `/v1/projects/${projectNumber}/databases`, {
            auth: true,
            origin: api.firedataOrigin,
        });
        if (response.status === 200) {
            return response.body.instance;
        }
        return _handleErrorResponse(response);
    });
}
exports.listDatabaseInstances = listDatabaseInstances;
