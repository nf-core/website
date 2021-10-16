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
function handleErrorResponse(response) {
    if (response.body && response.body.error) {
        return utils.reject(response.body.error, { code: 2 });
    }
    logger.debug("[rules] error:", response.status, response.body);
    return utils.reject("Unexpected error encountered with database.", {
        code: 2,
    });
}
function listAllRulesets(databaseName) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("GET", `/namespaces/${databaseName}/rulesets`, {
            auth: true,
            origin: api.rtdbMetadataOrigin,
            json: true,
        });
        if (response.status === 200) {
            return response.body.rulesets;
        }
        return handleErrorResponse(response);
    });
}
exports.listAllRulesets = listAllRulesets;
function getRuleset(databaseName, rulesetId) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("GET", `/namespaces/${databaseName}/rulesets/${rulesetId}`, {
            auth: true,
            origin: api.rtdbMetadataOrigin,
            json: true,
        });
        if (response.status === 200) {
            return response.body;
        }
        return handleErrorResponse(response);
    });
}
exports.getRuleset = getRuleset;
function getRulesetLabels(databaseName) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("GET", `/namespaces/${databaseName}/ruleset_labels`, {
            auth: true,
            origin: api.rtdbMetadataOrigin,
        });
        if (response.status === 200) {
            return response.body;
        }
        return handleErrorResponse(response);
    });
}
exports.getRulesetLabels = getRulesetLabels;
function createRuleset(databaseName, source) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("POST", `/.settings/rulesets.json`, {
            auth: true,
            origin: utils.addSubdomain(api.realtimeOrigin, databaseName),
            json: false,
            data: source,
        });
        if (response.status === 200) {
            return JSON.parse(response.body).id;
        }
        return handleErrorResponse(response);
    });
}
exports.createRuleset = createRuleset;
function setRulesetLabels(databaseName, labels) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("PUT", `/.settings/ruleset_labels.json`, {
            auth: true,
            origin: utils.addSubdomain(api.realtimeOrigin, databaseName),
            data: labels,
        });
        if (response.status === 200) {
            return response.body;
        }
        return handleErrorResponse(response);
    });
}
exports.setRulesetLabels = setRulesetLabels;
