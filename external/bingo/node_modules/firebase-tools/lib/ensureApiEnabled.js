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
const _ = require("lodash");
const api = require("./api");
const utils = require("./utils");
const error_1 = require("./error");
const POLL_INTERVAL = 10000;
const POLLS_BEFORE_RETRY = 12;
function check(projectId, apiName, prefix, silent = false) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("GET", `/v1/projects/${projectId}/services/${apiName}`, {
            auth: true,
            origin: api.serviceUsageOrigin,
        });
        const isEnabled = _.get(response.body, "state") === "ENABLED";
        if (isEnabled && !silent) {
            utils.logLabeledSuccess(prefix, "all necessary APIs are enabled");
        }
        return isEnabled;
    });
}
exports.check = check;
function enable(projectId, apiName) {
    return __awaiter(this, void 0, void 0, function* () {
        return api.request("POST", `/v1/projects/${projectId}/services/${apiName}:enable`, {
            auth: true,
            origin: api.serviceUsageOrigin,
        });
    });
}
exports.enable = enable;
function ensure(projectId, apiName, prefix, silent = false) {
    return __awaiter(this, void 0, void 0, function* () {
        if (!silent) {
            utils.logLabeledBullet(prefix, "ensuring necessary APIs are enabled...");
        }
        const isEnabled = yield check(projectId, apiName, prefix, silent);
        if (isEnabled) {
            return;
        }
        if (!silent) {
            utils.logLabeledWarning(prefix, "missing necessary APIs. Enabling now...");
        }
        return enableApiWithRetries(projectId, apiName, prefix, silent);
    });
}
exports.ensure = ensure;
function pollCheckEnabled(projectId, apiName, prefix, silent, enablementRetries, pollRetries = 0) {
    return __awaiter(this, void 0, void 0, function* () {
        if (pollRetries > POLLS_BEFORE_RETRY) {
            return enableApiWithRetries(projectId, apiName, prefix, silent, enablementRetries + 1);
        }
        yield new Promise((resolve) => {
            setTimeout(resolve, POLL_INTERVAL);
        });
        const isEnabled = yield check(projectId, apiName, prefix, silent);
        if (isEnabled) {
            return;
        }
        if (!silent) {
            utils.logLabeledBullet(prefix, "waiting for APIs to activate...");
        }
        return pollCheckEnabled(projectId, apiName, prefix, silent, enablementRetries, pollRetries + 1);
    });
}
function enableApiWithRetries(projectId, apiName, prefix, silent, enablementRetries = 0) {
    return __awaiter(this, void 0, void 0, function* () {
        if (enablementRetries > 1) {
            throw new error_1.FirebaseError("Timed out while waiting for APIs to enable. Please try again in a few minutes.");
        }
        yield enable(projectId, apiName);
        return pollCheckEnabled(projectId, apiName, prefix, silent, enablementRetries);
    });
}
