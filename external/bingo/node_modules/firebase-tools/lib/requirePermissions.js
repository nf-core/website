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
const cli_color_1 = require("cli-color");
const api_1 = require("./api");
const getProjectId = require("./getProjectId");
const requireAuth = require("./requireAuth");
const logger_1 = require("./logger");
const error_1 = require("./error");
const BASE_PERMISSIONS = ["firebase.projects.get"];
function requirePermissions(options, permissions = []) {
    return __awaiter(this, void 0, void 0, function* () {
        const projectId = getProjectId(options);
        const requiredPermissions = BASE_PERMISSIONS.concat(permissions).sort();
        yield requireAuth(options);
        logger_1.debug(`[iam] checking project ${projectId} for permissions ${JSON.stringify(requiredPermissions)}`);
        let response;
        try {
            response = yield api_1.request("POST", `/v1/projects/${projectId}:testIamPermissions`, {
                auth: true,
                data: { permissions: requiredPermissions },
                origin: api_1.resourceManagerOrigin,
            });
        }
        catch (err) {
            logger_1.debug(`[iam] error while checking permissions, command may fail: ${err}`);
            return;
        }
        const allowedPermissions = (response.body.permissions || []).sort();
        const missingPermissions = lodash_1.difference(requiredPermissions, allowedPermissions);
        if (missingPermissions.length) {
            throw new error_1.FirebaseError(`Authorization failed. This account is missing the following required permissions on project ${cli_color_1.bold(projectId)}:\n\n  ${missingPermissions.join("\n  ")}`);
        }
    });
}
exports.requirePermissions = requirePermissions;
