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
const api = require("../api");
const utils = require("../utils");
const gcp_1 = require("../gcp");
const utils_1 = require("./utils");
const API_VERSION = "v1";
function grantRoles(projectId, serviceAccountEmail, rolesToAdd, rolesToRemove = []) {
    rolesToAdd = rolesToAdd.map((role) => `roles/${role}`);
    rolesToRemove = rolesToRemove.map((role) => `roles/${role}`);
    return api
        .request("POST", utils.endpoint([API_VERSION, "projects", projectId, ":getIamPolicy"]), {
        data: {
            options: { requestedPolicyVersion: 3 },
        },
        auth: true,
        origin: api.resourceManagerOrigin,
    })
        .then((response) => {
        const policy = response.body;
        const bindings = policy.bindings;
        rolesToAdd.forEach((role) => {
            bindings.push({ role, members: [`serviceAccount:${serviceAccountEmail}`] });
        });
        rolesToRemove.forEach((role) => {
            const binding = _.find(bindings, (b) => {
                return b.role === role;
            });
            _.remove(binding.members, (member) => {
                return member === `serviceAccount:${serviceAccountEmail}`;
            });
        });
        return api.request("POST", utils.endpoint([API_VERSION, "projects", projectId, ":setIamPolicy"]), {
            auth: true,
            origin: api.resourceManagerOrigin,
            data: { policy },
        });
    });
}
exports.grantRoles = grantRoles;
function createServiceAccountAndSetRoles(projectId, roles, instanceId) {
    return __awaiter(this, void 0, void 0, function* () {
        let serviceAccount;
        const shortenedInstanceId = instanceId.length <= 26 ? instanceId : `${instanceId.slice(0, 21)}-${utils_1.getRandomString(4)}`;
        try {
            serviceAccount = yield gcp_1.iam.createServiceAccount(projectId, `ext-${shortenedInstanceId}`, `Runtime service account for Firebase Extension ${instanceId}`, `Firebase Extension ${instanceId} service account`);
        }
        catch (err) {
            if (err.status === 409) {
                return utils.reject(`A service account ext-${shortenedInstanceId} already exists in project ${projectId}. ` +
                    `Please delete it or choose a different extension instance id.`, {
                    exit: 1,
                    status: 409,
                });
            }
            throw err;
        }
        yield grantRoles(projectId, serviceAccount.email, roles.map((role) => role.role));
        return serviceAccount.email;
    });
}
exports.createServiceAccountAndSetRoles = createServiceAccountAndSetRoles;
