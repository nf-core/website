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
const operationPoller = require("../operation-poller");
const VERSION = "v1beta";
var ParamType;
(function (ParamType) {
    ParamType["STRING"] = "STRING";
    ParamType["SELECT"] = "SELECT";
    ParamType["MULTISELECT"] = "MULTISELECT";
})(ParamType = exports.ParamType || (exports.ParamType = {}));
function createInstance(projectId, instanceId, extensionSource, params, serviceAccountEmail) {
    return __awaiter(this, void 0, void 0, function* () {
        const createRes = yield api.request("POST", `/${VERSION}/projects/${projectId}/instances/`, {
            auth: true,
            origin: api.extensionsOrigin,
            data: {
                name: `projects/${projectId}/instances/${instanceId}`,
                config: {
                    source: { name: extensionSource.name },
                    params,
                },
                serviceAccountEmail,
            },
        });
        const pollRes = yield operationPoller.pollOperation({
            apiOrigin: api.extensionsOrigin,
            apiVersion: VERSION,
            operationResourceName: createRes.body.name,
            masterTimeout: 600000,
        });
        return pollRes;
    });
}
exports.createInstance = createInstance;
function deleteInstance(projectId, instanceId) {
    return __awaiter(this, void 0, void 0, function* () {
        const deleteRes = yield api.request("DELETE", `/${VERSION}/projects/${projectId}/instances/${instanceId}`, {
            auth: true,
            origin: api.extensionsOrigin,
        });
        const pollRes = yield operationPoller.pollOperation({
            apiOrigin: api.extensionsOrigin,
            apiVersion: VERSION,
            operationResourceName: deleteRes.body.name,
            masterTimeout: 600000,
        });
        return pollRes;
    });
}
exports.deleteInstance = deleteInstance;
function getInstance(projectId, instanceId, options = {}) {
    return __awaiter(this, void 0, void 0, function* () {
        const res = yield api.request("GET", `/${VERSION}/projects/${projectId}/instances/${instanceId}`, _.assign({
            auth: true,
            origin: api.extensionsOrigin,
        }, options));
        return res.body;
    });
}
exports.getInstance = getInstance;
function listInstances(projectId) {
    return __awaiter(this, void 0, void 0, function* () {
        const instances = [];
        const getNextPage = (pageToken) => __awaiter(this, void 0, void 0, function* () {
            const res = yield api.request("GET", `/${VERSION}/projects/${projectId}/instances`, {
                auth: true,
                origin: api.extensionsOrigin,
                query: {
                    pageSize: 100,
                    pageToken,
                },
            });
            instances.push(...res.body.instances);
            if (res.body.nextPageToken) {
                yield getNextPage(res.body.nextPageToken);
            }
        });
        yield getNextPage();
        return instances;
    });
}
exports.listInstances = listInstances;
function configureInstance(projectId, instanceId, params) {
    return __awaiter(this, void 0, void 0, function* () {
        const res = yield patchInstance(projectId, instanceId, "config.params", {
            config: {
                params,
            },
        });
        return res;
    });
}
exports.configureInstance = configureInstance;
function updateInstance(projectId, instanceId, extensionSource, params) {
    return __awaiter(this, void 0, void 0, function* () {
        const body = {
            config: {
                source: { name: extensionSource.name },
            },
        };
        let updateMask = "config.source.name";
        if (params) {
            body.params = params;
            updateMask += ",config.params";
        }
        return yield patchInstance(projectId, instanceId, updateMask, body);
    });
}
exports.updateInstance = updateInstance;
function patchInstance(projectId, instanceId, updateMask, data) {
    return __awaiter(this, void 0, void 0, function* () {
        const updateRes = yield api.request("PATCH", `/${VERSION}/projects/${projectId}/instances/${instanceId}`, {
            auth: true,
            origin: api.extensionsOrigin,
            query: {
                updateMask,
            },
            data,
        });
        const pollRes = yield operationPoller.pollOperation({
            apiOrigin: api.extensionsOrigin,
            apiVersion: VERSION,
            operationResourceName: updateRes.body.name,
            masterTimeout: 600000,
        });
        return pollRes;
    });
}
function createSource(projectId, packageUri, extensionRoot) {
    return __awaiter(this, void 0, void 0, function* () {
        const createRes = yield api.request("POST", `/${VERSION}/projects/${projectId}/sources/`, {
            auth: true,
            origin: api.extensionsOrigin,
            data: {
                packageUri,
                extensionRoot,
            },
        });
        const pollRes = yield operationPoller.pollOperation({
            apiOrigin: api.extensionsOrigin,
            apiVersion: VERSION,
            operationResourceName: createRes.body.name,
            masterTimeout: 600000,
        });
        return pollRes;
    });
}
exports.createSource = createSource;
function getSource(sourceName) {
    return api
        .request("GET", `/${VERSION}/${sourceName}`, {
        auth: true,
        origin: api.extensionsOrigin,
    })
        .then((res) => {
        return res.body;
    });
}
exports.getSource = getSource;
