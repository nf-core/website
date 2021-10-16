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
const logger = require("../logger");
const utils = require("../utils");
const API_VERSION = "v1";
function _handleErrorResponse(response) {
    if (response.body && response.body.error) {
        return utils.reject(response.body.error, { code: 2 });
    }
    logger.debug("[rules] error:", response.status, response.body);
    return utils.reject("Unexpected error encountered with rules.", {
        code: 2,
    });
}
function getLatestRulesetName(projectId, service) {
    return __awaiter(this, void 0, void 0, function* () {
        const releases = yield listAllReleases(projectId);
        const prefix = `projects/${projectId}/releases/${service}`;
        const release = _.find(releases, (r) => r.name.indexOf(prefix) === 0);
        if (!release) {
            return null;
        }
        return release.rulesetName;
    });
}
exports.getLatestRulesetName = getLatestRulesetName;
const MAX_RELEASES_PAGE_SIZE = 10;
function listReleases(projectId, pageToken) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("GET", `/${API_VERSION}/projects/${projectId}/releases`, {
            auth: true,
            origin: api.rulesOrigin,
            query: {
                pageSize: MAX_RELEASES_PAGE_SIZE,
                pageToken,
            },
        });
        if (response.status === 200) {
            return response.body;
        }
        return _handleErrorResponse(response);
    });
}
exports.listReleases = listReleases;
function listAllReleases(projectId) {
    return __awaiter(this, void 0, void 0, function* () {
        let pageToken;
        let releases = [];
        do {
            const response = yield listReleases(projectId, pageToken);
            if (response.releases && response.releases.length > 0) {
                releases = releases.concat(response.releases);
            }
            pageToken = response.nextPageToken;
        } while (pageToken);
        return _.orderBy(releases, ["createTime"], ["desc"]);
    });
}
exports.listAllReleases = listAllReleases;
function getRulesetContent(name) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("GET", `/${API_VERSION}/${name}`, {
            auth: true,
            origin: api.rulesOrigin,
        });
        if (response.status === 200) {
            const source = response.body.source;
            return source.files;
        }
        return _handleErrorResponse(response);
    });
}
exports.getRulesetContent = getRulesetContent;
const MAX_RULESET_PAGE_SIZE = 100;
function listRulesets(projectId, pageToken) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("GET", `/${API_VERSION}/projects/${projectId}/rulesets`, {
            auth: true,
            origin: api.rulesOrigin,
            query: {
                pageSize: MAX_RULESET_PAGE_SIZE,
                pageToken,
            },
        });
        if (response.status === 200) {
            return response.body;
        }
        return _handleErrorResponse(response);
    });
}
exports.listRulesets = listRulesets;
function listAllRulesets(projectId) {
    return __awaiter(this, void 0, void 0, function* () {
        let pageToken;
        let rulesets = [];
        do {
            const response = yield listRulesets(projectId, pageToken);
            if (response.rulesets) {
                rulesets = rulesets.concat(response.rulesets);
            }
            pageToken = response.nextPageToken;
        } while (pageToken);
        return _.orderBy(rulesets, ["createTime"], ["desc"]);
    });
}
exports.listAllRulesets = listAllRulesets;
function getRulesetId(ruleset) {
    return ruleset.name.split("/").pop();
}
exports.getRulesetId = getRulesetId;
function deleteRuleset(projectId, id) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield api.request("DELETE", `/${API_VERSION}/projects/${projectId}/rulesets/${id}`, {
            auth: true,
            origin: api.rulesOrigin,
        });
        if (response.status === 200) {
            return;
        }
        return _handleErrorResponse(response);
    });
}
exports.deleteRuleset = deleteRuleset;
function createRuleset(projectId, files) {
    return __awaiter(this, void 0, void 0, function* () {
        const payload = { source: { files } };
        const response = yield api.request("POST", `/${API_VERSION}/projects/${projectId}/rulesets`, {
            auth: true,
            data: payload,
            origin: api.rulesOrigin,
        });
        if (response.status === 200) {
            logger.debug("[rules] created ruleset", response.body.name);
            return response.body.name;
        }
        return _handleErrorResponse(response);
    });
}
exports.createRuleset = createRuleset;
function createRelease(projectId, rulesetName, releaseName) {
    return __awaiter(this, void 0, void 0, function* () {
        const payload = {
            name: `projects/${projectId}/releases/${releaseName}`,
            rulesetName,
        };
        const response = yield api.request("POST", `/${API_VERSION}/projects/${projectId}/releases`, {
            auth: true,
            data: payload,
            origin: api.rulesOrigin,
        });
        if (response.status === 200) {
            logger.debug("[rules] created release", response.body.name);
            return response.body.name;
        }
        return _handleErrorResponse(response);
    });
}
exports.createRelease = createRelease;
function updateRelease(projectId, rulesetName, releaseName) {
    return __awaiter(this, void 0, void 0, function* () {
        const payload = {
            release: {
                name: `projects/${projectId}/releases/${releaseName}`,
                rulesetName,
            },
        };
        const response = yield api.request("PATCH", `/${API_VERSION}/projects/${projectId}/releases/${releaseName}`, {
            auth: true,
            data: payload,
            origin: api.rulesOrigin,
        });
        if (response.status === 200) {
            logger.debug("[rules] updated release", response.body.name);
            return response.body.name;
        }
        return _handleErrorResponse(response);
    });
}
exports.updateRelease = updateRelease;
function updateOrCreateRelease(projectId, rulesetName, releaseName) {
    return __awaiter(this, void 0, void 0, function* () {
        logger.debug("[rules] releasing", releaseName, "with ruleset", rulesetName);
        return updateRelease(projectId, rulesetName, releaseName).catch(() => {
            logger.debug("[rules] ruleset update failed, attempting to create instead");
            return createRelease(projectId, rulesetName, releaseName);
        });
    });
}
exports.updateOrCreateRelease = updateOrCreateRelease;
function testRuleset(projectId, files) {
    return api.request("POST", `/${API_VERSION}/projects/${encodeURIComponent(projectId)}:test`, {
        origin: api.rulesOrigin,
        data: {
            source: { files },
        },
        auth: true,
    });
}
exports.testRuleset = testRuleset;
