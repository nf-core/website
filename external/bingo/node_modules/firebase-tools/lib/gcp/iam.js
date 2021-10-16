"use strict";
var api = require("../api");
var utils = require("../utils");
var API_VERSION = "v1";
function createServiceAccount(projectId, accountId, description, displayName) {
    return api
        .request("POST", `/${API_VERSION}/projects/${projectId}/serviceAccounts`, {
        auth: true,
        origin: api.iamOrigin,
        data: {
            accountId,
            serviceAccount: {
                displayName,
                description,
            },
        },
    })
        .then((res) => {
        return res.body;
    });
}
function deleteServiceAccount(projectId, accountEmail) {
    return api.request("DELETE", `/${API_VERSION}/projects/${projectId}/serviceAccounts/${accountEmail}`, {
        auth: true,
        origin: api.iamOrigin,
        resolveOnHTTPError: true,
    });
}
function getRole(role) {
    return api
        .request("GET", utils.endpoint([API_VERSION, "roles", role]), {
        auth: true,
        origin: api.iamOrigin,
        retryCodes: [500, 503],
    })
        .then(function (response) {
        return response.body;
    });
}
module.exports = {
    createServiceAccount,
    deleteServiceAccount,
    getRole,
};
