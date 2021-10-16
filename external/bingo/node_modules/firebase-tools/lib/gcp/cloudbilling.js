"use strict";
var api = require("../api");
var utils = require("../utils");
var API_VERSION = "v1";
function _checkBillingEnabled(projectId) {
    return api
        .request("GET", utils.endpoint([API_VERSION, "projects", projectId, "billingInfo"]), {
        auth: true,
        origin: api.cloudbillingOrigin,
        retryCodes: [500, 503],
    })
        .then(function (response) {
        return response.body.billingEnabled;
    });
}
function _setBillingAccount(projectId, billingAccount) {
    return api
        .request("PUT", utils.endpoint([API_VERSION, "projects", projectId, "billingInfo"]), {
        auth: true,
        origin: api.cloudbillingOrigin,
        retryCodes: [500, 503],
        data: {
            billingAccountName: billingAccount,
        },
    })
        .then(function (response) {
        return response.body.billingEnabled;
    });
}
function _listBillingAccounts() {
    return api
        .request("GET", utils.endpoint([API_VERSION, "billingAccounts"]), {
        auth: true,
        origin: api.cloudbillingOrigin,
        retryCodes: [500, 503],
    })
        .then(function (response) {
        return response.body.billingAccounts || [];
    });
}
module.exports = {
    checkBillingEnabled: _checkBillingEnabled,
    listBillingAccounts: _listBillingAccounts,
    setBillingAccount: _setBillingAccount,
};
