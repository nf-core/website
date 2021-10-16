"use strict";
var fs = require("fs-extra");
var path = require("path");
var api = require("./api");
var { configstore } = require("./configstore");
var logger = require("./logger");
var configDir = function () {
    if (process.platform === "win32") {
        return process.env.APPDATA;
    }
    return process.env.HOME && path.resolve(process.env.HOME, ".config");
};
module.exports = function () {
    if (!configDir()) {
        logger.debug("Cannot ensure default credentials, no home directory found.");
        return;
    }
    var GCLOUD_CREDENTIAL_DIR = path.resolve(configDir(), "gcloud");
    var GCLOUD_CREDENTIAL_PATH = path.join(GCLOUD_CREDENTIAL_DIR, "application_default_credentials.json");
    var tokens = configstore.get("tokens") || {};
    var credentials = {
        client_id: api.clientId,
        client_secret: api.clientSecret,
        type: "authorized_user",
        refresh_token: tokens.refresh_token || process.env.FIREBASE_TOKEN,
    };
    fs.ensureDirSync(GCLOUD_CREDENTIAL_DIR);
    fs.writeFileSync(GCLOUD_CREDENTIAL_PATH, JSON.stringify(credentials, null, 2));
    return;
};
