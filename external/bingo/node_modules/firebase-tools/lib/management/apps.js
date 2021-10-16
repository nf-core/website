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
const fs = require("fs");
const api = require("../api");
const error_1 = require("../error");
const logger = require("../logger");
const operation_poller_1 = require("../operation-poller");
const TIMEOUT_MILLIS = 30000;
const APP_LIST_PAGE_SIZE = 100;
const CREATE_APP_API_REQUEST_TIMEOUT_MILLIS = 15000;
const WEB_CONFIG_FILE_NAME = "google-config.js";
var AppPlatform;
(function (AppPlatform) {
    AppPlatform["PLATFORM_UNSPECIFIED"] = "PLATFORM_UNSPECIFIED";
    AppPlatform["IOS"] = "IOS";
    AppPlatform["ANDROID"] = "ANDROID";
    AppPlatform["WEB"] = "WEB";
    AppPlatform["ANY"] = "ANY";
})(AppPlatform = exports.AppPlatform || (exports.AppPlatform = {}));
function getAppPlatform(platform) {
    switch (platform.toUpperCase()) {
        case "IOS":
            return AppPlatform.IOS;
        case "ANDROID":
            return AppPlatform.ANDROID;
        case "WEB":
            return AppPlatform.WEB;
        case "":
            return AppPlatform.ANY;
        default:
            throw new error_1.FirebaseError("Unexpected platform. Only iOS, Android, and Web apps are supported");
    }
}
exports.getAppPlatform = getAppPlatform;
function createIosApp(projectId, options) {
    return __awaiter(this, void 0, void 0, function* () {
        try {
            const response = yield api.request("POST", `/v1beta1/projects/${projectId}/iosApps`, {
                auth: true,
                origin: api.firebaseApiOrigin,
                timeout: CREATE_APP_API_REQUEST_TIMEOUT_MILLIS,
                data: options,
            });
            const appData = yield operation_poller_1.pollOperation({
                pollerName: "Create iOS app Poller",
                apiOrigin: api.firebaseApiOrigin,
                apiVersion: "v1beta1",
                operationResourceName: response.body.name,
            });
            return appData;
        }
        catch (err) {
            logger.debug(err.message);
            throw new error_1.FirebaseError(`Failed to create iOS app for project ${projectId}. See firebase-debug.log for more info.`, { exit: 2, original: err });
        }
    });
}
exports.createIosApp = createIosApp;
function createAndroidApp(projectId, options) {
    return __awaiter(this, void 0, void 0, function* () {
        try {
            const response = yield api.request("POST", `/v1beta1/projects/${projectId}/androidApps`, {
                auth: true,
                origin: api.firebaseApiOrigin,
                timeout: CREATE_APP_API_REQUEST_TIMEOUT_MILLIS,
                data: options,
            });
            const appData = yield operation_poller_1.pollOperation({
                pollerName: "Create Android app Poller",
                apiOrigin: api.firebaseApiOrigin,
                apiVersion: "v1beta1",
                operationResourceName: response.body.name,
            });
            return appData;
        }
        catch (err) {
            logger.debug(err.message);
            throw new error_1.FirebaseError(`Failed to create Android app for project ${projectId}. See firebase-debug.log for more info.`, {
                exit: 2,
                original: err,
            });
        }
    });
}
exports.createAndroidApp = createAndroidApp;
function createWebApp(projectId, options) {
    return __awaiter(this, void 0, void 0, function* () {
        try {
            const response = yield api.request("POST", `/v1beta1/projects/${projectId}/webApps`, {
                auth: true,
                origin: api.firebaseApiOrigin,
                timeout: CREATE_APP_API_REQUEST_TIMEOUT_MILLIS,
                data: options,
            });
            const appData = yield operation_poller_1.pollOperation({
                pollerName: "Create Web app Poller",
                apiOrigin: api.firebaseApiOrigin,
                apiVersion: "v1beta1",
                operationResourceName: response.body.name,
            });
            return appData;
        }
        catch (err) {
            logger.debug(err.message);
            throw new error_1.FirebaseError(`Failed to create Web app for project ${projectId}. See firebase-debug.log for more info.`, { exit: 2, original: err });
        }
    });
}
exports.createWebApp = createWebApp;
function getListAppsResourceString(projectId, platform) {
    let resourceSuffix;
    switch (platform) {
        case AppPlatform.IOS:
            resourceSuffix = "/iosApps";
            break;
        case AppPlatform.ANDROID:
            resourceSuffix = "/androidApps";
            break;
        case AppPlatform.WEB:
            resourceSuffix = "/webApps";
            break;
        case AppPlatform.ANY:
            resourceSuffix = ":searchApps";
            break;
        default:
            throw new error_1.FirebaseError("Unexpected platform. Only support iOS, Android and Web apps");
    }
    return `/v1beta1/projects/${projectId}${resourceSuffix}`;
}
function listFirebaseApps(projectId, platform, pageSize = APP_LIST_PAGE_SIZE) {
    return __awaiter(this, void 0, void 0, function* () {
        const apps = [];
        try {
            let nextPageToken = "";
            do {
                const pageTokenQueryString = nextPageToken ? `&pageToken=${nextPageToken}` : "";
                const response = yield api.request("GET", getListAppsResourceString(projectId, platform) +
                    `?pageSize=${pageSize}${pageTokenQueryString}`, {
                    auth: true,
                    origin: api.firebaseApiOrigin,
                    timeout: TIMEOUT_MILLIS,
                });
                if (response.body.apps) {
                    const appsOnPage = response.body.apps.map((app) => (app.platform ? app : Object.assign(Object.assign({}, app), { platform })));
                    apps.push(...appsOnPage);
                }
                nextPageToken = response.body.nextPageToken;
            } while (nextPageToken);
            return apps;
        }
        catch (err) {
            logger.debug(err.message);
            throw new error_1.FirebaseError(`Failed to list Firebase ${platform === AppPlatform.ANY ? "" : platform + " "}` +
                "apps. See firebase-debug.log for more info.", {
                exit: 2,
                original: err,
            });
        }
    });
}
exports.listFirebaseApps = listFirebaseApps;
function getAppConfigResourceString(appId, platform) {
    let platformResource;
    switch (platform) {
        case AppPlatform.IOS:
            platformResource = "iosApps";
            break;
        case AppPlatform.ANDROID:
            platformResource = "androidApps";
            break;
        case AppPlatform.WEB:
            platformResource = "webApps";
            break;
        default:
            throw new error_1.FirebaseError("Unexpected app platform");
    }
    return `/v1beta1/projects/-/${platformResource}/${appId}/config`;
}
function parseConfigFromResponse(responseBody, platform) {
    if (platform === AppPlatform.WEB) {
        const JS_TEMPLATE = fs.readFileSync(__dirname + "/../../templates/setup/web.js", "utf8");
        return {
            fileName: WEB_CONFIG_FILE_NAME,
            fileContents: JS_TEMPLATE.replace("{/*--CONFIG--*/}", JSON.stringify(responseBody, null, 2)),
        };
    }
    else if (platform === AppPlatform.ANDROID || platform === AppPlatform.IOS) {
        return {
            fileName: responseBody.configFilename,
            fileContents: Buffer.from(responseBody.configFileContents, "base64").toString("utf8"),
        };
    }
    throw new error_1.FirebaseError("Unexpected app platform");
}
function getAppConfigFile(config, platform) {
    return parseConfigFromResponse(config, platform);
}
exports.getAppConfigFile = getAppConfigFile;
function getAppConfig(appId, platform) {
    return __awaiter(this, void 0, void 0, function* () {
        let response;
        try {
            response = yield api.request("GET", getAppConfigResourceString(appId, platform), {
                auth: true,
                origin: api.firebaseApiOrigin,
                timeout: TIMEOUT_MILLIS,
            });
        }
        catch (err) {
            logger.debug(err.message);
            throw new error_1.FirebaseError(`Failed to get ${platform} app configuration. See firebase-debug.log for more info.`, {
                exit: 2,
                original: err,
            });
        }
        return response.body;
    });
}
exports.getAppConfig = getAppConfig;
