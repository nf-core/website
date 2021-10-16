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
const fs = require("fs-extra");
const command_1 = require("../command");
const utils = require("../utils");
const requireAuth = require("../requireAuth");
const client_1 = require("../appdistribution/client");
const error_1 = require("../error");
const distribution_1 = require("../appdistribution/distribution");
function ensureFileExists(file, message = "") {
    if (!fs.existsSync(file)) {
        throw new error_1.FirebaseError(`File ${file} does not exist: ${message}`);
    }
}
function getAppId(appId) {
    if (!appId) {
        throw new error_1.FirebaseError("set the --app option to a valid Firebase app id and try again");
    }
    return appId;
}
function getReleaseNotes(releaseNotes, releaseNotesFile) {
    if (releaseNotes) {
        return releaseNotes.replace(/\\n/g, "\n");
    }
    else if (releaseNotesFile) {
        ensureFileExists(releaseNotesFile);
        return fs.readFileSync(releaseNotesFile, "utf8");
    }
    return "";
}
function getTestersOrGroups(value, file) {
    if (!value && file) {
        ensureFileExists(file);
        value = fs.readFileSync(file, "utf8");
    }
    if (value) {
        return value
            .split(/,|\n/)
            .map((entry) => entry.trim())
            .filter((entry) => !!entry);
    }
    return [];
}
module.exports = new command_1.Command("appdistribution:distribute <distribution-file>")
    .description("upload a distribution")
    .option("--app <app_id>", "the app id of your Firebase app")
    .option("--release-notes <string>", "release notes to include with this distribution")
    .option("--release-notes-file <file>", "path to file with release notes to include with this distribution")
    .option("--testers <string>", "a comma separated list of tester emails to distribute to")
    .option("--testers-file <file>", "path to file with a comma separated list of tester emails to distribute to")
    .option("--groups <string>", "a comma separated list of group aliases to distribute to")
    .option("--groups-file <file>", "path to file with a comma separated list of group aliases to distribute to")
    .before(requireAuth)
    .action((file, options) => __awaiter(void 0, void 0, void 0, function* () {
    const appId = getAppId(options.app);
    const distribution = new distribution_1.Distribution(file);
    const releaseNotes = getReleaseNotes(options.releaseNotes, options.releaseNotesFile);
    const testers = getTestersOrGroups(options.testers, options.testersFile);
    const groups = getTestersOrGroups(options.groups, options.groupsFile);
    const requests = new client_1.AppDistributionClient(appId);
    let app;
    try {
        app = yield requests.getApp();
    }
    catch (err) {
        if (err.status === 404) {
            throw new error_1.FirebaseError(`App Distribution could not find your app ${appId}. ` +
                `Make sure to onboard your app by pressing the "Get started" ` +
                "button on the App Distribution page in the Firebase console: " +
                "https://console.firebase.google.com/project/_/appdistribution", { exit: 1 });
        }
        throw new error_1.FirebaseError(`failed to fetch app information. ${err.message}`, { exit: 1 });
    }
    if (!app.contactEmail) {
        throw new error_1.FirebaseError(`We could not find a contact email for app ${appId}. Please visit App Distribution within ` +
            "the Firebase Console to set one up.", { exit: 1 });
    }
    const releaseHash = yield distribution.releaseHash();
    let releaseId;
    const uploadStatus = yield requests.getUploadStatus(releaseHash);
    if (uploadStatus.status === client_1.UploadStatus.SUCCESS) {
        utils.logWarning("this distribution has been uploaded before, skipping upload");
        releaseId = uploadStatus.release.id;
    }
    else {
        utils.logBullet("uploading distribution...");
        try {
            const token = yield requests.getJwtToken();
            const releaseEtag = yield requests.uploadDistribution(token, distribution);
            releaseId = yield requests.pollReleaseIdByHash(releaseEtag);
            utils.logSuccess("uploaded distribution successfully!");
        }
        catch (err) {
            throw new error_1.FirebaseError(`failed to upload distribution. ${err.message}`, { exit: 1 });
        }
    }
    yield requests.addReleaseNotes(releaseId, releaseNotes);
    yield requests.enableAccess(releaseId, testers, groups);
}));
