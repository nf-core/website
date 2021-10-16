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
const clc = require("cli-color");
const marked = require("marked");
const semver = require("semver");
const api = require("../api");
const error_1 = require("../error");
const updateHelper_1 = require("./updateHelper");
const logger = require("../logger");
const prompt_1 = require("../prompt");
const EXTENSIONS_REGISTRY_ENDPOINT = "/extensions.json";
const AUDIENCE_WARNING_MESSAGES = {
    "open-alpha": marked(`${clc.bold("Important")}: This extension is part of the ${clc.bold("preliminary-release program")} for extensions.\n Its functionality might change in backward-incompatible ways before its official release. Learn more: https://github.com/firebase/extensions/tree/master/.preliminary-release-extensions`),
    "closed-alpha": marked(`${clc.yellow.bold("Important")}: This extension is part of the ${clc.bold("Firebase Alpha program")}.\n This extension is strictly confidential, and its functionality might change in backward-incompatible ways before its official, public release. Learn more: https://dev-partners.googlesource.com/samples/firebase/extensions-alpha/+/refs/heads/master/README.md`),
};
function resolveSourceUrl(registryEntry, name, version) {
    const targetVersion = getTargetVersion(registryEntry, version);
    const sourceUrl = _.get(registryEntry, ["versions", targetVersion]);
    if (!sourceUrl) {
        throw new error_1.FirebaseError(`Could not find version ${clc.bold(version)} of extension ${clc.bold(name)}.`);
    }
    return sourceUrl;
}
exports.resolveSourceUrl = resolveSourceUrl;
function resolveRegistryEntry(name) {
    return __awaiter(this, void 0, void 0, function* () {
        const extensionsRegistry = yield getExtensionRegistry();
        const registryEntry = _.get(extensionsRegistry, name);
        if (!registryEntry) {
            throw new error_1.FirebaseError(`Unable to find extension source named ${clc.bold(name)}.`);
        }
        return registryEntry;
    });
}
exports.resolveRegistryEntry = resolveRegistryEntry;
function getTargetVersion(registryEntry, versionOrLabel) {
    const seekVersion = versionOrLabel || "latest";
    const versionFromLabel = _.get(registryEntry, ["labels", seekVersion]);
    return versionFromLabel || seekVersion;
}
exports.getTargetVersion = getTargetVersion;
function promptForUpdateWarnings(registryEntry, startVersion, endVersion) {
    return __awaiter(this, void 0, void 0, function* () {
        if (registryEntry.updateWarnings) {
            for (const targetRange in registryEntry.updateWarnings) {
                if (semver.satisfies(endVersion, targetRange)) {
                    const updateWarnings = registryEntry.updateWarnings[targetRange];
                    for (const updateWarning of updateWarnings) {
                        if (semver.satisfies(startVersion, updateWarning.from)) {
                            yield updateHelper_1.confirmUpdateWarning(updateWarning);
                            break;
                        }
                    }
                    break;
                }
            }
        }
    });
}
exports.promptForUpdateWarnings = promptForUpdateWarnings;
function promptForAudienceConsent(registryEntry) {
    return __awaiter(this, void 0, void 0, function* () {
        let consent = true;
        if (registryEntry.audience && AUDIENCE_WARNING_MESSAGES[registryEntry.audience]) {
            logger.info(AUDIENCE_WARNING_MESSAGES[registryEntry.audience]);
            consent = yield prompt_1.promptOnce({
                type: "confirm",
                message: "Do you acknowledge the status of this extension?",
                default: true,
            });
        }
        return consent;
    });
}
exports.promptForAudienceConsent = promptForAudienceConsent;
function getExtensionRegistry() {
    return __awaiter(this, void 0, void 0, function* () {
        const res = yield api.request("GET", EXTENSIONS_REGISTRY_ENDPOINT, {
            origin: api.firebaseExtensionsRegistryOrigin,
        });
        return res.body.mods;
    });
}
exports.getExtensionRegistry = getExtensionRegistry;
