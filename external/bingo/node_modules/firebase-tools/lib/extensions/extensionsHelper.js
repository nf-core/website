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
const ora = require("ora");
const api_1 = require("../api");
const archiveDirectory_1 = require("../archiveDirectory");
const utils_1 = require("./utils");
const functionsConfig_1 = require("../functionsConfig");
const resolveSource_1 = require("./resolveSource");
const error_1 = require("../error");
const askUserForParam_1 = require("./askUserForParam");
const ensureApiEnabled_1 = require("../ensureApiEnabled");
const storage_1 = require("../gcp/storage");
const getProjectId = require("../getProjectId");
const extensionsApi_1 = require("./extensionsApi");
const prompt_1 = require("../prompt");
const logger = require("../logger");
const utils_2 = require("../utils");
var SpecParamType;
(function (SpecParamType) {
    SpecParamType["SELECT"] = "select";
    SpecParamType["MULTISELECT"] = "multiselect";
    SpecParamType["STRING"] = "string";
})(SpecParamType = exports.SpecParamType || (exports.SpecParamType = {}));
exports.logPrefix = "extensions";
const urlRegex = /^http[s]?:\/\/.*\.zip$/;
exports.EXTENSIONS_BUCKET_NAME = utils_2.envOverride("FIREBASE_EXTENSIONS_UPLOAD_BUCKET", "firebase-ext-eap-uploads");
exports.resourceTypeToNiceName = {
    "firebaseextensions.v1beta.scheduledFunction": "Scheduled Function",
    "firebaseextensions.v1beta.function": "Cloud Function",
};
function getDBInstanceFromURL(databaseUrl = "") {
    const instanceRegex = new RegExp("(?:https://)(.*)(?:.firebaseio.com)");
    const matches = databaseUrl.match(instanceRegex);
    if (matches && matches.length > 1) {
        return matches[1];
    }
    return "";
}
exports.getDBInstanceFromURL = getDBInstanceFromURL;
function getFirebaseProjectParams(projectId) {
    return __awaiter(this, void 0, void 0, function* () {
        const body = yield functionsConfig_1.getFirebaseConfig({ project: projectId });
        const FIREBASE_CONFIG = JSON.stringify({
            projectId: body.projectId,
            databaseURL: body.databaseURL,
            storageBucket: body.storageBucket,
        });
        return {
            PROJECT_ID: body.projectId,
            DATABASE_URL: body.databaseURL,
            STORAGE_BUCKET: body.storageBucket,
            FIREBASE_CONFIG,
            DATABASE_INSTANCE: getDBInstanceFromURL(body.databaseURL),
        };
    });
}
exports.getFirebaseProjectParams = getFirebaseProjectParams;
function substituteParams(original, params) {
    const startingString = JSON.stringify(original);
    const reduceFunction = (intermediateResult, paramVal, paramKey) => {
        const regex = new RegExp("\\$\\{" + paramKey + "\\}", "g");
        return intermediateResult.replace(regex, paramVal);
    };
    return JSON.parse(_.reduce(params, reduceFunction, startingString));
}
exports.substituteParams = substituteParams;
function populateDefaultParams(paramVars, paramSpec) {
    const newParams = paramVars;
    _.forEach(paramSpec, (env) => {
        if (!paramVars[env.param]) {
            if (env.default) {
                newParams[env.param] = env.default;
            }
            else {
                throw new error_1.FirebaseError(`${env.param} has not been set in the given params file` +
                    " and there is no default available. Please set this variable before installing again.");
            }
        }
    });
    return newParams;
}
exports.populateDefaultParams = populateDefaultParams;
function validateCommandLineParams(envVars, paramSpec) {
    if (_.size(envVars) > _.size(paramSpec)) {
        const paramList = _.map(paramSpec, (param) => {
            return param.param;
        });
        const misnamedParams = Object.keys(envVars).filter((key) => {
            return paramList.indexOf(key) === -1;
        });
        logger.info("Warning: The following params were specified in your env file but do not exist in the extension spec: " +
            `${misnamedParams.join(", ")}.`);
    }
    let allParamsValid = true;
    _.forEach(paramSpec, (param) => {
        if (!askUserForParam_1.checkResponse(envVars[param.param], param)) {
            allParamsValid = false;
        }
    });
    if (!allParamsValid) {
        throw new error_1.FirebaseError(`Some param values are not valid. Please check your params file.`);
    }
}
exports.validateCommandLineParams = validateCommandLineParams;
function validateSpec(spec) {
    const errors = [];
    if (!spec.name) {
        errors.push("extension.yaml is missing required field: name");
    }
    if (!spec.specVersion) {
        errors.push("extension.yaml is missing required field: specVersion");
    }
    if (!spec.version) {
        errors.push("extension.yaml; is missing required field: version");
    }
    for (let resource of spec.resources) {
        if (!resource.name) {
            errors.push("Resource is missing required field: name");
        }
        if (!resource.type) {
            errors.push(`Resource${resource.name ? ` ${resource.name}` : ""} is missing required field: type`);
        }
    }
    for (let api of spec.apis || []) {
        if (!api.apiName) {
            errors.push("API is missing required field: apiName");
        }
    }
    for (let role of spec.roles || []) {
        if (!role.role) {
            errors.push("Role is missing required field: role");
        }
    }
    for (let param of spec.params || []) {
        if (!param.param) {
            errors.push("Param is missing required field: param");
        }
        if (!param.label) {
            errors.push(`Param${param.param ? ` ${param.param}` : ""} is missing required field: label`);
        }
        if (param.type && !_.includes(SpecParamType, param.type)) {
            errors.push(`Invalid type ${param.type} for param${param.param ? ` ${param.param}` : ""}. Valid types are ${_.values(extensionsApi_1.ParamType).join(", ")}`);
        }
        if (!param.type || param.type == SpecParamType.STRING) {
            if (param.options) {
                errors.push(`Param${param.param ? ` ${param.param}` : ""} cannot have options because it is type STRING`);
            }
            if (param.default &&
                param.validationRegex &&
                !RegExp(param.validationRegex).test(param.default)) {
                errors.push(`Param${param.param ? ` ${param.param}` : ""} has default value '${param.default}', which does not pass the validationRegex ${param.validationRegex}`);
            }
        }
        if (param.type &&
            (param.type == SpecParamType.SELECT || param.type == SpecParamType.MULTISELECT)) {
            if (param.validationRegex) {
                errors.push(`Param${param.param ? ` ${param.param}` : ""} cannot have validationRegex because it is type ${param.type}`);
            }
            if (!param.options) {
                errors.push(`Param${param.param ? ` ${param.param}` : ""} requires options because it is type ${param.type}`);
            }
            for (let opt of param.options || []) {
                if (opt.value == undefined) {
                    errors.push(`Option for param${param.param ? ` ${param.param}` : ""} is missing required field: value`);
                }
            }
        }
    }
    if (errors.length) {
        const message = `The extension.yaml has the following errors: \n${errors.join("\n")}`;
        throw new error_1.FirebaseError(message);
    }
}
exports.validateSpec = validateSpec;
function promptForValidInstanceId(instanceId) {
    return __awaiter(this, void 0, void 0, function* () {
        let instanceIdIsValid = false;
        let newInstanceId;
        const instanceIdRegex = /^[a-z][a-z\d\-]*[a-z\d]$/;
        while (!instanceIdIsValid) {
            newInstanceId = yield prompt_1.promptOnce({
                type: "input",
                default: instanceId,
                message: `Please enter a new name for this instance:`,
            });
            if (newInstanceId.length <= 6 || 45 <= newInstanceId.length) {
                logger.info("Invalid instance ID. Instance ID must be between 6 and 45 characters.");
            }
            else if (!instanceIdRegex.test(newInstanceId)) {
                logger.info("Invalid instance ID. Instance ID must start with a lowercase letter, " +
                    "end with a lowercase letter or number, and only contain lowercase letters, numbers, or -");
            }
            else {
                instanceIdIsValid = true;
            }
        }
        return newInstanceId;
    });
}
exports.promptForValidInstanceId = promptForValidInstanceId;
function ensureExtensionsApiEnabled(options) {
    return __awaiter(this, void 0, void 0, function* () {
        const projectId = getProjectId(options);
        return yield ensureApiEnabled_1.ensure(projectId, "firebaseextensions.googleapis.com", "extensions", options.markdown);
    });
}
exports.ensureExtensionsApiEnabled = ensureExtensionsApiEnabled;
function archiveAndUploadSource(extPath, bucketName) {
    return __awaiter(this, void 0, void 0, function* () {
        const zippedSource = yield archiveDirectory_1.archiveDirectory(extPath, { type: "zip", ignore: ["node_modules"] });
        return yield storage_1.uploadObject(zippedSource, bucketName);
    });
}
function createSourceFromLocation(projectId, sourceUri) {
    return __awaiter(this, void 0, void 0, function* () {
        let packageUri;
        let extensionRoot;
        let objectPath = "";
        if (!urlRegex.test(sourceUri)) {
            const uploadSpinner = ora.default(" Archiving and uploading extension source code");
            try {
                uploadSpinner.start();
                objectPath = yield archiveAndUploadSource(sourceUri, exports.EXTENSIONS_BUCKET_NAME);
                uploadSpinner.succeed(" Uploaded extension source code");
                packageUri = api_1.firebaseStorageOrigin + objectPath + "?alt=media";
                extensionRoot = "/";
            }
            catch (err) {
                uploadSpinner.fail();
                throw err;
            }
        }
        else {
            [packageUri, extensionRoot] = sourceUri.split("#");
        }
        const res = yield extensionsApi_1.createSource(projectId, packageUri, extensionRoot);
        if (objectPath.length) {
            try {
                yield storage_1.deleteObject(objectPath);
                logger.debug("Cleaned up uploaded source archive");
            }
            catch (err) {
                logger.debug("Unable to clean up uploaded source archive");
            }
        }
        return res;
    });
}
exports.createSourceFromLocation = createSourceFromLocation;
function getExtensionSourceFromName(extensionName) {
    return __awaiter(this, void 0, void 0, function* () {
        const officialExtensionRegex = /^[a-zA-Z\-]+[0-9@.]*$/;
        const existingSourceRegex = /projects\/.+\/sources\/.+/;
        if (officialExtensionRegex.test(extensionName)) {
            const [name, version] = extensionName.split("@");
            const registryEntry = yield resolveSource_1.resolveRegistryEntry(name);
            const sourceUrl = yield resolveSource_1.resolveSourceUrl(registryEntry, name, version);
            return yield extensionsApi_1.getSource(sourceUrl);
        }
        else if (existingSourceRegex.test(extensionName)) {
            logger.info(`Fetching the source "${extensionName}"...`);
            return yield extensionsApi_1.getSource(extensionName);
        }
        throw new error_1.FirebaseError(`Could not find an extension named '${extensionName}'. `);
    });
}
exports.getExtensionSourceFromName = getExtensionSourceFromName;
function promptForOfficialExtension(message) {
    return __awaiter(this, void 0, void 0, function* () {
        const officialExts = yield resolveSource_1.getExtensionRegistry();
        return yield prompt_1.promptOnce({
            name: "input",
            type: "list",
            message,
            choices: utils_1.convertOfficialExtensionsToList(officialExts),
            pageSize: _.size(officialExts),
        });
    });
}
exports.promptForOfficialExtension = promptForOfficialExtension;
function promptForRepeatInstance(projectName, extensionName) {
    return __awaiter(this, void 0, void 0, function* () {
        const message = `An extension with the ID ${extensionName} already exists in the project ${projectName}.\n` +
            `Do you want to proceed with installing another instance of ${extensionName} in this project?`;
        return yield prompt_1.promptOnce({
            type: "confirm",
            message,
        });
    });
}
exports.promptForRepeatInstance = promptForRepeatInstance;
function instanceIdExists(projectId, instanceId) {
    return __awaiter(this, void 0, void 0, function* () {
        const instanceRes = yield extensionsApi_1.getInstance(projectId, instanceId, {
            resolveOnHTTPError: true,
        });
        if (instanceRes.error) {
            if (_.get(instanceRes, "error.code") === 404) {
                return false;
            }
            const msg = "Unexpected error when checking if instance ID exists: " +
                _.get(instanceRes, "error.message");
            throw new error_1.FirebaseError(msg, {
                original: instanceRes.error,
            });
        }
        return true;
    });
}
exports.instanceIdExists = instanceIdExists;
