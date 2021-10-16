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
const yaml = require("js-yaml");
const path = require("path");
const fs = require("fs-extra");
const fsutils_1 = require("../../fsutils");
const error_1 = require("../../error");
const paramHelper_1 = require("./paramHelper");
const SPEC_FILE = "extension.yaml";
const validFunctionTypes = [
    "firebaseextensions.v1beta.function",
    "firebaseextensions.v1beta.scheduledFunction",
];
function wrappedSafeLoad(source) {
    try {
        return yaml.safeLoad(source);
    }
    catch (err) {
        if (err instanceof yaml.YAMLException) {
            throw new error_1.FirebaseError(`YAML Error: ${err.message}`, { original: err });
        }
        throw err;
    }
}
function findExtensionYaml(directory) {
    return __awaiter(this, void 0, void 0, function* () {
        while (!fsutils_1.fileExistsSync(path.resolve(directory, SPEC_FILE))) {
            const parentDir = path.dirname(directory);
            if (parentDir === directory) {
                throw new error_1.FirebaseError("Couldn't find an extension.yaml file. Check that you are in the root directory of your extension.");
            }
            directory = parentDir;
        }
        return directory;
    });
}
exports.findExtensionYaml = findExtensionYaml;
function readExtensionYaml(directory) {
    return __awaiter(this, void 0, void 0, function* () {
        const extensionYaml = yield readFileFromDirectory(directory, SPEC_FILE);
        const source = extensionYaml.source;
        return wrappedSafeLoad(source);
    });
}
exports.readExtensionYaml = readExtensionYaml;
function readFileFromDirectory(directory, file) {
    return __awaiter(this, void 0, void 0, function* () {
        return new Promise((resolve, reject) => {
            fs.readFile(path.resolve(directory, file), "utf8", (err, data) => {
                if (err) {
                    if (err.code === "ENOENT") {
                        return reject(new error_1.FirebaseError(`Could not find "${file}" in "${directory}"`, { original: err }));
                    }
                    reject(new error_1.FirebaseError(`Failed to read file "${file}" in "${directory}"`, { original: err }));
                }
                else {
                    resolve(data);
                }
            });
        }).then((source) => {
            return {
                source,
                sourceDirectory: directory,
            };
        });
    });
}
exports.readFileFromDirectory = readFileFromDirectory;
function getFunctionResourcesWithParamSubstitution(extensionSpec, params) {
    const rawResources = extensionSpec.resources.filter((resource) => validFunctionTypes.includes(resource.type));
    return paramHelper_1.substituteParams(rawResources, params);
}
exports.getFunctionResourcesWithParamSubstitution = getFunctionResourcesWithParamSubstitution;
function getFunctionProperties(resources) {
    return resources.map((r) => r.properties);
}
exports.getFunctionProperties = getFunctionProperties;
