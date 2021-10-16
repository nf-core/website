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
const path = require("path");
const yaml = require("js-yaml");
const error_1 = require("../error");
const logger = require("../logger");
const EXTENSIONS_SPEC_FILE = "extension.yaml";
const EXTENSIONS_PREINSTALL_FILE = "PREINSTALL.md";
function getLocalExtensionSpec(directory) {
    return __awaiter(this, void 0, void 0, function* () {
        const spec = yield parseYAML(yield readFile(path.resolve(directory, EXTENSIONS_SPEC_FILE)));
        try {
            const preinstall = yield readFile(path.resolve(directory, EXTENSIONS_PREINSTALL_FILE));
            spec.preinstallContent = preinstall;
        }
        catch (err) {
            logger.debug(`No PREINSTALL.md found in directory ${directory}.`);
        }
        return spec;
    });
}
exports.getLocalExtensionSpec = getLocalExtensionSpec;
function readFile(pathToFile) {
    return __awaiter(this, void 0, void 0, function* () {
        try {
            return fs.readFileSync(pathToFile, "utf8");
        }
        catch (err) {
            if (err.code === "ENOENT") {
                throw new error_1.FirebaseError(`Could not find "${pathToFile}""`, { original: err });
            }
            throw new error_1.FirebaseError(`Failed to read file at "${pathToFile}"`, { original: err });
        }
    });
}
exports.readFile = readFile;
function isLocalExtension(extensionName) {
    try {
        fs.readdirSync(extensionName);
    }
    catch (err) {
        return false;
    }
    return true;
}
exports.isLocalExtension = isLocalExtension;
function parseYAML(source) {
    try {
        return yaml.safeLoad(source);
    }
    catch (err) {
        if (err instanceof yaml.YAMLException) {
            throw new error_1.FirebaseError(`YAML Error: ${err.message}`, { original: err });
        }
        throw new error_1.FirebaseError(err.message);
    }
}
