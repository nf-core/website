"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const _ = require("lodash");
const path = require("path");
const fs = require("fs-extra");
const dotenv = require("dotenv");
const error_1 = require("../../error");
function readParamsFile(envFilePath) {
    try {
        const buf = fs.readFileSync(path.resolve(envFilePath), "utf8");
        return dotenv.parse(buf.toString().trim(), { debug: true });
    }
    catch (err) {
        throw new error_1.FirebaseError(`Error reading --test-params file: ${err.message}\n`, {
            original: err,
        });
    }
}
exports.readParamsFile = readParamsFile;
function substituteParams(original, params) {
    const startingString = JSON.stringify(original);
    const reduceFunction = (intermediateResult, paramVal, paramKey) => {
        const regex = new RegExp("\\$\\{" + paramKey + "\\}", "g");
        return intermediateResult.replace(regex, paramVal);
    };
    return JSON.parse(_.reduce(params, reduceFunction, startingString));
}
exports.substituteParams = substituteParams;
