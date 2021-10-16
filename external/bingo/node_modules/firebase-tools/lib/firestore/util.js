"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const error_1 = require("../error");
const INDEX_NAME_REGEX = /projects\/([^\/]+?)\/databases\/\(default\)\/collectionGroups\/([^\/]+?)\/indexes\/([^\/]*)/;
const FIELD_NAME_REGEX = /projects\/([^\/]+?)\/databases\/\(default\)\/collectionGroups\/([^\/]+?)\/fields\/([^\/]*)/;
function parseIndexName(name) {
    if (!name) {
        throw new error_1.FirebaseError(`Cannot parse undefined index name.`);
    }
    const m = name.match(INDEX_NAME_REGEX);
    if (!m || m.length < 4) {
        throw new error_1.FirebaseError(`Error parsing index name: ${name}`);
    }
    return {
        projectId: m[1],
        collectionGroupId: m[2],
        indexId: m[3],
    };
}
exports.parseIndexName = parseIndexName;
function parseFieldName(name) {
    const m = name.match(FIELD_NAME_REGEX);
    if (!m || m.length < 4) {
        throw new error_1.FirebaseError(`Error parsing field name: ${name}`);
    }
    return {
        projectId: m[1],
        collectionGroupId: m[2],
        fieldPath: m[3],
    };
}
exports.parseFieldName = parseFieldName;
