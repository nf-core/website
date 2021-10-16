"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const fs_1 = require("fs");
function fileExistsSync(path) {
    try {
        return fs_1.statSync(path).isFile();
    }
    catch (e) {
        return false;
    }
}
exports.fileExistsSync = fileExistsSync;
function dirExistsSync(path) {
    try {
        return fs_1.statSync(path).isDirectory();
    }
    catch (e) {
        return false;
    }
}
exports.dirExistsSync = dirExistsSync;
