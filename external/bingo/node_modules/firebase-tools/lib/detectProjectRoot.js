"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const fsutils_1 = require("./fsutils");
const path_1 = require("path");
function detectProjectRoot(cwd) {
    let projectRootDir = cwd || process.cwd();
    while (!fsutils_1.fileExistsSync(path_1.resolve(projectRootDir, "./firebase.json"))) {
        const parentDir = path_1.dirname(projectRootDir);
        if (parentDir === projectRootDir) {
            return null;
        }
        projectRootDir = parentDir;
    }
    return projectRootDir;
}
exports.detectProjectRoot = detectProjectRoot;
