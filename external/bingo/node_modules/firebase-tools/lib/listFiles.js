"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const glob_1 = require("glob");
function listFiles(cwd, ignore = []) {
    return glob_1.sync("**/*", {
        cwd,
        dot: true,
        follow: true,
        ignore: ["**/firebase-debug.log", ".firebase/*"].concat(ignore),
        nodir: true,
        nosort: true,
    });
}
exports.listFiles = listFiles;
