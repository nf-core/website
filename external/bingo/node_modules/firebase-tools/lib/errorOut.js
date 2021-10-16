"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const logError = require("./logError");
const error_1 = require("./error");
function errorOut(error) {
    let fbError;
    if (error instanceof error_1.FirebaseError) {
        fbError = error;
    }
    else {
        fbError = new error_1.FirebaseError("An unexpected error has occurred.", {
            original: error,
            exit: 2,
        });
    }
    logError(fbError);
    process.exitCode = fbError.exit || 2;
    setTimeout(function () {
        process.exit();
    }, 250);
}
exports.errorOut = errorOut;
