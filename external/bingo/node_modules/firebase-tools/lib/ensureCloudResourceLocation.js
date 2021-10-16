"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const error_1 = require("./error");
function ensureLocationSet(location, feature) {
    if (!location) {
        throw new error_1.FirebaseError(`Cloud resource location is not set for this project but the operation ` +
            `you are attempting to perform in ${feature} requires it. ` +
            `Please see this documentation for more details: https://firebase.google.com/docs/projects/locations`);
    }
}
exports.ensureLocationSet = ensureLocationSet;
