"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const _ = require("lodash");
const error_1 = require("../error");
function isPlainObject(input) {
    return (typeof input === "object" &&
        input !== null &&
        _.isEqual(Object.getPrototypeOf(input), Object.prototype));
}
function encodeHelper(val) {
    if (_.isString(val)) {
        return { stringValue: val };
    }
    if (_.isBoolean(val)) {
        return { booleanValue: val };
    }
    if (_.isInteger(val)) {
        return { integerValue: val };
    }
    if (_.isNumber(val)) {
        return { doubleValue: val };
    }
    if (_.isDate(val)) {
        return { timestampValue: val.toISOString() };
    }
    if (_.isArray(val)) {
        const encodedElements = [];
        for (const v of val) {
            const enc = encodeHelper(v);
            if (enc) {
                encodedElements.push(enc);
            }
        }
        return {
            arrayValue: { values: encodedElements },
        };
    }
    if (_.isNull(val)) {
        return { nullValue: "NULL_VALUE" };
    }
    if (val instanceof Buffer || val instanceof Uint8Array) {
        return { bytesValue: val };
    }
    if (isPlainObject(val)) {
        return {
            mapValue: { fields: encodeFirestoreValue(val) },
        };
    }
    throw new error_1.FirebaseError(`Cannot encode ${val} to a Firestore Value. ` +
        "The emulator does not yet support Firestore document reference values or geo points.");
}
function encodeFirestoreValue(data) {
    return _.mapValues(data, encodeHelper);
}
exports.encodeFirestoreValue = encodeFirestoreValue;
