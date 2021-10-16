"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const clc = require("cli-color");
const error_1 = require("../error");
function assertHas(obj, prop) {
    const objString = clc.cyan(JSON.stringify(obj));
    if (!obj[prop]) {
        throw new error_1.FirebaseError(`Must contain "${prop}": ${objString}`);
    }
}
exports.assertHas = assertHas;
function assertHasOneOf(obj, props) {
    const objString = clc.cyan(JSON.stringify(obj));
    let count = 0;
    props.forEach((prop) => {
        if (obj[prop]) {
            count++;
        }
    });
    if (count !== 1) {
        throw new error_1.FirebaseError(`Must contain exactly one of "${props.join(",")}": ${objString}`);
    }
}
exports.assertHasOneOf = assertHasOneOf;
function assertEnum(obj, prop, valid) {
    const objString = clc.cyan(JSON.stringify(obj));
    if (valid.indexOf(obj[prop]) < 0) {
        throw new error_1.FirebaseError(`Field "${prop}" must be one of  ${valid.join(", ")}: ${objString}`);
    }
}
exports.assertEnum = assertEnum;
