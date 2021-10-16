"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const _ = require("lodash");
function populatePostinstall(instructions, params) {
    return _.reduce(params, (content, value, key) => {
        const regex = new RegExp("\\$\\{param:" + key + "\\}", "g");
        return _.replace(content, regex, value);
    }, instructions);
}
exports.populatePostinstall = populatePostinstall;
