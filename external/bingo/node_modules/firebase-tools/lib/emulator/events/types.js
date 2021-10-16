"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const _ = require("lodash");
class EventUtils {
    static isEvent(proto) {
        return _.has(proto, "context") && _.has(proto, "data");
    }
    static isLegacyEvent(proto) {
        return _.has(proto, "data") && _.has(proto, "resource");
    }
}
exports.EventUtils = EventUtils;
