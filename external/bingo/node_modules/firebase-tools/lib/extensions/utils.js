"use strict";
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
Object.defineProperty(exports, "__esModule", { value: true });
const _ = require("lodash");
const prompt_1 = require("../prompt");
function onceWithJoin(question) {
    return __awaiter(this, void 0, void 0, function* () {
        const response = yield prompt_1.promptOnce(question);
        if (Array.isArray(response)) {
            return response.join(",");
        }
        return response;
    });
}
exports.onceWithJoin = onceWithJoin;
function convertExtensionOptionToLabeledList(options) {
    return options.map((option) => {
        return {
            checked: false,
            name: option.label,
            value: option.value,
        };
    });
}
exports.convertExtensionOptionToLabeledList = convertExtensionOptionToLabeledList;
function convertOfficialExtensionsToList(officialExts) {
    return _.map(officialExts, (entry, key) => {
        return {
            checked: false,
            value: key,
        };
    });
}
exports.convertOfficialExtensionsToList = convertOfficialExtensionsToList;
function getRandomString(length) {
    const SUFFIX_CHAR_SET = "abcdefghijklmnopqrstuvwxyz0123456789";
    let result = "";
    for (let i = 0; i < length; i++) {
        result += SUFFIX_CHAR_SET.charAt(Math.floor(Math.random() * SUFFIX_CHAR_SET.length));
    }
    return result;
}
exports.getRandomString = getRandomString;
