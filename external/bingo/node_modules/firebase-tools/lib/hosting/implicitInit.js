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
const clc = require("cli-color");
const fs = require("fs");
const fetchWebSetup_1 = require("../fetchWebSetup");
const utils = require("../utils");
const logger = require("../logger");
const INIT_TEMPLATE = fs.readFileSync(__dirname + "/../../templates/hosting/init.js", "utf8");
function implicitInit(options) {
    return __awaiter(this, void 0, void 0, function* () {
        let config;
        try {
            config = yield fetchWebSetup_1.fetchWebSetup(options);
        }
        catch (e) {
            logger.debug("fetchWebSetup error: " + e);
            const statusCode = _.get(e, "context.response.statusCode");
            if (statusCode === 403) {
                utils.logLabeledWarning("hosting", `Authentication error when trying to fetch your current web app configuration, have you run ${clc.bold("firebase login")}?`);
            }
        }
        if (!config) {
            config = fetchWebSetup_1.getCachedWebSetup(options);
            if (config) {
                utils.logLabeledWarning("hosting", "Using web app configuration from cache.");
            }
        }
        if (!config) {
            config = undefined;
            utils.logLabeledWarning("hosting", "Could not fetch web app configuration and there is no cached configuration on this machine. " +
                "Check your internet connection and make sure you are authenticated. " +
                "To continue, you must call firebase.initializeApp({...}) in your code before using Firebase.");
        }
        const configJson = JSON.stringify(config, null, 2);
        return {
            js: INIT_TEMPLATE.replace("/*--CONFIG--*/", `var firebaseConfig = ${configJson};`),
            json: configJson,
        };
    });
}
exports.implicitInit = implicitInit;
