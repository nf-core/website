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
const clc = require("cli-color");
const fs = require("fs");
const logger = require("../../logger");
const prompt_1 = require("../../prompt");
const ensureCloudResourceLocation_1 = require("../../ensureCloudResourceLocation");
const RULES_TEMPLATE = fs.readFileSync(__dirname + "/../../../templates/init/storage/storage.rules", "utf8");
function doSetup(setup, config) {
    return __awaiter(this, void 0, void 0, function* () {
        setup.config.storage = {};
        ensureCloudResourceLocation_1.ensureLocationSet(setup.projectLocation, "Cloud Storage");
        logger.info();
        logger.info("Firebase Storage Security Rules allow you to define how and when to allow");
        logger.info("uploads and downloads. You can keep these rules in your project directory");
        logger.info("and publish them with " + clc.bold("firebase deploy") + ".");
        logger.info();
        const storageRulesFile = yield prompt_1.promptOnce({
            type: "input",
            name: "rules",
            message: "What file should be used for Storage Rules?",
            default: "storage.rules",
        });
        setup.config.storage.rules = storageRulesFile;
        config.writeProjectFile(setup.config.storage.rules, RULES_TEMPLATE);
    });
}
exports.doSetup = doSetup;
