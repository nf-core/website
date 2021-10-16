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
const gcp = require("../../../gcp");
const fsutils = require("../../../fsutils");
const prompt_1 = require("../../../prompt");
const logger = require("../../../logger");
const utils = require("../../../utils");
const DEFAULT_RULES_FILE = "firestore.rules";
const RULES_TEMPLATE = fs.readFileSync(__dirname + "/../../../../templates/init/firestore/firestore.rules", "utf8");
function initRules(setup, config) {
    return __awaiter(this, void 0, void 0, function* () {
        logger.info();
        logger.info("Firestore Security Rules allow you to define how and when to allow");
        logger.info("requests. You can keep these rules in your project directory");
        logger.info("and publish them with " + clc.bold("firebase deploy") + ".");
        logger.info();
        return prompt_1.prompt(setup.config.firestore, [
            {
                type: "input",
                name: "rules",
                message: "What file should be used for Firestore Rules?",
                default: DEFAULT_RULES_FILE,
            },
        ])
            .then(() => {
            const filename = setup.config.firestore.rules;
            if (fsutils.fileExistsSync(filename)) {
                const msg = "File " +
                    clc.bold(filename) +
                    " already exists." +
                    " Do you want to overwrite it with the Firestore Rules from the Firebase Console?";
                return prompt_1.promptOnce({
                    type: "confirm",
                    message: msg,
                    default: false,
                });
            }
            return Promise.resolve(true);
        })
            .then((overwrite) => {
            if (!overwrite) {
                return Promise.resolve();
            }
            return getRulesFromConsole(setup.projectId).then((contents) => {
                return config.writeProjectFile(setup.config.firestore.rules, contents);
            });
        });
    });
}
exports.initRules = initRules;
function getRulesFromConsole(projectId) {
    return __awaiter(this, void 0, void 0, function* () {
        return gcp.rules
            .getLatestRulesetName(projectId, "cloud.firestore")
            .then((name) => {
            if (!name) {
                logger.debug("No rulesets found, using default.");
                return [{ name: DEFAULT_RULES_FILE, content: RULES_TEMPLATE }];
            }
            logger.debug("Found ruleset: " + name);
            return gcp.rules.getRulesetContent(name);
        })
            .then((rules) => {
            if (rules.length <= 0) {
                return utils.reject("Ruleset has no files", { exit: 1 });
            }
            if (rules.length > 1) {
                return utils.reject("Ruleset has too many files: " + rules.length, { exit: 1 });
            }
            return rules[0].content;
        });
    });
}
