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
const loadCJSON = require("../../loadCJSON");
const rulesDeploy_1 = require("../../rulesDeploy");
const utils = require("../../utils");
function prepareRules(context, options) {
    return __awaiter(this, void 0, void 0, function* () {
        const rulesFile = options.config.get("firestore.rules");
        if (context.firestoreRules && rulesFile) {
            const rulesDeploy = new rulesDeploy_1.RulesDeploy(options, rulesDeploy_1.RulesetServiceType.CLOUD_FIRESTORE);
            _.set(context, "firestore.rulesDeploy", rulesDeploy);
            rulesDeploy.addFile(rulesFile);
            yield rulesDeploy.compile();
        }
    });
}
function prepareIndexes(context, options) {
    if (!context.firestoreIndexes || !options.config.get("firestore.indexes")) {
        return;
    }
    const indexesFileName = options.config.get("firestore.indexes");
    const indexesPath = options.config.path(indexesFileName);
    const parsedSrc = loadCJSON(indexesPath);
    utils.logBullet(`${clc.bold.cyan("firestore:")} reading indexes from ${clc.bold(indexesFileName)}...`);
    context.firestore = context.firestore || {};
    context.firestore.indexes = {
        name: indexesFileName,
        content: parsedSrc,
    };
}
function default_1(context, options) {
    return __awaiter(this, void 0, void 0, function* () {
        if (options.only) {
            const targets = options.only.split(",");
            const onlyIndexes = targets.indexOf("firestore:indexes") >= 0;
            const onlyRules = targets.indexOf("firestore:rules") >= 0;
            const onlyFirestore = targets.indexOf("firestore") >= 0;
            context.firestoreIndexes = onlyIndexes || onlyFirestore;
            context.firestoreRules = onlyRules || onlyFirestore;
        }
        else {
            context.firestoreIndexes = true;
            context.firestoreRules = true;
        }
        prepareIndexes(context, options);
        yield prepareRules(context, options);
    });
}
exports.default = default_1;
