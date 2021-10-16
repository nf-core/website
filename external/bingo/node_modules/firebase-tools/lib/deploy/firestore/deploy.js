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
const error_1 = require("../../error");
const indexes_1 = require("../../firestore/indexes");
const logger = require("../../logger");
const utils = require("../../utils");
const rulesDeploy_1 = require("../../rulesDeploy");
function deployRules(context) {
    return __awaiter(this, void 0, void 0, function* () {
        const rulesDeploy = _.get(context, "firestore.rulesDeploy");
        if (!context.firestoreRules || !rulesDeploy) {
            return;
        }
        yield rulesDeploy.createRulesets(rulesDeploy_1.RulesetServiceType.CLOUD_FIRESTORE);
    });
}
function deployIndexes(context, options) {
    return __awaiter(this, void 0, void 0, function* () {
        if (!context.firestoreIndexes) {
            return;
        }
        const indexesFileName = _.get(context, "firestore.indexes.name");
        const indexesSrc = _.get(context, "firestore.indexes.content");
        if (!indexesSrc) {
            logger.debug("No Firestore indexes present.");
            return;
        }
        const indexes = indexesSrc.indexes;
        if (!indexes) {
            throw new error_1.FirebaseError(`Index file must contain an "indexes" property.`);
        }
        const fieldOverrides = indexesSrc.fieldOverrides || [];
        yield new indexes_1.FirestoreIndexes().deploy(options.project, indexes, fieldOverrides);
        utils.logSuccess(`${clc.bold.green("firestore:")} deployed indexes in ${clc.bold(indexesFileName)} successfully`);
    });
}
function default_1(context, options) {
    return __awaiter(this, void 0, void 0, function* () {
        yield Promise.all([deployRules(context), deployIndexes(context, options)]);
    });
}
exports.default = default_1;
