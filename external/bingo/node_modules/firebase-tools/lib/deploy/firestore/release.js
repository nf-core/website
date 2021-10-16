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
const rulesDeploy_1 = require("../../rulesDeploy");
function default_1(context, options) {
    return __awaiter(this, void 0, void 0, function* () {
        const rulesDeploy = _.get(context, "firestore.rulesDeploy");
        if (!context.firestoreRules || !rulesDeploy) {
            return;
        }
        yield rulesDeploy.release(options.config.get("firestore.rules"), rulesDeploy_1.RulesetServiceType.CLOUD_FIRESTORE);
    });
}
exports.default = default_1;
