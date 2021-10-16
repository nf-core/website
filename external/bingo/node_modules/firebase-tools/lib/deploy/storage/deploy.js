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
const lodash_1 = require("lodash");
const rulesDeploy_1 = require("../../rulesDeploy");
function default_1(context) {
    return __awaiter(this, void 0, void 0, function* () {
        const rulesDeploy = lodash_1.get(context, "storage.rulesDeploy");
        if (!rulesDeploy) {
            return;
        }
        yield rulesDeploy.createRulesets(rulesDeploy_1.RulesetServiceType.FIREBASE_STORAGE);
    });
}
exports.default = default_1;
