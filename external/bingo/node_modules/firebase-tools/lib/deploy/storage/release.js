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
function default_1(context, options) {
    return __awaiter(this, void 0, void 0, function* () {
        const rules = lodash_1.get(context, "storage.rules", []);
        const rulesDeploy = lodash_1.get(context, "storage.rulesDeploy");
        if (!rules.length || !rulesDeploy) {
            return;
        }
        const toRelease = [];
        for (const ruleConfig of rules) {
            if (ruleConfig.target) {
                options.rc
                    .target(options.project, "storage", ruleConfig.target)
                    .forEach(function (bucket) {
                    toRelease.push({ bucket: bucket, rules: ruleConfig.rules });
                });
            }
            else {
                toRelease.push({ bucket: ruleConfig.bucket, rules: ruleConfig.rules });
            }
        }
        yield Promise.all(toRelease.map((release) => {
            return rulesDeploy.release(release.rules, rulesDeploy_1.RulesetServiceType.FIREBASE_STORAGE, release.bucket);
        }));
    });
}
exports.default = default_1;
