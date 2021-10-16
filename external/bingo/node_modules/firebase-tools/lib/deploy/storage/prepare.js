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
const gcp = require("../../gcp");
const rulesDeploy_1 = require("../../rulesDeploy");
function default_1(context, options) {
    return __awaiter(this, void 0, void 0, function* () {
        let rulesConfig = options.config.get("storage");
        if (!rulesConfig) {
            return;
        }
        _.set(context, "storage.rules", rulesConfig);
        const rulesDeploy = new rulesDeploy_1.RulesDeploy(options, rulesDeploy_1.RulesetServiceType.FIREBASE_STORAGE);
        _.set(context, "storage.rulesDeploy", rulesDeploy);
        if (_.isPlainObject(rulesConfig)) {
            const defaultBucket = yield gcp.storage.getDefaultBucket(options.project);
            rulesConfig = [_.assign(rulesConfig, { bucket: defaultBucket })];
            _.set(context, "storage.rules", rulesConfig);
        }
        rulesConfig.forEach((ruleConfig) => {
            if (ruleConfig.target) {
                options.rc.requireTarget(context.projectId, "storage", ruleConfig.target);
            }
            rulesDeploy.addFile(ruleConfig.rules);
        });
        yield rulesDeploy.compile();
    });
}
exports.default = default_1;
