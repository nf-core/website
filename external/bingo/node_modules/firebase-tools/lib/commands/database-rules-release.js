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
const command_1 = require("../command");
const requireInstance = require("../requireInstance");
const requirePermissions_1 = require("../requirePermissions");
const metadata = require("../database/metadata");
exports.default = new command_1.Command("database:rules:release <rulesetId>")
    .description("mark a staged ruleset as the stable ruleset")
    .option("--instance <instance>", "use the database <instance>.firebaseio.com (if omitted, uses default database instance)")
    .before(requirePermissions_1.requirePermissions, ["firebasedatabase.instances.update"])
    .before(requireInstance)
    .action((rulesetId, options) => __awaiter(void 0, void 0, void 0, function* () {
    const oldLabels = yield metadata.getRulesetLabels(options.instance);
    const newLabels = {
        stable: rulesetId,
        canary: oldLabels.canary,
    };
    yield metadata.setRulesetLabels(options.instance, newLabels);
    return newLabels;
}));
