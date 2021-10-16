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
const logger = require("../logger");
const requireInstance = require("../requireInstance");
const requirePermissions_1 = require("../requirePermissions");
const metadata = require("../database/metadata");
exports.default = new command_1.Command("database:rules:list")
    .description("list realtime database rulesets")
    .option("--instance <instance>", "use the database <instance>.firebaseio.com (if omitted, uses default database instance)")
    .before(requirePermissions_1.requirePermissions, ["firebasedatabase.instances.get"])
    .before(requireInstance)
    .action((options) => __awaiter(void 0, void 0, void 0, function* () {
    const labeled = yield metadata.getRulesetLabels(options.instance);
    const rulesets = yield metadata.listAllRulesets(options.instance);
    for (const ruleset of rulesets) {
        const labels = [];
        if (ruleset.id == labeled.stable) {
            labels.push("stable");
        }
        if (ruleset.id == labeled.canary) {
            labels.push("canary");
        }
        logger.info(`${ruleset.id}  ${ruleset.createdAt}  ${labels.join(",")}`);
    }
    logger.info("Labels:");
    logger.info(`  stable: ${labeled.stable}`);
    if (labeled.canary) {
        logger.info(`  canary: ${labeled.canary}`);
    }
    return { rulesets, labeled };
}));
