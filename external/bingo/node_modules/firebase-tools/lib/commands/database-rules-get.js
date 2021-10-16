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
exports.default = new command_1.Command("database:rules:get <rulesetId>")
    .description("get a realtime database ruleset by id")
    .option("--instance <instance>", "use the database <instance>.firebaseio.com (if omitted, uses default database instance)")
    .before(requirePermissions_1.requirePermissions, ["firebasedatabase.instances.get"])
    .before(requireInstance)
    .action((rulesetId, options) => __awaiter(void 0, void 0, void 0, function* () {
    const ruleset = yield metadata.getRuleset(options.instance, rulesetId);
    logger.info(`Ruleset ${ruleset.id} was created at ${ruleset.createdAt}`);
    logger.info(ruleset.source);
    return ruleset;
}));
