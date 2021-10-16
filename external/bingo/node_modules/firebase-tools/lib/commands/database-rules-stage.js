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
const fs = require("fs-extra");
const path = require("path");
exports.default = new command_1.Command("database:rules:stage")
    .description("create a new realtime database ruleset")
    .option("--instance <instance>", "use the database <instance>.firebaseio.com (if omitted, uses default database instance)")
    .before(requirePermissions_1.requirePermissions, ["firebasedatabase.instances.update"])
    .before(requireInstance)
    .action((options) => __awaiter(void 0, void 0, void 0, function* () {
    const filepath = options.config.data.database.rules;
    logger.info(`staging ruleset from ${filepath}`);
    const source = fs.readFileSync(path.resolve(filepath), "utf8");
    const rulesetId = yield metadata.createRuleset(options.instance, source);
    logger.info(`staged ruleset ${rulesetId}`);
    return rulesetId;
}));
