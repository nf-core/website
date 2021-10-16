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
const requirePermissions_1 = require("../requirePermissions");
const getProjectNumber = require("../getProjectNumber");
const firedata = require("../gcp/firedata");
exports.default = new command_1.Command("database:instances:list")
    .description("list realtime database instances")
    .before(requirePermissions_1.requirePermissions, ["firebasedatabase.instances.list"])
    .action((options) => __awaiter(void 0, void 0, void 0, function* () {
    const projectNumber = yield getProjectNumber(options);
    const instances = yield firedata.listDatabaseInstances(projectNumber);
    for (const instance of instances) {
        logger.info(instance.instance);
    }
    logger.info(`Project ${options.project} has ${instances.length} database instances`);
    return instances;
}));
