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
const command_1 = require("../command");
const getProjectId = require("../getProjectId");
const extensionsHelper_1 = require("../extensions/extensionsHelper");
const listExtensions_1 = require("../extensions/listExtensions");
const requirePermissions_1 = require("../requirePermissions");
const logger = require("../logger");
const utils = require("../utils");
module.exports = new command_1.Command("ext")
    .description("display information on how to use ext commands and extensions installed to your project")
    .action((options) => __awaiter(void 0, void 0, void 0, function* () {
    utils.logLabeledBullet(extensionsHelper_1.logPrefix, "list of extensions commands:");
    const firebaseTools = require("../");
    const commandNames = [
        "ext:install",
        "ext:info",
        "ext:list",
        "ext:configure",
        "ext:update",
        "ext:uninstall",
    ];
    _.forEach(commandNames, (commandName) => {
        const command = firebaseTools.getCommand(commandName);
        logger.info(clc.bold("\n" + command.name()));
        command.outputHelp();
    });
    logger.info();
    try {
        yield requirePermissions_1.requirePermissions(options, ["firebasemods.instances.list"]);
        const projectId = getProjectId(options);
        return listExtensions_1.listExtensions(projectId);
    }
    catch (err) {
        return;
    }
}));
