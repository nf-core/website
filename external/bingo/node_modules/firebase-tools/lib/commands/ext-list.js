"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const command_1 = require("../command");
const getProjectId = require("../getProjectId");
const listExtensions_1 = require("../extensions/listExtensions");
const extensionsHelper_1 = require("../extensions/extensionsHelper");
const requirePermissions_1 = require("../requirePermissions");
module.exports = new command_1.Command("ext:list")
    .description("list all the extensions that are installed in your Firebase project")
    .before(requirePermissions_1.requirePermissions, ["firebasemods.instances.list"])
    .before(extensionsHelper_1.ensureExtensionsApiEnabled)
    .action((options) => {
    const projectId = getProjectId(options);
    return listExtensions_1.listExtensions(projectId);
});
