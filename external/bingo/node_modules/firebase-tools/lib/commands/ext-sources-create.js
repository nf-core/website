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
const getProjectId = require("../getProjectId");
const logger = require("../logger");
const extensionsHelper_1 = require("../extensions/extensionsHelper");
const requirePermissions_1 = require("../requirePermissions");
exports.default = new command_1.Command("ext:sources:create <sourceLocation>")
    .description(`create a extension source from sourceLocation`)
    .help("sourceLocation can be a local directory containing an extension, or URL pointing to a zipped extension. " +
    'If using a URL, you can specify a root folder for the extension by adding "#<extensionRoot>". ' +
    "For example, if your extension.yaml is in the my/extension directory of the archive, " +
    "you should use sourceUrl#my/extension. If no extensionRoot is specified, / is assumed.")
    .before(requirePermissions_1.requirePermissions, ["firebasemods.sources.create"])
    .before(extensionsHelper_1.ensureExtensionsApiEnabled)
    .action((sourceLocation, options) => __awaiter(void 0, void 0, void 0, function* () {
    const projectId = getProjectId(options);
    const res = yield extensionsHelper_1.createSourceFromLocation(projectId, sourceLocation);
    logger.info(`Extension source creation successful for ${res.spec.name}! Your new source is ${res.name}`);
    return res;
}));
