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
const clc = require("cli-color");
const _ = require("lodash");
const command_1 = require("../command");
const resolveSource_1 = require("../extensions/resolveSource");
const extensionsApi = require("../extensions/extensionsApi");
const extensionsHelper_1 = require("../extensions/extensionsHelper");
const localHelper_1 = require("../extensions/localHelper");
const logger = require("../logger");
const requirePermissions_1 = require("../requirePermissions");
const utils = require("../utils");
const marked = require("marked");
const TerminalRenderer = require("marked-terminal");
const FUNCTION_TYPE_REGEX = /\..+\.function/;
exports.default = new command_1.Command("ext:info <extensionName>")
    .description("display information about an extension by name (extensionName@x.y.z for a specific version)")
    .option("--markdown", "output info in Markdown suitable for constructing a README file")
    .action((extensionName, options) => __awaiter(void 0, void 0, void 0, function* () {
    let spec;
    if (localHelper_1.isLocalExtension(extensionName)) {
        if (!options.markdown) {
            utils.logLabeledBullet(extensionsHelper_1.logPrefix, `reading extension from directory: ${extensionName}`);
        }
        spec = yield localHelper_1.getLocalExtensionSpec(extensionName);
    }
    else {
        yield requirePermissions_1.requirePermissions(options, ["firebasemods.sources.get"]);
        yield extensionsHelper_1.ensureExtensionsApiEnabled(options);
        const [name, version] = extensionName.split("@");
        const registryEntry = yield resolveSource_1.resolveRegistryEntry(name);
        const sourceUrl = yield resolveSource_1.resolveSourceUrl(registryEntry, name, version);
        const source = yield extensionsApi.getSource(sourceUrl);
        spec = source.spec;
    }
    if (!options.markdown) {
        utils.logLabeledBullet(extensionsHelper_1.logPrefix, `information about ${extensionName}:\n`);
    }
    const lines = [];
    if (options.markdown) {
        lines.push(`# ${spec.displayName}`);
    }
    else {
        lines.push(`**Name**: ${spec.displayName}`);
    }
    if (spec.description) {
        lines.push(`**Description**: ${spec.description}`);
    }
    if (spec.preinstallContent) {
        lines.push("", `**Details**: ${spec.preinstallContent}`);
    }
    if (spec.params && Array.isArray(spec.params) && spec.params.length > 0) {
        lines.push("", "**Configuration Parameters:**");
        _.forEach(spec.params, (param) => {
            lines.push(`* ${param.label}` + (param.description ? `: ${param.description}` : ""));
        });
    }
    const functions = [];
    const otherResources = [];
    _.forEach(spec.resources, (resource) => {
        if (FUNCTION_TYPE_REGEX.test(resource.type)) {
            functions.push(resource);
        }
        else {
            otherResources.push(resource);
        }
    });
    if (functions.length > 0) {
        lines.push("", "**Cloud Functions:**");
        _.forEach(functions, (func) => {
            lines.push(`* **${func.name}:** ${func.description}`);
        });
    }
    if (otherResources.length > 0) {
        lines.push("", "**Other Resources**:");
        _.forEach(otherResources, (resource) => {
            lines.push(`* ${resource.name} (${resource.type})`);
        });
    }
    if (spec.apis) {
        lines.push("", "**APIs Used**:");
        _.forEach(spec.apis, (api) => {
            lines.push(`* ${api.apiName}` + (api.reason ? ` (Reason: ${api.reason})` : ""));
        });
    }
    if (spec.roles) {
        lines.push("", "**Access Required**:");
        lines.push("", "This extension will operate with the following project IAM roles:");
        _.forEach(spec.roles, (role) => {
            lines.push(`* ${role.role}` + (role.reason ? ` (Reason: ${role.reason})` : ""));
        });
    }
    if (options.markdown) {
        logger.info(lines.join("\n\n"));
    }
    else {
        marked.setOptions({
            renderer: new TerminalRenderer(),
        });
        logger.info(marked(lines.join("\n")));
        utils.logLabeledBullet(extensionsHelper_1.logPrefix, `to install this extension, type ` +
            clc.bold(`firebase ext:install ${extensionName} --project=YOUR_PROJECT`));
    }
}));
