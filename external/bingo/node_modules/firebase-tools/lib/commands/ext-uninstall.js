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
const ora = require("ora");
const command_1 = require("../command");
const error_1 = require("../error");
const getProjectId = require("../getProjectId");
const gcp_1 = require("../gcp");
const extensionsApi = require("../extensions/extensionsApi");
const extensionsHelper_1 = require("../extensions/extensionsHelper");
const prompt_1 = require("../prompt");
const requirePermissions_1 = require("../requirePermissions");
const utils = require("../utils");
const logger = require("../logger");
exports.default = new command_1.Command("ext:uninstall <extensionInstanceId>")
    .description("uninstall an extension that is installed in your Firebase project by instance ID")
    .option("-f, --force", "No confirmation. Otherwise, a confirmation prompt will appear.")
    .before(requirePermissions_1.requirePermissions, ["firebasemods.instances.delete"])
    .before(extensionsHelper_1.ensureExtensionsApiEnabled)
    .action((instanceId, options) => __awaiter(void 0, void 0, void 0, function* () {
    const projectId = getProjectId(options);
    let instance;
    try {
        instance = yield extensionsApi.getInstance(projectId, instanceId);
    }
    catch (err) {
        if (err.status === 404) {
            return utils.reject(`No extension instance ${instanceId} in project ${projectId}.`, {
                exit: 1,
            });
        }
        throw err;
    }
    if (!options.force) {
        const serviceAccountMessage = `Uninstalling deletes the service account used by this extension instance:\n${clc.bold(instance.serviceAccountEmail)}\n\n`;
        const resourcesMessage = _.get(instance, "config.source.spec.resources", []).length
            ? "Uninstalling deletes all extension resources created for this extension instance:\n" +
                instance.config.source.spec.resources
                    .map((resource) => clc.bold(`- ${extensionsHelper_1.resourceTypeToNiceName[resource.type] || resource.type}: ${resource.name} \n`))
                    .join("") +
                "\n"
            : "";
        const artifactsMessage = `The following ${clc.bold("will not")} be deleted:\n` +
            "Any artifacts (for example, stored images) created by this extension instance.\n" +
            "Any other project resources with which this extension instance interacted.\n";
        const extensionDeletionMessage = `Here's what will happen when you uninstall ${clc.bold(instanceId)} from project ${clc.bold(projectId)}. Be aware that this cannot be undone.\n\n` +
            `${serviceAccountMessage}` +
            `${resourcesMessage}` +
            `${artifactsMessage}`;
        logger.info(extensionDeletionMessage);
        const confirmedExtensionDeletion = yield prompt_1.promptOnce({
            type: "confirm",
            default: true,
            message: "Are you sure that you want to uninstall this extension?",
        });
        if (!confirmedExtensionDeletion) {
            return utils.reject("Command aborted.", { exit: 1 });
        }
    }
    const spinner = ora.default(` ${clc.green.bold(extensionsHelper_1.logPrefix)}: uninstalling ${clc.bold(instanceId)}. This usually takes 1 to 2 minutes...`);
    spinner.start();
    try {
        spinner.info();
        spinner.text = ` ${clc.green.bold(extensionsHelper_1.logPrefix)}: deleting your extension instance's resources.`;
        spinner.start();
        yield extensionsApi.deleteInstance(projectId, instanceId);
        spinner.succeed(` ${clc.green.bold(extensionsHelper_1.logPrefix)}: deleted your extension instance's resources.`);
        spinner.text = ` ${clc.green.bold(extensionsHelper_1.logPrefix)}: deleting your extension instance's service account.`;
        spinner.start();
        const saDeletionRes = yield gcp_1.iam.deleteServiceAccount(projectId, instance.serviceAccountEmail);
        if (_.get(saDeletionRes, "body.error")) {
            if (_.get(saDeletionRes, "body.error.code") === 404) {
                spinner.succeed(` ${clc.green.bold(extensionsHelper_1.logPrefix)}: service account ${clc.bold(instance.serviceAccountEmail)} was previously deleted.`);
            }
            else {
                throw new error_1.FirebaseError("Unable to delete service account", {
                    original: saDeletionRes.body.error,
                });
            }
        }
        else {
            spinner.succeed(` ${clc.green.bold(extensionsHelper_1.logPrefix)}: deleted your extension instance's service account ${clc.bold(instance.serviceAccountEmail)}`);
        }
    }
    catch (err) {
        if (spinner.isSpinning) {
            spinner.fail();
        }
        if (err instanceof error_1.FirebaseError) {
            throw err;
        }
        return utils.reject(`Error occurred uninstalling extension ${instanceId}`, { original: err });
    }
    utils.logLabeledSuccess(extensionsHelper_1.logPrefix, `uninstalled ${clc.bold(instanceId)}`);
}));
