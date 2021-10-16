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
const marked = require("marked");
const ora = require("ora");
const command_1 = require("../command");
const error_1 = require("../error");
const extensionsApi = require("../extensions/extensionsApi");
const extensionsHelper_1 = require("../extensions/extensionsHelper");
const paramHelper = require("../extensions/paramHelper");
const resolveSource = require("../extensions/resolveSource");
const updateHelper_1 = require("../extensions/updateHelper");
const getProjectId = require("../getProjectId");
const requirePermissions_1 = require("../requirePermissions");
const utils = require("../utils");
const TerminalRenderer = require("marked-terminal");
marked.setOptions({
    renderer: new TerminalRenderer(),
});
exports.default = new command_1.Command("ext:update <extensionInstanceId>")
    .description("update an existing extension instance to the latest version")
    .before(requirePermissions_1.requirePermissions, ["firebasemods.instances.update", "firebasemods.instances.get"])
    .before(extensionsHelper_1.ensureExtensionsApiEnabled)
    .action((instanceId, options) => __awaiter(void 0, void 0, void 0, function* () {
    const spinner = ora.default(`Updating ${clc.bold(instanceId)}. This usually takes 3 to 5 minutes...`);
    try {
        const projectId = getProjectId(options, false);
        let existingInstance;
        try {
            existingInstance = yield extensionsApi.getInstance(projectId, instanceId);
        }
        catch (err) {
            if (err.status === 404) {
                return utils.reject(`No extension instance ${instanceId} found in project ${projectId}.`, {
                    exit: 1,
                });
            }
            throw err;
        }
        const currentSpec = _.get(existingInstance, "config.source.spec");
        const currentParams = _.get(existingInstance, "config.params");
        const registryEntry = yield resolveSource.resolveRegistryEntry(currentSpec.name);
        const targetVersion = resolveSource.getTargetVersion(registryEntry, "latest");
        utils.logLabeledBullet(extensionsHelper_1.logPrefix, `Updating ${instanceId} from version ${clc.bold(currentSpec.version)} to version ${clc.bold(targetVersion)}`);
        yield resolveSource.promptForUpdateWarnings(registryEntry, currentSpec.version, targetVersion);
        const sourceUrl = resolveSource.resolveSourceUrl(registryEntry, currentSpec.name, targetVersion);
        const newSource = yield extensionsApi.getSource(sourceUrl);
        const newSpec = newSource.spec;
        if (currentSpec.version === newSpec.version) {
            utils.logLabeledBullet(extensionsHelper_1.logPrefix, `${clc.bold(instanceId)} is already up to date. Its version is ${clc.bold(currentSpec.version)}.`);
            return;
        }
        yield updateHelper_1.displayChanges(currentSpec, newSpec);
        const newParams = yield paramHelper.promptForNewParams(currentSpec, newSpec, currentParams, projectId);
        const rolesToRemove = _.differenceWith(currentSpec.roles, _.get(newSpec, "roles", []), _.isEqual);
        spinner.start();
        const updateOptions = {
            projectId,
            instanceId,
            source: newSource,
            rolesToAdd: _.get(newSpec, "roles", []),
            rolesToRemove,
            serviceAccountEmail: existingInstance.serviceAccountEmail,
            billingRequired: newSpec.billingRequired,
        };
        if (!_.isEqual(newParams, currentParams)) {
            updateOptions.params = newParams;
        }
        yield updateHelper_1.update(updateOptions);
        spinner.stop();
        utils.logLabeledSuccess(extensionsHelper_1.logPrefix, `successfully updated ${clc.bold(instanceId)}.`);
        utils.logLabeledBullet(extensionsHelper_1.logPrefix, marked(`You can view your updated instance in the Firebase console: ${utils.consoleUrl(projectId, `/extensions/instances/${instanceId}?tab=usage`)}`));
    }
    catch (err) {
        if (spinner.isSpinning) {
            spinner.fail();
        }
        if (!(err instanceof error_1.FirebaseError)) {
            throw new error_1.FirebaseError(`Error occurred while updating the instance: ${err.message}`, {
                original: err,
            });
        }
        throw err;
    }
}));
