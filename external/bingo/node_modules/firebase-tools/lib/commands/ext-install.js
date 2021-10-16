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
const marked = require("marked");
const ora = require("ora");
const TerminalRenderer = require("marked-terminal");
const askUserForConsent = require("../extensions/askUserForConsent");
const checkProjectBilling = require("../extensions/checkProjectBilling");
const command_1 = require("../command");
const error_1 = require("../error");
const getProjectId = require("../getProjectId");
const rolesHelper_1 = require("../extensions/rolesHelper");
const extensionsApi = require("../extensions/extensionsApi");
const resolveSource_1 = require("../extensions/resolveSource");
const paramHelper = require("../extensions/paramHelper");
const extensionsHelper_1 = require("../extensions/extensionsHelper");
const utils_1 = require("../extensions/utils");
const requirePermissions_1 = require("../requirePermissions");
const utils = require("../utils");
const logger = require("../logger");
const prompt_1 = require("../prompt");
const previews = require("../previews");
marked.setOptions({
    renderer: new TerminalRenderer(),
});
function installExtension(options) {
    return __awaiter(this, void 0, void 0, function* () {
        const { projectId, source, paramFilePath } = options;
        const spec = source.spec;
        const spinner = ora.default("Installing your extension instance. This usually takes 3 to 5 minutes...");
        try {
            yield checkProjectBilling(projectId, spec.displayName || spec.name, spec.billingRequired);
            const roles = spec.roles ? spec.roles.map((role) => role.role) : [];
            yield askUserForConsent.prompt(spec.displayName || spec.name, projectId, roles);
            const params = yield paramHelper.getParams(projectId, _.get(spec, "params", []), paramFilePath);
            let instanceId = spec.name;
            const anotherInstanceExists = yield extensionsHelper_1.instanceIdExists(projectId, instanceId);
            if (anotherInstanceExists) {
                const consent = yield extensionsHelper_1.promptForRepeatInstance(projectId, spec.name);
                if (!consent) {
                    logger.info(marked("Installation cancelled. For a list of all available Firebase Extensions commands, run `firebase ext`."));
                    return;
                }
                instanceId = yield extensionsHelper_1.promptForValidInstanceId(`${instanceId}-${utils_1.getRandomString(4)}`);
            }
            spinner.start();
            let serviceAccountEmail;
            while (!serviceAccountEmail) {
                try {
                    serviceAccountEmail = yield rolesHelper_1.createServiceAccountAndSetRoles(projectId, _.get(spec, "roles", []), instanceId);
                }
                catch (err) {
                    if (err.status === 409) {
                        spinner.stop();
                        logger.info(err.message);
                        instanceId = yield extensionsHelper_1.promptForValidInstanceId(`${instanceId}-${utils_1.getRandomString(4)}`);
                        spinner.start();
                    }
                    else {
                        throw err;
                    }
                }
            }
            const response = yield extensionsApi.createInstance(projectId, instanceId, source, params, serviceAccountEmail);
            spinner.stop();
            utils.logLabeledSuccess(extensionsHelper_1.logPrefix, `Successfully installed your instance of ${clc.bold(spec.displayName || spec.name)}! ` +
                `Its Instance ID is ${clc.bold(instanceId)}.`);
            utils.logLabeledBullet(extensionsHelper_1.logPrefix, marked("Go to the Firebase console to view instructions for using your extension, " +
                `which may include some required post-installation tasks: ${utils.consoleUrl(projectId, `/extensions/instances/${instanceId}?tab=usage`)}`));
            logger.info(marked("You can run `firebase ext` to view available Firebase Extensions commands, " +
                "including those to update, reconfigure, or delete your installed extension."));
        }
        catch (err) {
            if (spinner.isSpinning) {
                spinner.fail();
            }
            if (err instanceof error_1.FirebaseError) {
                throw err;
            }
            throw new error_1.FirebaseError(`Error occurred installing extension: ${err.message}`, {
                original: err,
            });
        }
    });
}
exports.default = new command_1.Command("ext:install [extensionName]")
    .description("install an official extension if [extensionName] or [extensionName@version] is provided;" +
    previews.extdev
    ? "install a local extension if [localPathOrUrl] or [url#root] is provided;"
    : "" + "or run with `-i` to see all available extensions.")
    .option("--params <paramsFile>", "name of params variables file with .env format.")
    .before(requirePermissions_1.requirePermissions, ["firebasemods.instances.create"])
    .before(extensionsHelper_1.ensureExtensionsApiEnabled)
    .action((extensionName, options) => __awaiter(void 0, void 0, void 0, function* () {
    const projectId = getProjectId(options, false);
    const paramFilePath = options.params;
    let learnMore = false;
    if (!extensionName) {
        if (options.interactive) {
            learnMore = true;
            extensionName = yield extensionsHelper_1.promptForOfficialExtension("Which official extension do you want to install?\n" +
                "  Select an extension, then press Enter to learn more.");
        }
        else {
            throw new error_1.FirebaseError(`Please provide an extension name, or run ${clc.bold("firebase ext:install -i")} to select from the list of all available official extensions.`);
        }
    }
    const [name, version] = extensionName.split("@");
    let source;
    try {
        const registryEntry = yield resolveSource_1.resolveRegistryEntry(name);
        const audienceConsent = yield resolveSource_1.promptForAudienceConsent(registryEntry);
        if (!audienceConsent) {
            logger.info("Install cancelled.");
            return;
        }
        const sourceUrl = resolveSource_1.resolveSourceUrl(registryEntry, name, version);
        source = yield extensionsApi.getSource(sourceUrl);
    }
    catch (err) {
        if (previews.extdev) {
            try {
                source = yield extensionsHelper_1.createSourceFromLocation(projectId, extensionName);
            }
            catch (err) {
                throw new error_1.FirebaseError(`Unable to find official extension named ${clc.bold(extensionName)} ` +
                    `or local extension at ${clc.bold(extensionName)}` +
                    `Run ${clc.bold("firebase ext:install -i")} to select from the list of all available official extensions.`, { original: err });
            }
        }
        else {
            throw new error_1.FirebaseError(`Unable to find offical extension source named ${clc.bold(extensionName)}. ` +
                `Run ${clc.bold("firebase ext:install -i")} to select from the list of all available official extensions.`, { original: err });
        }
    }
    try {
        if (learnMore) {
            utils.logLabeledBullet(extensionsHelper_1.logPrefix, `You selected: ${clc.bold(source.spec.displayName)}.\n` +
                `${source.spec.description}\n` +
                `View details: https://firebase.google.com/products/extensions/${name}\n`);
            const confirm = yield prompt_1.promptOnce({
                type: "confirm",
                default: true,
                message: "Do you want to install this extension?",
            });
            if (!confirm) {
                return;
            }
        }
        return installExtension({
            paramFilePath,
            projectId,
            source,
        });
    }
    catch (err) {
        if (!(err instanceof error_1.FirebaseError)) {
            throw new error_1.FirebaseError(`Error occurred installing the extension: ${err.message}`, {
                original: err,
            });
        }
        throw err;
    }
}));
