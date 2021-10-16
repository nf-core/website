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
const ora = require("ora");
const fs = require("fs-extra");
const command_1 = require("../command");
const apps_1 = require("../management/apps");
const projects_1 = require("../management/projects");
const error_1 = require("../error");
const requireAuth = require("../requireAuth");
const logger = require("../logger");
const prompt_1 = require("../prompt");
function selectAppInteractively(projectId, appPlatform) {
    return __awaiter(this, void 0, void 0, function* () {
        if (!projectId) {
            throw new error_1.FirebaseError("Project ID must not be empty.");
        }
        const apps = yield apps_1.listFirebaseApps(projectId, appPlatform);
        if (apps.length === 0) {
            throw new error_1.FirebaseError(`There are no ${appPlatform === apps_1.AppPlatform.ANY ? "" : appPlatform + " "}apps ` +
                "associated with this Firebase project");
        }
        const choices = apps.map((app) => {
            return {
                name: `${app.displayName || app.bundleId || app.packageName}` +
                    ` - ${app.appId} (${app.platform})`,
                value: app,
            };
        });
        return yield prompt_1.promptOnce({
            type: "list",
            name: "id",
            message: `Select the ${appPlatform === apps_1.AppPlatform.ANY ? "" : appPlatform + " "}` +
                "app to get the configuration data:",
            choices,
        });
    });
}
module.exports = new command_1.Command("apps:sdkconfig [platform] [appId]")
    .description("print the Google Services config of a Firebase app. " +
    "[platform] can be IOS, ANDROID or WEB (case insensitive)")
    .option("-o, --out [file]", "(optional) write config output to a file")
    .before(requireAuth)
    .action((platform = "", appId = "", options) => __awaiter(void 0, void 0, void 0, function* () {
    let appPlatform = apps_1.getAppPlatform(platform);
    if (!appId) {
        if (options.nonInteractive) {
            throw new error_1.FirebaseError("App ID must not be empty.");
        }
        const { projectId } = yield projects_1.getOrPromptProject(options);
        const appMetadata = yield selectAppInteractively(projectId, appPlatform);
        appId = appMetadata.appId;
        appPlatform = appMetadata.platform;
    }
    let configData;
    const spinner = ora(`Downloading configuration data of your Firebase ${appPlatform} app`).start();
    try {
        configData = yield apps_1.getAppConfig(appId, appPlatform);
    }
    catch (err) {
        spinner.fail();
        throw err;
    }
    spinner.succeed();
    const fileInfo = apps_1.getAppConfigFile(configData, appPlatform);
    if (appPlatform == apps_1.AppPlatform.WEB) {
        fileInfo.sdkConfig = configData;
    }
    if (options.out === undefined) {
        logger.info(fileInfo.fileContents);
        return fileInfo;
    }
    const shouldUseDefaultFilename = options.out === true || options.out === "";
    const filename = shouldUseDefaultFilename ? configData.fileName : options.out;
    if (fs.existsSync(filename)) {
        if (options.nonInteractive) {
            throw new error_1.FirebaseError(`${filename} already exists`);
        }
        const overwrite = yield prompt_1.promptOnce({
            type: "confirm",
            default: false,
            message: `${filename} already exists. Do you want to overwrite?`,
        });
        if (!overwrite) {
            return configData;
        }
    }
    fs.writeFileSync(filename, configData.fileContents);
    logger.info(`App configuration is written in ${filename}`);
    return configData;
}));
