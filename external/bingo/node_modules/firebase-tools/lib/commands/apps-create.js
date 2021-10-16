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
const ora = require("ora");
const command_1 = require("../command");
const getProjectId = require("../getProjectId");
const error_1 = require("../error");
const apps_1 = require("../management/apps");
const prompt_1 = require("../prompt");
const requireAuth = require("../requireAuth");
const logger = require("../logger");
const DISPLAY_NAME_QUESTION = {
    type: "input",
    name: "displayName",
    default: "",
    message: "What would you like to call your app?",
};
function logPostAppCreationInformation(appMetadata, appPlatform) {
    logger.info("");
    logger.info(`🎉🎉🎉 Your Firebase ${appPlatform} App is ready! 🎉🎉🎉`);
    logger.info("");
    logger.info("App information:");
    logger.info(`  - App ID: ${appMetadata.appId}`);
    if (appMetadata.displayName) {
        logger.info(`  - Display name: ${appMetadata.displayName}`);
    }
    if (appPlatform === apps_1.AppPlatform.IOS) {
        const iosAppMetadata = appMetadata;
        logger.info(`  - Bundle ID: ${iosAppMetadata.bundleId}`);
        if (iosAppMetadata.appStoreId) {
            logger.info(`  - App Store ID: ${iosAppMetadata.appStoreId}`);
        }
    }
    else if (appPlatform === apps_1.AppPlatform.ANDROID) {
        logger.info(`  - Package name: ${appMetadata.packageName}`);
    }
    logger.info("");
    logger.info("You can run this command to print out your new app's Google Services config:");
    logger.info(`  firebase apps:sdkconfig ${appPlatform} ${appMetadata.appId}`);
}
function initiateIosAppCreation(options) {
    return __awaiter(this, void 0, void 0, function* () {
        if (!options.nonInteractive) {
            yield prompt_1.prompt(options, [
                DISPLAY_NAME_QUESTION,
                {
                    type: "input",
                    default: "",
                    name: "bundleId",
                    message: "Please specify your iOS app bundle ID:",
                },
                {
                    type: "input",
                    default: "",
                    name: "appStoreId",
                    message: "Please specify your iOS app App Store ID:",
                },
            ]);
        }
        if (!options.bundleId) {
            throw new error_1.FirebaseError("Bundle ID for iOS app cannot be empty");
        }
        const spinner = ora("Creating your iOS app").start();
        try {
            const appData = yield apps_1.createIosApp(options.project, {
                displayName: options.displayName,
                bundleId: options.bundleId,
                appStoreId: options.appStoreId,
            });
            spinner.succeed();
            return appData;
        }
        catch (err) {
            spinner.fail();
            throw err;
        }
    });
}
function initiateAndroidAppCreation(options) {
    return __awaiter(this, void 0, void 0, function* () {
        if (!options.nonInteractive) {
            yield prompt_1.prompt(options, [
                DISPLAY_NAME_QUESTION,
                {
                    type: "input",
                    default: "",
                    name: "packageName",
                    message: "Please specify your Android app package name:",
                },
            ]);
        }
        if (!options.packageName) {
            throw new error_1.FirebaseError("Package name for Android app cannot be empty");
        }
        const spinner = ora("Creating your Android app").start();
        try {
            const appData = yield apps_1.createAndroidApp(options.project, {
                displayName: options.displayName,
                packageName: options.packageName,
            });
            spinner.succeed();
            return appData;
        }
        catch (err) {
            spinner.fail();
            throw err;
        }
    });
}
function initiateWebAppCreation(options) {
    return __awaiter(this, void 0, void 0, function* () {
        if (!options.nonInteractive) {
            yield prompt_1.prompt(options, [DISPLAY_NAME_QUESTION]);
        }
        if (!options.displayName) {
            throw new error_1.FirebaseError("Display name for Web app cannot be empty");
        }
        const spinner = ora("Creating your Web app").start();
        try {
            const appData = yield apps_1.createWebApp(options.project, { displayName: options.displayName });
            spinner.succeed();
            return appData;
        }
        catch (err) {
            spinner.fail();
            throw err;
        }
    });
}
module.exports = new command_1.Command("apps:create [platform] [displayName]")
    .description("create a new Firebase app. [platform] can be IOS, ANDROID or WEB (case insensitive).")
    .option("-a, --package-name <packageName>", "required package name for the Android app")
    .option("-b, --bundle-id <bundleId>", "required bundle id for the iOS app")
    .option("-s, --app-store-id <appStoreId>", "(optional) app store id for the iOS app")
    .before(requireAuth)
    .action((platform = "", displayName, options) => __awaiter(void 0, void 0, void 0, function* () {
    const projectId = getProjectId(options);
    if (!options.nonInteractive && !platform) {
        platform = yield prompt_1.promptOnce({
            type: "list",
            message: "Please choose the platform of the app:",
            choices: [
                { name: "iOS", value: apps_1.AppPlatform.IOS },
                { name: "Android", value: apps_1.AppPlatform.ANDROID },
                { name: "Web", value: apps_1.AppPlatform.WEB },
            ],
        });
    }
    const appPlatform = apps_1.getAppPlatform(platform);
    if (appPlatform === apps_1.AppPlatform.ANY) {
        throw new error_1.FirebaseError("App platform must be provided");
    }
    logger.info(`Create your ${appPlatform} app in project ${clc.bold(projectId)}:`);
    options.displayName = displayName;
    let appData;
    switch (appPlatform) {
        case apps_1.AppPlatform.IOS:
            appData = yield initiateIosAppCreation(options);
            break;
        case apps_1.AppPlatform.ANDROID:
            appData = yield initiateAndroidAppCreation(options);
            break;
        case apps_1.AppPlatform.WEB:
            appData = yield initiateWebAppCreation(options);
            break;
        default:
            throw new error_1.FirebaseError("Unexpected error. This should not happen");
    }
    logPostAppCreationInformation(appData, appPlatform);
    return appData;
}));
