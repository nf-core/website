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
const fs = require("fs");
const path = require("path");
const rimraf = require("rimraf");
const api = require("../api");
const command_1 = require("../command");
const commandUtils = require("../emulator/commandUtils");
const utils = require("../utils");
const hub_1 = require("../emulator/hub");
const error_1 = require("../error");
const hubExport_1 = require("../emulator/hubExport");
const prompt_1 = require("../prompt");
module.exports = new command_1.Command("emulators:export <path>")
    .description("export data from running emulators")
    .option(commandUtils.FLAG_ONLY, commandUtils.DESC_ONLY)
    .option("--force", "Overwrite any export data in the target directory.")
    .action((exportPath, options) => __awaiter(void 0, void 0, void 0, function* () {
    const projectId = options.project;
    if (!projectId) {
        throw new error_1.FirebaseError("Could not determine project ID, make sure you're running in a Firebase project directory or add the --project flag.", { exit: 1 });
    }
    const locator = hub_1.EmulatorHub.readLocatorFile(projectId);
    if (!locator) {
        throw new error_1.FirebaseError(`Did not find any running emulators for project ${clc.bold(projectId)}.`, { exit: 1 });
    }
    const hubOrigin = `http://${locator.host}:${locator.port}`;
    try {
        yield api.request("GET", "/", {
            origin: hubOrigin,
        });
    }
    catch (e) {
        const filePath = hub_1.EmulatorHub.getLocatorFilePath(projectId);
        throw new error_1.FirebaseError(`The emulator hub at ${hubOrigin} did not respond to a status check. If this error continues try shutting down all running emulators and deleting the file ${filePath}`, { exit: 1 });
    }
    utils.logBullet(`Found running emulator hub for project ${clc.bold(projectId)} at ${hubOrigin}`);
    const exportAbsPath = path.resolve(exportPath);
    if (!fs.existsSync(exportAbsPath)) {
        utils.logBullet(`Creating export directory ${exportAbsPath}`);
        fs.mkdirSync(exportAbsPath);
    }
    const existingMetadata = hubExport_1.HubExport.readMetadata(exportAbsPath);
    if (existingMetadata && !options.force) {
        if (options.noninteractive) {
            throw new error_1.FirebaseError("Export already exists in the target directory, re-run with --force to overwrite.", { exit: 1 });
        }
        const prompt = yield prompt_1.promptOnce({
            type: "confirm",
            message: `The directory ${exportAbsPath} already contains export data. Exporting again to the same directory will overwrite all data. Do you want to continue?`,
            default: false,
        });
        if (!prompt) {
            throw new error_1.FirebaseError("Command aborted", { exit: 1 });
        }
    }
    if (existingMetadata) {
        if (existingMetadata.firestore) {
            const firestorePath = path.join(exportAbsPath, existingMetadata.firestore.path);
            utils.logBullet(`Deleting directory ${firestorePath}`);
            rimraf.sync(firestorePath);
        }
    }
    utils.logBullet(`Exporting data to: ${exportAbsPath}`);
    yield api
        .request("POST", hub_1.EmulatorHub.PATH_EXPORT, {
        origin: hubOrigin,
        json: true,
        data: {
            path: exportAbsPath,
        },
    })
        .catch((e) => {
        throw new error_1.FirebaseError("Export request failed, see emulator logs for more information.", {
            exit: 1,
            original: e,
        });
    });
    utils.logSuccess("Export complete");
}));
