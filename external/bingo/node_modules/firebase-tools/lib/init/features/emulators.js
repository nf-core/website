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
const utils = require("../../utils");
const prompt_1 = require("../../prompt");
const types_1 = require("../../emulator/types");
const constants_1 = require("../../emulator/constants");
const downloadableEmulators_1 = require("../../emulator/downloadableEmulators");
function doSetup(setup, config) {
    return __awaiter(this, void 0, void 0, function* () {
        const choices = types_1.ALL_SERVICE_EMULATORS.map((e) => {
            return {
                value: e,
                name: _.capitalize(e),
                checked: config && config.has(e),
            };
        });
        const selections = {};
        yield prompt_1.prompt(selections, [
            {
                type: "checkbox",
                name: "emulators",
                message: "Which Firebase emulators do you want to set up? " +
                    "Press Space to select emulators, then Enter to confirm your choices.",
                choices: choices,
            },
        ]);
        if (!selections.emulators) {
            return;
        }
        setup.config.emulators = setup.config.emulators || {};
        for (const selected of selections.emulators) {
            setup.config.emulators[selected] = setup.config.emulators[selected] || {};
            const currentPort = setup.config.emulators[selected].port;
            if (currentPort) {
                utils.logBullet(`Port for ${selected} already configured: ${clc.cyan(currentPort)}`);
            }
            else {
                yield prompt_1.prompt(setup.config.emulators[selected], [
                    {
                        type: "input",
                        name: "port",
                        message: `Which port do you want to use for the ${clc.underline(selected)} emulator?`,
                        default: constants_1.Constants.getDefaultPort(selected),
                    },
                ]);
            }
        }
        if (selections.emulators.length) {
            yield prompt_1.prompt(selections, [
                {
                    name: "download",
                    type: "confirm",
                    message: "Would you like to download the emulators now?",
                    default: false,
                },
            ]);
        }
        if (selections.download) {
            for (const selected of selections.emulators) {
                if (types_1.isDownloadableEmulator(selected)) {
                    yield downloadableEmulators_1.downloadIfNecessary(selected);
                }
            }
        }
    });
}
exports.doSetup = doSetup;
