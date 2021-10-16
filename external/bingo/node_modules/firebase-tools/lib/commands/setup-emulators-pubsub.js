"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const command_1 = require("../command");
const downloadEmulator = require("../emulator/download");
const EMULATOR_NAME = "pubsub";
module.exports = new command_1.Command(`setup:emulators:${EMULATOR_NAME}`)
    .description(`downloads the ${EMULATOR_NAME} emulator`)
    .action(() => {
    return downloadEmulator(EMULATOR_NAME);
});
