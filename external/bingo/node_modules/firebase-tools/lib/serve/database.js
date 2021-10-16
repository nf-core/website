"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const databaseEmulator_1 = require("../emulator/databaseEmulator");
const emulatorServer_1 = require("../emulator/emulatorServer");
module.exports = new emulatorServer_1.EmulatorServer(new databaseEmulator_1.DatabaseEmulator({}));
