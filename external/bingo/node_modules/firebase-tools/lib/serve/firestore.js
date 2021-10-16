"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const firestoreEmulator_1 = require("../emulator/firestoreEmulator");
const emulatorServer_1 = require("../emulator/emulatorServer");
module.exports = new emulatorServer_1.EmulatorServer(new firestoreEmulator_1.FirestoreEmulator({}));
