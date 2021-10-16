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
const command_1 = require("../command");
const commandUtils = require("../emulator/commandUtils");
const optionsHelper = require("../extensions/emulator/optionsHelper");
module.exports = new command_1.Command("ext:dev:emulators:exec <script>")
    .description("emulate an extension, run a test script, then shut down the emulators")
    .option(commandUtils.FLAG_INSPECT_FUNCTIONS, commandUtils.DESC_INSPECT_FUNCTIONS)
    .option(commandUtils.FLAG_TEST_CONFIG, commandUtils.DESC_TEST_CONFIG)
    .option(commandUtils.FLAG_TEST_PARAMS, commandUtils.DESC_TEST_PARAMS)
    .action((script, options) => __awaiter(void 0, void 0, void 0, function* () {
    const emulatorOptions = yield optionsHelper.buildOptions(options);
    commandUtils.beforeEmulatorCommand(emulatorOptions);
    yield commandUtils.emulatorExec(script, emulatorOptions);
}));
