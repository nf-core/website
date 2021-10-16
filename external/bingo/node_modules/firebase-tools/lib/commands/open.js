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
const open = require("open");
const error_1 = require("../error");
const api = require("../api");
const command_1 = require("../command");
const logger = require("../logger");
const prompt_1 = require("../prompt");
const requirePermissions_1 = require("../requirePermissions");
const requireInstance = require("../requireInstance");
const utils = require("../utils");
const LINKS = [
    { name: "Analytics", arg: "analytics", consolePath: "/analytics" },
    { name: "Authentication: Providers", arg: "auth", consolePath: "/authentication/providers" },
    { name: "Authentication: Users", arg: "auth:users", consolePath: "/authentication/users" },
    { name: "Crash Reporting", arg: "crash", consolePath: "/monitoring" },
    { name: "Database: Data", arg: "database", consolePath: "/database/data" },
    { name: "Database: Rules", arg: "database:rules", consolePath: "/database/rules" },
    { name: "Docs", arg: "docs", url: "https://firebase.google.com/docs" },
    { name: "Dynamic Links", arg: "links", consolePath: "/durablelinks" },
    { name: "Functions", arg: "functions", consolePath: "/functions/list" },
    { name: "Functions Log", arg: "functions:log" },
    { name: "Hosting: Deployed Site", arg: "hosting:site" },
    { name: "Hosting", arg: "hosting", consolePath: "/hosting/main" },
    { name: "Notifications", arg: "notifications", consolePath: "/notification" },
    { name: "Project Dashboard", arg: "dashboard", consolePath: "/overview" },
    { name: "Project Settings", arg: "settings", consolePath: "/settings/general" },
    {
        name: "Remote Config: Conditions",
        arg: "config:conditions",
        consolePath: "/config/conditions",
    },
    { name: "Remote Config", arg: "config", consolePath: "/config" },
    { name: "Storage: Files", arg: "storage", consolePath: "/storage/files" },
    { name: "Storage: Rules", arg: "storage:rules", consolePath: "/storage/rules" },
    { name: "Test Lab", arg: "testlab", consolePath: "/testlab/histories/" },
];
const CHOICES = _.map(LINKS, "name");
exports.default = new command_1.Command("open [link]")
    .description("quickly open a browser to relevant project resources")
    .before(requirePermissions_1.requirePermissions)
    .before(requireInstance)
    .action((linkName, options) => __awaiter(void 0, void 0, void 0, function* () {
    let link = _.find(LINKS, { arg: linkName });
    if (linkName && !link) {
        throw new error_1.FirebaseError("Unrecognized link name. Valid links are:\n\n" + _.map(LINKS, "arg").join("\n"));
    }
    if (!link) {
        const name = yield prompt_1.promptOnce({
            type: "list",
            message: "What link would you like to open?",
            choices: CHOICES,
        });
        link = _.find(LINKS, { name });
    }
    if (!link) {
        throw new error_1.FirebaseError("Unrecognized link name. Valid links are:\n\n" + _.map(LINKS, "arg").join("\n"));
    }
    let url;
    if (link.consolePath) {
        url = utils.consoleUrl(options.project, link.consolePath);
    }
    else if (link.url) {
        url = link.url;
    }
    else if (link.arg === "hosting:site") {
        url = utils.addSubdomain(api.hostingOrigin, options.instance);
    }
    else if (link.arg === "functions:log") {
        url = `https://console.developers.google.com/logs/viewer?resource=cloudfunctions.googleapis.com&project=${options.project}`;
    }
    else {
        throw new error_1.FirebaseError(`Unable to determine URL for link: ${link}`);
    }
    if (link.arg !== linkName) {
        logger.info(`${clc.bold.cyan("Tip:")} You can also run ${clc.bold.underline(`firebase open ${link.arg}`)}`);
        logger.info();
    }
    logger.info(`Opening ${clc.bold(link.name)} link in your default browser:`);
    logger.info(clc.bold.underline(url));
    open(url);
}));
