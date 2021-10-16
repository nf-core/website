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
const cli_color_1 = require("cli-color");
const lodash_1 = require("lodash");
const error_1 = require("./error");
const utils_1 = require("./utils");
const rc_1 = require("./rc");
const config_1 = require("./config");
const configstore_1 = require("./configstore");
const detectProjectRoot_1 = require("./detectProjectRoot");
const logger = require("./logger");
const track = require("./track");
class Command {
    constructor(cmd) {
        this.cmd = cmd;
        this.name = "";
        this.descriptionText = "";
        this.options = [];
        this.actionFn = () => { };
        this.befores = [];
        this.helpText = "";
        this.name = lodash_1.first(cmd.split(" ")) || "";
    }
    description(t) {
        this.descriptionText = t;
        return this;
    }
    option(...args) {
        this.options.push(args);
        return this;
    }
    before(fn, ...args) {
        this.befores.push({ fn: fn, args: args });
        return this;
    }
    help(t) {
        this.helpText = t;
        return this;
    }
    action(fn) {
        this.actionFn = fn;
        return this;
    }
    register(client) {
        this.client = client;
        const program = client.cli;
        const cmd = program.command(this.cmd);
        if (this.descriptionText) {
            cmd.description(this.descriptionText);
        }
        this.options.forEach((args) => {
            const flags = args.shift();
            cmd.option(flags, ...args);
        });
        if (this.helpText) {
            cmd.on("--help", () => {
                console.log(this.helpText);
            });
        }
        cmd.action((...args) => {
            const runner = this.runner();
            const start = new Date().getTime();
            const options = lodash_1.last(args);
            if (args.length - 1 > cmd._args.length) {
                client.errorOut(new error_1.FirebaseError(`Too many arguments. Run ${cli_color_1.bold("firebase help " + this.name)} for usage instructions`, { exit: 1 }));
                return;
            }
            runner(...args)
                .then((result) => {
                if (utils_1.getInheritedOption(options, "json")) {
                    console.log(JSON.stringify({
                        status: "success",
                        result: result,
                    }, null, 2));
                }
                const duration = new Date().getTime() - start;
                track(this.name, "success", duration).then(process.exit);
            })
                .catch((err) => __awaiter(this, void 0, void 0, function* () {
                if (utils_1.getInheritedOption(options, "json")) {
                    console.log(JSON.stringify({
                        status: "error",
                        error: err.message,
                    }, null, 2));
                }
                const duration = Date.now() - start;
                const errorEvent = err.exit === 1 ? "Error (User)" : "Error (Unexpected)";
                yield Promise.all([track(this.name, "error", duration), track(errorEvent, "", duration)]);
                client.errorOut(err);
            }));
        });
    }
    prepare(options) {
        options = options || {};
        options.project = utils_1.getInheritedOption(options, "project");
        if (!process.stdin.isTTY || utils_1.getInheritedOption(options, "nonInteractive")) {
            options.nonInteractive = true;
        }
        if (utils_1.getInheritedOption(options, "interactive")) {
            options.interactive = true;
            options.nonInteractive = false;
        }
        if (utils_1.getInheritedOption(options, "debug")) {
            logger.transports.console.level = "debug";
            options.debug = true;
        }
        if (utils_1.getInheritedOption(options, "json")) {
            options.nonInteractive = true;
            logger.transports.console.level = "none";
        }
        try {
            options.config = config_1.load(options);
        }
        catch (e) {
            options.configError = e;
        }
        options.projectRoot = detectProjectRoot_1.detectProjectRoot(options.cwd);
        this.applyRC(options);
    }
    applyRC(options) {
        const rc = rc_1.load(options.cwd);
        options.rc = rc;
        options.project =
            options.project || (configstore_1.configstore.get("activeProjects") || {})[options.projectRoot];
        if (options.config && !options.project) {
            options.project = options.config.defaults.project;
        }
        const aliases = rc.projects;
        const rcProject = lodash_1.get(aliases, options.project);
        if (rcProject) {
            options.projectAlias = options.project;
            options.project = rcProject;
        }
        else if (!options.project && lodash_1.size(aliases) === 1) {
            options.projectAlias = lodash_1.head(lodash_1.keys(aliases));
            options.project = lodash_1.head(lodash_1.values(aliases));
        }
    }
    runner() {
        return (...args) => __awaiter(this, void 0, void 0, function* () {
            if (args.length === 0) {
                args.push({});
            }
            const options = lodash_1.last(args);
            this.prepare(options);
            for (const before of this.befores) {
                yield before.fn(options, ...before.args);
            }
            return this.actionFn(...args);
        });
    }
}
exports.Command = Command;
