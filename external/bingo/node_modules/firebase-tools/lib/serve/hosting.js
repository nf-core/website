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
const superstatic = require("superstatic").server;
const detectProjectRoot_1 = require("../detectProjectRoot");
const error_1 = require("../error");
const implicitInit_1 = require("../hosting/implicitInit");
const initMiddleware_1 = require("../hosting/initMiddleware");
const normalizedHostingConfigs_1 = require("../hosting/normalizedHostingConfigs");
const utils = require("../utils");
const cloudRunProxy_1 = require("../hosting/cloudRunProxy");
const functionsProxy_1 = require("../hosting/functionsProxy");
const MAX_PORT_ATTEMPTS = 10;
let attempts = 0;
let server;
function startServer(options, config, port, init) {
    server = superstatic({
        debug: true,
        port: port,
        host: options.host,
        config: config,
        cwd: detectProjectRoot_1.detectProjectRoot(options.cwd),
        stack: "strict",
        before: {
            files: initMiddleware_1.initMiddleware(init),
        },
        rewriters: {
            function: functionsProxy_1.default(options),
            run: cloudRunProxy_1.default(options),
        },
    }).listen(() => {
        const siteName = config.target || config.site;
        const label = siteName ? "hosting[" + siteName + "]" : "hosting";
        if (config.public && config.public !== ".") {
            utils.logLabeledBullet(label, "Serving hosting files from: " + clc.bold(config.public));
        }
        utils.logLabeledSuccess(label, "Local server: " + clc.underline(clc.bold("http://" + options.host + ":" + port)));
    });
    server.on("error", (err) => {
        if (err.code === "EADDRINUSE") {
            const message = "Port " + options.port + " is not available.";
            if (attempts < MAX_PORT_ATTEMPTS) {
                utils.logWarning(clc.yellow("hosting: ") + message + " Trying another port...");
                attempts++;
                startServer(options, config, port + 5, init);
            }
            else {
                utils.logWarning(message);
                throw new error_1.FirebaseError("Could not find an open port for hosting development server.", {
                    exit: 1,
                });
            }
        }
        else {
            throw new error_1.FirebaseError("An error occurred while starting the hosting development server:\n\n" + err.toString(), { exit: 1 });
        }
    });
}
function stop() {
    return __awaiter(this, void 0, void 0, function* () {
        if (server) {
            yield server.close();
        }
    });
}
exports.stop = stop;
function start(options) {
    return __awaiter(this, void 0, void 0, function* () {
        const init = yield implicitInit_1.implicitInit(options);
        const configs = normalizedHostingConfigs_1.normalizedHostingConfigs(options);
        for (let i = 0; i < configs.length; i++) {
            const port = i === 0 ? options.port : options.port + 4 + i;
            startServer(options, configs[i], port, init);
        }
    });
}
exports.start = start;
function connect() {
    return __awaiter(this, void 0, void 0, function* () {
        yield Promise.resolve();
    });
}
exports.connect = connect;
