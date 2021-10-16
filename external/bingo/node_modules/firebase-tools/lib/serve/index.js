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
const logger = require("../logger");
const TARGETS = {
    hosting: require("./hosting"),
    functions: require("./functions"),
    database: require("./database"),
    firestore: require("./firestore"),
};
function serve(options) {
    return __awaiter(this, void 0, void 0, function* () {
        const targetNames = options.targets;
        options.port = parseInt(options.port, 10);
        yield Promise.all(_.map(targetNames, (targetName) => {
            return TARGETS[targetName].start(options);
        }));
        yield Promise.all(_.map(targetNames, (targetName) => {
            return TARGETS[targetName].connect();
        }));
        yield new Promise((resolve) => {
            process.on("SIGINT", () => {
                logger.info("Shutting down...");
                return Promise.all(_.map(targetNames, (targetName) => {
                    return TARGETS[targetName].stop(options);
                }))
                    .then(resolve)
                    .catch(resolve);
            });
        });
    });
}
exports.serve = serve;
