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
const error_1 = require("../../../error");
const iv2 = require("../../../firestore/indexes");
const fsutils = require("../../../fsutils");
const prompt_1 = require("../../../prompt");
const logger = require("../../../logger");
const indexes = new iv2.FirestoreIndexes();
const INDEXES_TEMPLATE = fs.readFileSync(__dirname + "/../../../../templates/init/firestore/firestore.indexes.json", "utf8");
function initIndexes(setup, config) {
    return __awaiter(this, void 0, void 0, function* () {
        logger.info();
        logger.info("Firestore indexes allow you to perform complex queries while");
        logger.info("maintaining performance that scales with the size of the result");
        logger.info("set. You can keep index definitions in your project directory");
        logger.info("and publish them with " + clc.bold("firebase deploy") + ".");
        logger.info();
        return prompt_1.prompt(setup.config.firestore, [
            {
                type: "input",
                name: "indexes",
                message: "What file should be used for Firestore indexes?",
                default: "firestore.indexes.json",
            },
        ])
            .then(() => {
            const filename = setup.config.firestore.indexes;
            if (fsutils.fileExistsSync(filename)) {
                const msg = "File " +
                    clc.bold(filename) +
                    " already exists." +
                    " Do you want to overwrite it with the Firestore Indexes from the Firebase Console?";
                return prompt_1.promptOnce({
                    type: "confirm",
                    message: msg,
                    default: false,
                });
            }
            return Promise.resolve(true);
        })
            .then((overwrite) => {
            if (!overwrite) {
                return Promise.resolve();
            }
            return getIndexesFromConsole(setup.projectId).then((contents) => {
                return config.writeProjectFile(setup.config.firestore.indexes, contents);
            });
        });
    });
}
exports.initIndexes = initIndexes;
function getIndexesFromConsole(projectId) {
    return __awaiter(this, void 0, void 0, function* () {
        const indexesPromise = indexes.listIndexes(projectId);
        const fieldOverridesPromise = indexes.listFieldOverrides(projectId);
        return Promise.all([indexesPromise, fieldOverridesPromise])
            .then((res) => {
            return indexes.makeIndexSpec(res[0], res[1]);
        })
            .catch((e) => {
            if (e.message.indexOf("is not a Cloud Firestore enabled project") >= 0) {
                return INDEXES_TEMPLATE;
            }
            throw new error_1.FirebaseError("Error fetching Firestore indexes", {
                original: e,
            });
        });
    });
}
