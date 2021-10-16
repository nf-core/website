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
const clc = require("cli-color");
const fsi = require("../firestore/indexes");
const logger = require("../logger");
const requirePermissions_1 = require("../requirePermissions");
module.exports = new command_1.Command("firestore:indexes")
    .description("List indexes in your project's Cloud Firestore database.")
    .option("--pretty", "Pretty print. When not specified the indexes are printed in the " +
    "JSON specification format.")
    .before(requirePermissions_1.requirePermissions, ["datastore.indexes.list"])
    .action((options) => __awaiter(void 0, void 0, void 0, function* () {
    const indexApi = new fsi.FirestoreIndexes();
    const indexes = yield indexApi.listIndexes(options.project);
    const fieldOverrides = yield indexApi.listFieldOverrides(options.project);
    const indexSpec = indexApi.makeIndexSpec(indexes, fieldOverrides);
    if (options.pretty) {
        logger.info(clc.bold.white("Compound Indexes"));
        indexApi.prettyPrintIndexes(indexes);
        if (fieldOverrides) {
            logger.info();
            logger.info(clc.bold.white("Field Overrides"));
            indexApi.printFieldOverrides(fieldOverrides);
        }
    }
    else {
        logger.info(JSON.stringify(indexSpec, undefined, 2));
    }
    return indexSpec;
}));
